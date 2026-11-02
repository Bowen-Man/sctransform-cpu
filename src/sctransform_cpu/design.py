"""R-compatible SCT model-matrix construction for supported v2 designs."""
from __future__ import annotations

from dataclasses import dataclass
import re
import numpy as np

from .port_native import platform_log10_native
from .model import r45_offset_from_library


@dataclass(frozen=True)
class SCTDesign:
    fit: np.ndarray
    full: np.ndarray
    coefficient_names: tuple[str, ...]
    offset: np.ndarray
    model_str: str
    attributes: dict[str, np.ndarray]
    nonreg_full: np.ndarray | None
    nonreg_coefficient_names: tuple[str, ...]
    model_str_nonreg: str
    batch_levels: tuple[object, ...] = ()


def _cell_attributes(counts, library, supplied, required):
    attributes = {} if supplied is None else {
        str(name): np.asarray(values) for name, values in supplied.items()
    }
    cells = counts.shape[1]
    for name, values in attributes.items():
        if values.ndim != 1 or len(values) != cells:
            raise ValueError(f'cell attribute "{name}" must have one value per cell')
    if "log_umi" in attributes:
        expected = platform_log10_native(np.asarray(library, dtype=np.float64))
        try:
            provided = np.asarray(attributes["log_umi"], dtype=np.float64)
        except (TypeError, ValueError) as error:
            raise ValueError('cell attribute "log_umi" must be numeric') from error
        if not np.array_equal(provided.view(np.uint64), expected.view(np.uint64)):
            raise ValueError('supplied "log_umi" must equal log10 of library size')

    missing = set(required) - set(attributes)
    known = {
        "umi", "gene", "log_umi", "log_gene", "umi_per_gene",
        "log_umi_per_gene",
    }
    unknown = missing - known
    if unknown:
        raise ValueError(f"unknown cell attributes: {', '.join(sorted(unknown))}")
    if missing & {"umi", "log_umi", "umi_per_gene", "log_umi_per_gene"}:
        attributes["umi"] = np.asarray(library, dtype=np.float64)
        attributes["log_umi"] = platform_log10_native(attributes["umi"])
    if missing & {"gene", "log_gene", "umi_per_gene", "log_umi_per_gene"}:
        attributes["gene"] = np.asarray(counts.getnnz(axis=0), dtype=np.float64)
        attributes["log_gene"] = platform_log10_native(attributes["gene"])
    if missing & {"umi_per_gene", "log_umi_per_gene"}:
        attributes["umi_per_gene"] = attributes["umi"] / attributes["gene"]
        attributes["log_umi_per_gene"] = platform_log10_native(
            attributes["umi_per_gene"])
    return attributes


def build_design(counts, library, latent_variables, batch_variable,
                 latent_variables_nonreg, cell_attributes, fit_cells) -> SCTDesign:
    """Build the model matrices used by ``vst(..., method='glmGamPoi_offset')``."""
    latent = tuple(latent_variables)
    nonreg = tuple(latent_variables_nonreg)
    if (not latent or not all(isinstance(name, str) for name in latent)
            or len(set(latent)) != len(latent) or "log_umi" not in latent):
        raise ValueError("latent_variables must contain distinct names including log_umi")
    if (not all(isinstance(name, str) for name in nonreg)
            or len(set(nonreg)) != len(nonreg)
            or set(latent) & set(nonreg)):
        raise ValueError(
            "latent_variables_nonreg must contain distinct names not in latent_variables")
    if batch_variable is not None and not isinstance(batch_variable, str):
        raise ValueError("batch_variable must be a string name")
    if any(name in {"(Intercept)", "Intercept"} for name in latent):
        raise ValueError("reserved intercept names cannot be latent variables")
    if batch_variable in {"(Intercept)", "Intercept"}:
        raise ValueError("reserved intercept names cannot be batch variables")
    identifiers = latent + nonreg + (() if batch_variable is None else (batch_variable,))
    if any(not re.fullmatch(r"[A-Za-z.][A-Za-z0-9._]*", name)
           for name in identifiers):
        raise ValueError("model variables must use simple R-compatible names")
    required = latent + nonreg + (() if batch_variable is None else (batch_variable,))
    attributes = _cell_attributes(counts, library, cell_attributes, required)
    fit_cells = np.asarray(fit_cells, dtype=int)

    numeric = {}
    for name in latent + nonreg:
        try:
            values = np.asarray(attributes[name], dtype=np.float64)
        except (TypeError, ValueError) as error:
            raise ValueError(f'latent variable "{name}" must be numeric') from error
        if not np.isfinite(values).all():
            raise ValueError(f'cell attribute "{name}" contains a non-finite value')
        numeric[name] = values

    if batch_variable is None:
        fitted_names = tuple(name for name in latent if name != "log_umi")
        fit = np.column_stack(
            (np.ones(len(fit_cells)), *(numeric[name][fit_cells] for name in fitted_names)))
        names = ("(Intercept)", "log_umi", *fitted_names)
        full = np.column_stack(
            (np.ones(counts.shape[1]), numeric["log_umi"],
             *(numeric[name] for name in fitted_names)))
        offset = r45_offset_from_library(
            np.asarray(library, dtype=np.float64)[fit_cells])
        return SCTDesign(
            np.asarray(fit, dtype=np.float64), np.asarray(full, dtype=np.float64),
            names, offset, f"y ~ {' + '.join(latent)}", attributes,
            (np.column_stack((np.ones(counts.shape[1]),
                              *(numeric[name] for name in nonreg)))
             if nonreg else None),
            (("(Intercept)", *nonreg) if nonreg else ()),
            (f"y ~ {' + '.join(nonreg)}" if nonreg else ""))

    batch = np.asarray(attributes[batch_variable])
    if batch.dtype.kind not in "US":
        raise ValueError("batch variable must contain string labels")
    if any(not str(value).isascii() or not str(value) for value in batch):
        raise ValueError("batch labels must be nonempty ASCII strings")
    levels = tuple(np.unique(batch).tolist())
    fitted_levels = set(batch[fit_cells].tolist())
    if fitted_levels != set(levels):
        raise ValueError("fit cells must retain every batch level")

    columns = []
    names = []
    for level in levels:
        indicator = (batch == level).astype(np.float64)
        columns.append(indicator)
        names.append(f"{batch_variable}{level}")
    for variable in latent:
        for level in levels:
            columns.append(numeric[variable] * (batch == level))
            names.append(f"{variable}:{batch_variable}{level}")
    full = np.column_stack(columns)
    fit = full[fit_cells]
    model = f"y ~ ({' + '.join(latent)}) : {batch_variable} + {batch_variable} + 0"
    nonreg_columns = []
    nonreg_names = []
    if nonreg:
        for level in levels:
            nonreg_columns.append((batch == level).astype(np.float64))
            nonreg_names.append(f"{batch_variable}{level}")
        for variable in nonreg:
            for level in levels:
                nonreg_columns.append(numeric[variable] * (batch == level))
                nonreg_names.append(f"{variable}:{batch_variable}{level}")
    nonreg_model = (
        f"y ~ ({' + '.join(nonreg)}) : {batch_variable} + {batch_variable} + 0"
        if nonreg else "")
    return SCTDesign(
        np.asarray(fit, dtype=np.float64), np.asarray(full, dtype=np.float64),
        tuple(names), np.zeros(len(fit_cells)), model, attributes,
        (np.asarray(np.column_stack(nonreg_columns), dtype=np.float64)
         if nonreg else None),
        tuple(nonreg_names), nonreg_model, levels)


def correction_design(design: SCTDesign, latent_variables,
                      batch_variable, scale_factor) -> np.ndarray:
    """Return the official ``correct()`` target model matrix."""
    latent = tuple(latent_variables)
    numeric = {}
    for name in latent:
        values = np.asarray(design.attributes[name], dtype=np.float64)
        numeric[name] = np.full(len(values), np.median(values))
    if scale_factor is not None:
        if not np.isfinite(scale_factor) or scale_factor <= 0:
            raise ValueError("scale_factor must be finite and positive")
        numeric["log_umi"][:] = platform_log10_native(
            np.asarray([scale_factor], dtype=np.float64))[0]

    if batch_variable is None:
        fitted_names = tuple(name for name in latent if name != "log_umi")
        return np.asarray(np.column_stack(
            (np.ones(len(numeric["log_umi"])), numeric["log_umi"],
             *(numeric[name] for name in fitted_names))), dtype=np.float64)

    batch = np.asarray(design.attributes[batch_variable])
    columns = [(batch == level).astype(np.float64)
               for level in design.batch_levels]
    for variable in latent:
        for level in design.batch_levels:
            columns.append(numeric[variable] * (batch == level))
    return np.asarray(np.column_stack(columns), dtype=np.float64)
