"""R 4.5 compatible sampling used by the targeted ``vst`` v2 path.

The implementation is intentionally narrow: R's default Mersenne-Twister
``RNGkind`` and rejection ``sample.kind`` are reproduced for uniform and
weighted sampling without replacement.  This keeps fit selections explicit
and testable without making R a runtime dependency.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy as np


_UINT32_MASK = (1 << 32) - 1
_MATRIX_A = 0x9908B0DF
_UPPER_MASK = 0x80000000
_LOWER_MASK = 0x7FFFFFFF


class RMersenneTwister:
    """R 4.5's default uniform RNG after ``set.seed(seed)``.

    This follows ``src/main/RNG.c`` rather than NumPy's MT19937 seeding,
    because equal integer seeds do not otherwise create equal streams.
    """

    def __init__(self, seed: int):
        value = int(seed) & _UINT32_MASK
        for _ in range(50):
            value = (69069 * value + 1) & _UINT32_MASK
        state = []
        for _ in range(625):
            value = (69069 * value + 1) & _UINT32_MASK
            state.append(value)
        self._state = state[1:]
        self._index = 624

    def _next_uint32(self) -> int:
        if self._index >= 624:
            state = self._state
            for index in range(624 - 397):
                value = ((state[index] & _UPPER_MASK)
                         | (state[index + 1] & _LOWER_MASK))
                state[index] = (state[index + 397] ^ (value >> 1)
                                ^ (_MATRIX_A if value & 1 else 0))
            for index in range(624 - 397, 623):
                value = ((state[index] & _UPPER_MASK)
                         | (state[index + 1] & _LOWER_MASK))
                state[index] = (state[index - (624 - 397)] ^ (value >> 1)
                                ^ (_MATRIX_A if value & 1 else 0))
            value = ((state[623] & _UPPER_MASK)
                     | (state[0] & _LOWER_MASK))
            state[623] = (state[396] ^ (value >> 1)
                          ^ (_MATRIX_A if value & 1 else 0))
            self._index = 0
        value = self._state[self._index]
        self._index += 1
        value ^= value >> 11
        value ^= (value << 7) & 0x9D2C5680
        value ^= (value << 15) & 0xEFC60000
        value ^= value >> 18
        return value & _UINT32_MASK

    def uniform(self) -> float:
        value = self._next_uint32() * 2.3283064365386963e-10
        if value <= 0.0:
            return 0.5 * 2.328306437080797e-10
        if 1.0 - value <= 0.0:
            return 1.0 - 0.5 * 2.328306437080797e-10
        return value

    def _bits(self, bits: int) -> int:
        value = 0
        for _ in range(0, bits + 1, 16):
            value = 65536 * value + int(math.floor(self.uniform() * 65536.0))
        return value & ((1 << bits) - 1)

    def index(self, size: int) -> int:
        if size <= 0:
            return 0
        bits = int(math.ceil(math.log2(size)))
        while True:
            value = self._bits(bits)
            if value < size:
                return value

    def sample_uniform(self, population_size: int, sample_size: int) -> np.ndarray:
        if not 0 <= sample_size <= population_size:
            raise ValueError("sample_size must be between zero and population_size")
        population = np.arange(population_size, dtype=np.int64)
        result = np.empty(sample_size, dtype=np.int64)
        remaining = population_size
        for offset in range(sample_size):
            selected = self.index(remaining)
            result[offset] = population[selected]
            remaining -= 1
            population[selected] = population[remaining]
        return result

    def sample_weighted(self, probabilities, sample_size: int) -> np.ndarray:
        probability = np.asarray(probabilities, dtype=np.float64).copy()
        if probability.ndim != 1 or not np.isfinite(probability).all():
            raise ValueError("probabilities must be a finite vector")
        if np.any(probability < 0.0) or np.count_nonzero(probability) < sample_size:
            raise ValueError("too few positive probabilities")
        probability /= probability.sum()
        identity = np.arange(len(probability), dtype=np.int64)
        _r_revsort(probability, identity)

        result = np.empty(sample_size, dtype=np.int64)
        total_mass = 1.0
        remaining = len(probability) - 1
        for offset in range(sample_size):
            target = total_mass * self.uniform()
            mass = 0.0
            selected = remaining
            for candidate in range(remaining):
                mass += probability[candidate]
                if target <= mass:
                    selected = candidate
                    break
            result[offset] = identity[selected]
            total_mass -= probability[selected]
            probability[selected:remaining] = probability[selected + 1:remaining + 1]
            identity[selected:remaining] = identity[selected + 1:remaining + 1]
            remaining -= 1
        return result


def _r_revsort(values: np.ndarray, identity: np.ndarray) -> None:
    """In-place translation of R's descending ``revsort`` heapsort."""
    size = len(values)
    if size <= 1:
        return
    # The R implementation is one-indexed. Padding keeps its tie behaviour.
    val = np.empty(size + 1, dtype=np.float64)
    idx = np.empty(size + 1, dtype=np.int64)
    val[1:], idx[1:] = values, identity
    left = (size >> 1) + 1
    right = size
    while True:
        if left > 1:
            left -= 1
            held_value, held_index = val[left], idx[left]
        else:
            held_value, held_index = val[right], idx[right]
            val[right], idx[right] = val[1], idx[1]
            right -= 1
            if right == 1:
                val[1], idx[1] = held_value, held_index
                break
        parent = left
        child = left << 1
        while child <= right:
            if child < right and val[child] > val[child + 1]:
                child += 1
            if held_value > val[child]:
                val[parent], idx[parent] = val[child], idx[child]
                parent = child
                child += parent
            else:
                child = right + 1
        val[parent], idx[parent] = held_value, held_index
    values[:], identity[:] = val[1:], idx[1:]


def r_nrd_density_weights(values) -> np.ndarray:
    """Return the inverse density weights used by ``vst`` gene sampling.

    The calculation mirrors R 4.5 ``density(..., bw='nrd')`` followed by
    ``approx`` at the original values.  The FFT result is used only to make a
    discrete sampling decision; the selected titles are the compatibility
    contract and are covered by a live-R gate.
    """
    data = np.asarray(values, dtype=np.float64)
    if data.ndim != 1 or len(data) < 2 or not np.isfinite(data).all():
        raise ValueError("density sampling needs at least two finite values")
    scale = min(float(np.std(data, ddof=1)),
                float((np.quantile(data, 0.75) - np.quantile(data, 0.25)) / 1.34))
    bandwidth = 1.06 * scale * len(data) ** -0.2
    if not math.isfinite(bandwidth) or bandwidth <= 0.0:
        raise ValueError("R's nrd bandwidth is not positive")

    size = 512
    lower_output = float(np.min(data) - 3.0 * bandwidth)
    upper_output = float(np.max(data) + 3.0 * bandwidth)
    lower = lower_output - 4.0 * bandwidth
    upper = upper_output + 4.0 * bandwidth
    delta = (upper - lower) / (size - 1)
    binned = np.zeros(2 * size, dtype=np.float64)
    weight = 1.0 / len(data)
    for value in data:
        position = (float(value) - lower) / delta
        index = math.floor(position)
        fraction = position - index
        if 0 <= index <= size - 2:
            binned[index] += (1.0 - fraction) * weight
            binned[index + 1] += fraction * weight
        elif index == -1:
            binned[0] += fraction * weight
        elif index == size - 1:
            binned[index] += (1.0 - fraction) * weight

    kernel_coordinates = np.linspace(
        0.0, ((2 * size - 1) / (size - 1)) * (upper - lower), 2 * size)
    kernel_coordinates[size + 1:] = -kernel_coordinates[size - 1:0:-1]
    kernel = (np.exp(-0.5 * (kernel_coordinates / bandwidth) ** 2)
              / (math.sqrt(2.0 * math.pi) * bandwidth))
    density = np.maximum(
        0.0,
        np.fft.ifft(np.fft.fft(binned) * np.conj(np.fft.fft(kernel))).real[:size],
    )
    density_grid = np.linspace(lower, upper, size)
    output_grid = np.linspace(lower_output, upper_output, size)
    output_density = np.interp(output_grid, density_grid, density)
    sampled_density = np.interp(data, output_grid, output_density)
    return 1.0 / (sampled_density + np.finfo(np.float64).eps)


@dataclass(frozen=True)
class SelectionManifest:
    """Ordered titles selected for the first fitting step."""

    seed: int
    cells: tuple[str, ...]
    genes: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return {"seed": self.seed, "cells": list(self.cells), "genes": list(self.genes)}

    def write(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2) + "\n", encoding="utf-8")

    @classmethod
    def read(cls, path: str | Path) -> "SelectionManifest":
        value = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls(int(value["seed"]), tuple(value["cells"]), tuple(value["genes"]))


def make_selection_manifest(
    counts,
    *,
    gene_names,
    cell_names,
    seed: int,
    n_cells: int = 2000,
    n_genes: int = 2000,
    min_cells: int = 5,
) -> SelectionManifest:
    """Select ordered step-1 titles without running model fitting.

    The function implements the same filtering, RNG consumption, and weighted
    sampling order as the targeted v2 pipeline.  Saving its result provides a
    language-neutral manifest for Python/R comparisons and repeated runs.
    """
    from scipy import sparse
    from .r_compat_stats import row_statistics_r45

    matrix = sparse.csr_matrix(counts, dtype=np.float64, copy=True)
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    matrix.sort_indices()
    genes = np.asarray(gene_names, dtype=str)
    cells = np.asarray(cell_names, dtype=str)
    if matrix.shape != (len(genes), len(cells)):
        raise ValueError("matrix shape does not match gene and cell titles")
    if len(set(genes)) != len(genes) or len(set(cells)) != len(cells):
        raise ValueError("gene and cell titles must be unique")
    if min_cells < 1 or n_cells < 1 or n_genes < 1:
        raise ValueError("selection sizes and min_cells must be positive")
    if (not np.isfinite(matrix.data).all() or np.any(matrix.data < 0.0)
            or np.any(matrix.data != np.floor(matrix.data))):
        raise ValueError("counts must be finite non-negative integers")

    keep = np.asarray(matrix.getnnz(axis=1)).ravel() >= min_cells
    matrix, genes = matrix[keep].tocsr(), genes[keep]
    means, variances, log_gmean = row_statistics_r45(matrix)
    rng = RMersenneTwister(seed)
    if n_cells < matrix.shape[1]:
        selected_cells = rng.sample_uniform(matrix.shape[1], n_cells)
    else:
        selected_cells = np.arange(matrix.shape[1])
    detected = np.asarray(matrix[:, selected_cells].getnnz(axis=1)).ravel()
    eligible = np.flatnonzero((detected >= min_cells) & (variances > means))
    if n_genes < len(eligible):
        chosen = rng.sample_weighted(r_nrd_density_weights(log_gmean[eligible]), n_genes)
        eligible = eligible[chosen]
    if len(eligible) < 2:
        raise ValueError("fewer than two overdispersed genes are available")
    return SelectionManifest(
        int(seed),
        tuple(cells[selected_cells]),
        tuple(genes[eligible]),
    )
