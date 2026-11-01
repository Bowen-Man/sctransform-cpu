"""Install one wheel in a temporary environment and load its backend."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import venv

import numpy
import scipy


def main() -> None:
    directory = Path(sys.argv[1] if len(sys.argv) > 1 else "dist").resolve()
    wheels = list(directory.glob("*.whl"))
    if len(wheels) != 1:
        raise SystemExit(f"expected one wheel in {directory}, found {len(wheels)}")
    with tempfile.TemporaryDirectory(prefix="sct-wheel-") as temporary:
        root = Path(temporary)
        venv.EnvBuilder(with_pip=True).create(root)
        python = root / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
        subprocess.run(
            [
                str(python), "-m", "pip", "install", "--no-deps",
                "--disable-pip-version-check", str(wheels[0]),
            ],
            check=True,
        )
        dependency_paths = {
            str(Path(module.__file__).resolve().parent.parent)
            for module in (numpy, scipy)
        }
        environment = os.environ.copy()
        environment["PYTHONPATH"] = os.pathsep.join(sorted(dependency_paths))
        subprocess.run(
            [
                str(python),
                "-c",
                "import site, sys; sys.path.insert(0, site.getsitepackages()[0]); "
                "from sctransform_cpu import __version__; "
                "from sctransform_cpu.port_native import native_library_path; "
                "print(__version__, native_library_path())",
            ],
            check=True,
            env=environment,
        )


if __name__ == "__main__":
    main()
