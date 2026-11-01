"""Install the single wheel in an isolated environment and load its backend."""
from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile
import venv


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
            [str(python), "-m", "pip", "install", "--disable-pip-version-check", str(wheels[0])],
            check=True,
        )
        subprocess.run(
            [
                str(python),
                "-c",
                "from sctransform_cpu import __version__; "
                "from sctransform_cpu.port_native import native_library_path; "
                "print(__version__, native_library_path())",
            ],
            check=True,
        )


if __name__ == "__main__":
    main()
