"""Fail if a release tag and project metadata version disagree."""
from __future__ import annotations

from pathlib import Path
import sys
import tomllib


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_release_version.py vX.Y.Z")
    tag = sys.argv[1]
    version = tomllib.loads(Path("pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
    if tag != f"v{version}":
        raise SystemExit(f"tag {tag!r} does not match project version {version!r}")
    print(f"release version OK: {tag}")


if __name__ == "__main__":
    main()
