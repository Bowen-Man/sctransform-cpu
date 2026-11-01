# Release process

Releases are built from immutable `vX.Y.Z` tags. The tag must equal the
`project.version` in `pyproject.toml`; CI rejects a mismatch.

## Required gates

1. Open a release branch and update `pyproject.toml`, `CHANGELOG.md` and the
   installation examples.
2. Pass Windows/Linux CI on Python 3.10, 3.11, 3.12 and 3.13.
3. Pass the Python 3.12 live R bit-identity gate in each platform's pinned
   reference environment. Functional CI on other Python versions is not a
   bit-identity claim.
4. Build and install each wheel in a new virtual environment. Confirm the
   native library loads from `site-packages` without an override variable.
5. Merge the reviewed pull request, create the matching annotated tag and push
   it. The release workflow builds two platform wheels plus one source archive,
   runs `twine check`, writes `SHA256SUMS` and publishes a GitHub Release.
6. Download the published artifacts and repeat the isolated install smoke test.

The release workflow packages the committed, oracle-validated native binaries.
The CI native-source job separately proves that Linux source still builds; it
does not silently replace the validated release binary.

## PyPI

PyPI publication is intentionally disabled until the Bowen-Man PyPI project
has an owner-controlled Trusted Publisher linked to this repository and a
protected `pypi` GitHub environment. Do not store a long-lived PyPI API token
in the repository. Once configured, add an OIDC publishing job that consumes
the already-tested release artifacts; do not rebuild them.
