from setuptools import Distribution, setup


class WindowsBinaryDistribution(Distribution):
    """The package is platform-specific because it bundles a ctypes DLL."""

    def has_ext_modules(self):
        return True


setup(distclass=WindowsBinaryDistribution)
