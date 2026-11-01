import os

from setuptools import Distribution, setup
from setuptools.command.bdist_wheel import bdist_wheel


class NativeBinaryDistribution(Distribution):
    """The package is platform-specific because it bundles a native library."""

    def has_ext_modules(self):
        return True


class PlatformIndependentAbiWheel(bdist_wheel):
    """Tag the ctypes-only wheel as Python/ABI independent, but platform bound."""

    def get_tag(self):
        _, _, platform_tag = super().get_tag()
        return "py3", "none", platform_tag


native_pattern = "native/*.dll" if os.name == "nt" else "native/*.so"
setup(
    distclass=NativeBinaryDistribution,
    cmdclass={"bdist_wheel": PlatformIndependentAbiWheel},
    package_data={"sctransform_cpu": [native_pattern]},
)
