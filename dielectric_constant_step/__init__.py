# -*- coding: utf-8 -*-

"""
dielectric_constant_step
A SEAMM plug-in for Dielectric Constant
"""

# Bring up the classes so that they appear to be directly in
# the dielectric_constant_step package.

from .dielectric_constant import DielectricConstant  # noqa: F401, E501
from .dielectric_constant_parameters import DielectricConstantParameters  # noqa: F401
from .dielectric_constant_step import DielectricConstantStep  # noqa: F401, E501
from .tk_dielectric_constant import TkDielectricConstant  # noqa: F401, E501

from .metadata import metadata  # noqa: F401

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
