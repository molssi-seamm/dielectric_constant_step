# -*- coding: utf-8 -*-

"""This file contains metadata describing the results from DielectricConstant"""

from seamm_md_analysis.step_support import state_results

metadata = {}

"""Description of the computational models for DielectricConstant.

Hamiltonians, approximations, and basis set or parameterizations,
only if appropriate for this code. For example::

    metadata["computational models"] = {
        "Hartree-Fock": {
            "models": {
                "PM7": {
                    "parameterizations": {
                        "PM7": {
                            "elements": "1-60,62-83",
                            "periodic": True,
                            "reactions": True,
                            "optimization": True,
                            "code": "mopac",
                        },
                        "PM7-TS": {
                            "elements": "1-60,62-83",
                            "periodic": True,
                            "reactions": True,
                            "optimization": False,
                            "code": "mopac",
                        },
                    },
                },
            },
        },
    }
"""
# metadata["computational models"] = {
# }

"""Description of the DielectricConstant keywords.

(Only needed if this code uses keywords)

Fields
------
description : str
    A human readable description of the keyword.
takes values : int (optional)
    Number of values the keyword takes. If missing the keyword takes no values.
default : str (optional)
    The default value(s) if the keyword takes values.
format : str (optional)
    How the keyword is formatted in the MOPAC input.

For example::
    metadata["keywords"] = {
        "0SCF": {
            "description": "Read in data, then stop",
        },
        "ALT_A": {
            "description": "In PDB files with alternative atoms, select atoms A",
            "takes values": 1,
            "default": "A",
            "format": "{}={}",
        },
    }
"""
# metadata["keywords"] = {
# }

"""Properties that DielectricConstant produces.
`metadata["results"]` describes the results that this step can produce. It is a
dictionary where the keys are the internal names of the results within this step, and
the values are a dictionary describing the result. For example::

    metadata["results"] = {
        "total_energy": {
            "calculation": [
                "energy",
                "optimization",
            ],
            "description": "The total energy",
            "dimensionality": "scalar",
            "methods": [
                "ccsd",
                "ccsd(t)",
                "dft",
                "hf",
            ],
            "property": "total energy#Psi4#{model}",
            "type": "float",
            "units": "E_h",
        },
    }

Fields
______

calculation : [str]
    Optional metadata describing what subtype of the step produces this result.
    The subtypes are completely arbitrary, but often they are types of calculations
    which is why this is name `calculation`. To use this, the step or a substep
    define `self._calculation` as a value. That value is used to select only the
    results with that value in this field.

description : str
    A human-readable description of the result.

dimensionality : str
    The dimensions of the data. The value can be "scalar" or an array definition
    of the form "[dim1, dim2,...]". Symmetric tringular matrices are denoted
    "triangular[n,n]". The dimensions can be integers, other scalar
    results, or standard parameters such as `n_atoms`. For example, '[3]',
    [3, n_atoms], or "triangular[n_aos, n_aos]".

methods : str
    Optional metadata like the `calculation` data. `methods` provides a second
    level of filtering, often used for the Hamiltionian for *ab initio* calculations
    where some properties may or may not be calculated depending on the type of
    theory.

property : str
    An optional definition of the property for storing this result. Must be one of
    the standard properties defined either in SEAMM or in this steps property
    metadata in `data/properties.csv`.

type : str
    The type of the data: string, integer, or float.

units : str
    Optional units for the result. If present, the value should be in these units.
"""
metadata["results"] = {
    "epsilon": {
        "description": "The static dielectric constant",
        "dimensionality": "scalar",
        "property": "dielectric constant#Dielectric Constant#{model}",
        "type": "float",
        "units": "",
    },
    "epsilon,stderr": {
        "description": "The standard error of the static dielectric constant",
        "dimensionality": "scalar",
        "property": "dielectric constant,stderr#Dielectric Constant#{model}",
        "type": "float",
        "units": "",
    },
    "dipole correlation time": {
        "description": "The correlation time of the dipole of the cell",
        "dimensionality": "scalar",
        "type": "float",
        "units": "ps",
    },
    "dipole moment": {
        "description": "The mean dipole moment of a molecule",
        "dimensionality": "scalar",
        "type": "float",
        "units": "debye",
    },
    "Kirkwood factor": {
        "description": "The Kirkwood factor g_K = (⟨M²⟩ − ⟨M⟩²) / (N μ²)",
        "dimensionality": "scalar",
        "type": "float",
        "units": "",
    },
    "Kirkwood factor,stderr": {
        "description": "The standard error of the Kirkwood factor",
        "dimensionality": "scalar",
        "type": "float",
        "units": "",
    },
    **state_results,
    "nruns": {
        "description": "The number of MD runs analyzed",
        "dimensionality": "scalar",
        "type": "integer",
    },
    "n_frames": {
        "description": "The number of frames analyzed",
        "dimensionality": "scalar",
        "type": "integer",
    },
}
