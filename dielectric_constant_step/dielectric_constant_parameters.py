# -*- coding: utf-8 -*-
"""
Control parameters for the Dielectric Constant step in a SEAMM flowchart
"""

import logging
import seamm
from seamm_md_analysis.step_support import source_parameters

logger = logging.getLogger(__name__)


class DielectricConstantParameters(seamm.Parameters):
    """
    The control parameters for Dielectric Constant.

    The MD to analyze is run earlier in the flowchart -- in a Loop for independent
    runs -- or in other jobs (see seamm_md_analysis.step_support.source_parameters).
    """

    # Dropped with the own-MD route (phase 8): ignored in older flowcharts
    retired = ("nruns", "errors")

    parameters = {
        **source_parameters("the dielectric constant"),
        "frame interval": {
            "default": 0.5,
            "kind": "float",
            "default_units": "ps",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Longest interval between frames:",
            "help_text": "The positions must be sampled at least this often.",
        },
        "dipole of": {
            "default": "automatic",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ("automatic", "the neutral molecules", "every molecule"),
            "format_string": "s",
            "description": "Dipole of:",
            "help_text": (
                "The molecules whose dipole fluctuations give the dielectric "
                "constant. 'automatic' takes every molecule of a neutral liquid and "
                "the neutral molecules of an electrolyte: the ions carry a current, "
                "not a bounded dipole."
            ),
        },
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("DielectricConstantParameters.__init__")

        super().__init__(
            defaults={**DielectricConstantParameters.parameters, **defaults}, data=data
        )
