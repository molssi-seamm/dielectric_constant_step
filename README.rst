=================================
SEAMM Dielectric Constant Plug-in
=================================

.. image:: https://img.shields.io/github/issues-pr-raw/molssi-seamm/dielectric_constant_step
   :target: https://github.com/molssi-seamm/dielectric_constant_step/pulls
   :alt: GitHub pull requests

.. image:: https://github.com/molssi-seamm/dielectric_constant_step/workflows/CI/badge.svg
   :target: https://github.com/molssi-seamm/dielectric_constant_step/actions
   :alt: Build Status

.. image:: https://codecov.io/gh/molssi-seamm/dielectric_constant_step/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/molssi-seamm/dielectric_constant_step
   :alt: Code Coverage

.. image:: https://github.com/molssi-seamm/dielectric_constant_step/workflows/CodeQL/badge.svg
   :target: https://github.com/molssi-seamm/dielectric_constant_step/security/code-scanning
   :alt: Code Quality

.. image:: https://github.com/molssi-seamm/dielectric_constant_step/workflows/Release/badge.svg
   :target: https://molssi-seamm.github.io/dielectric_constant_step/index.html
   :alt: Documentation Status

.. image:: https://img.shields.io/pypi/v/dielectric_constant_step.svg
   :target: https://pypi.python.org/pypi/dielectric_constant_step
   :alt: PyPi VERSION

A SEAMM plug-in for the static dielectric constant of liquids and electrolytes from MD.

* Free software: BSD-3-Clause
* Documentation: https://molssi-seamm.github.io/dielectric_constant_step/index.html
* Code: https://github.com/molssi-seamm/dielectric_constant_step

Features
--------

* The static dielectric constant of liquids, and of the solvent of electrolytes,
  from the fluctuations of the dipole of the cell, with errors from runs or blocks,
  the running estimate and the dipole's correlation time.
* Analyzes the MD run before it in the flowchart -- the runs of a Loop around the MD,
  as independent runs -- or md.h5md files of earlier runs.

Acknowledgements
----------------

This package was created with the `molssi-seamm/cookiecutter-seamm-plugin`_ tool, which
is based on the excellent Cookiecutter_.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`molssi-seamm/cookiecutter-seamm-plugin`: https://github.com/molssi-seamm/cookiecutter-seamm-plugin

Developed by the Molecular Sciences Software Institute (MolSSI_),
which receives funding from the `National Science Foundation`_ under
award CHE-2136142.

.. _MolSSI: https://molssi.org
.. _`National Science Foundation`: https://www.nsf.gov
