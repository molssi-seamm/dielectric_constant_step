.. _user-guide:

**********
User Guide
**********

The Dielectric Constant step calculates the static dielectric constant ε of a liquid,
or of the solvent of an electrolyte, from the fluctuations of the dipole of the cell
in equilibrium molecular dynamics:

.. math::

   \varepsilon = 1 + \frac{\langle \mathbf{M}^2 \rangle - \langle \mathbf{M} \rangle^2}
                        {3 V \varepsilon_0 k_B T}

with M = Σ qᵢ rᵢ the dipole of the molecules, V the volume and T the temperature.
The formula holds for Ewald or PPPM electrostatics with conducting boundaries, as
LAMMPS uses; the step warns when the MD used anything else.

Which molecules
===============

"Dipole of" chooses the molecules whose dipole is followed:

* **automatic** (the default): every molecule of a neutral liquid; for an
  electrolyte, the **neutral molecules** only. The ions carry a current, not a bounded
  dipole, so including them makes M grow without limit. The result is the dielectric
  constant of the solvent in the presence of the ions -- the quantity usually reported
  for electrolyte solutions (the dielectric decrement). It leaves out the ion-solvent
  cross term, which a rigorous separation of the ionic current would include.
* **the neutral molecules** or **every molecule**, explicitly.

Molecules are sorted by their total charge: neutral molecules, ions (monatomic like
Na⁺ or molecular like BF₄⁻), and any molecule whose charge is not an integer, which
is a problem with the charges and is left out with a note.

The charges are the forcefield's. A machine-learned forcefield has no fixed charges,
so the step stops with a message for such MD.

Where the MD comes from
=======================

Like the other analysis steps, the step runs its own MD in its sub-flowchart a number
of times, or analyzes the MD that ran earlier in the flowchart, or the ``md.h5md``
files of earlier runs. The MD must sample the atoms' positions ("Sample the atomic
positions" in the LAMMPS NVE, NVT and NPT sub-steps) at least every 0.5 ps ("Longest
interval between frames"); the step turns this on for its own MD, and checking the
flowchart reports it for MD earlier in the flowchart. The positions are unwrapped, so
molecules stay whole.

Convergence and errors
======================

ε converges slowly: the dipole of the cell decorrelates over picoseconds (the step
reports this correlation time) and its fluctuations need nanoseconds. The error is the
standard error over independent runs, or over blocks of a single run ("Blocks for
errors from one run"). The step notes when a run is shorter than 50 correlation times,
and graphs the running estimate of ε of each run, with the autocorrelation function
of the dipole (``dielectric_constant.graph``).

Results
=======

ε and its standard error (which can be stored as properties of the configuration),
the correlation time of the dipole, the mean temperature and volume, and the numbers
of runs and frames.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
