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

The step analyzes MD that ran before it in the flowchart, or the ``md.h5md`` files of
earlier runs. For independent runs, put the MD in a Loop -- in parallel, with each
run starting from the same structure ("Discard the structure" in the MD sub-steps)
and new velocities -- and analyze "the last MD run of each loop iteration", the
default. A typical flowchart:

.. code-block:: text

    Model Chemistry (FF@...), Packmol, minimization
    NPT equilibration; Set Cell to its mean volume; NVT equilibration
    Loop over the replicas (parallel)
        LAMMPS: velocities, NVT decorrelation, NVT production
    Dielectric Constant

The MD must sample the atoms' positions ("Sample the atomic positions" in the LAMMPS
NVE, NVT and NPT sub-steps) at least every 0.5 ps ("Longest interval between
frames"). The positions are unwrapped, so molecules stay whole. Checking the flowchart
before it runs reports what is missing, with a fix that turns the sampling on, and
also warns when the analysis is inside the loop holding its MD, when the loop changes
the temperature or pressure (the runs are then not replicas of one state), or when
the thermostat does not give canonical fluctuations (Berendsen, plain rescaling).
Before the MD runs, the step warns if the system is small for the analysis (fewer than
about 200 molecules).

Each run must be long compared with the correlation time of the dipole, about 11 ps
for water: 300 ps or more.

Convergence and errors
======================

ε converges slowly: the dipole of the cell decorrelates over picoseconds (the step
reports this correlation time) and its fluctuations need nanoseconds. The error is
from independent runs, or from blocks of a single run ("Blocks for errors from one
run"); the table gives the 95% confidence interval (Student's t times the standard
error) and says so. The step notes when a run is shorter than 50 correlation times,
and graphs the running estimate of ε of each run, with the autocorrelation function
of the dipole (``dielectric_constant.graph``).

For an isotropic run -- a cubic cell with nothing directional in it -- ε is one value,
a third of the trace of the dipole fluctuations, and a line checks that the diagonal
components agree within their errors. Otherwise the components are given too. The
values of each run or block are in the detailed output.

Results
=======

* ε with its error, and the correlation time of the dipole;
* the mean dipole of a molecule μ, the Kirkwood factor
  g_K = (⟨M²⟩ − ⟨M⟩²) / (N μ²), and the mean dipole of the cell ⟨M⟩, which should
  be close to zero;
* the state of the runs: temperature, pressure, volume and density with their errors,
  and for runs at constant volume the equilibrium density at the target pressure when
  an NPT run before them gives the bulk modulus;
* the numbers of runs and frames.

ε, its error and the state can be stored as properties or in variables and tables.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
