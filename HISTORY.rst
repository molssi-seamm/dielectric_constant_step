=======
History
=======
2026.10.17.1 -- Analyze MD run before the step; more results; labelled errors
    * The step analyzes MD run before it in the flowchart -- the runs of a Loop around
      the MD, as independent replicas ("the last MD run of each loop iteration", the
      default) -- or md.h5md files; it no longer runs its own MD in a sub-flowchart.
      Older flowcharts that stored its settings for that can still be read; the check
      before the job runs says what replaces them.
    * The check before the job runs also warns when the analysis is inside the loop
      holding its MD, when the loop changes the temperature or pressure, and when the
      thermostat does not suit the analysis, with fixes; the step warns before the MD
      when the system is small.
    * Errors are the 95% confidence interval, Student's t for the number of runs or
      blocks, and the output says so; stored errors are standard errors.
    * The state of the runs -- temperature, pressure, volume and density with their
      errors, and the equilibrium density where an NPT run before them gives the bulk
      modulus -- is reported and stored. The values of each run or block are in
      step.out.
    * New results: the mean dipole of a molecule, the Kirkwood factor g_K, and the mean
      dipole of the cell, which should be near zero; a statement that the formula
      assumes conducting (tin-foil) Ewald boundaries.
    * For an isotropic run ε is one value, with a check that the diagonal components
      agree; otherwise the components are given too.


2026.10.17 -- First release: the static dielectric constant
    * The static dielectric constant of liquids, and of the solvent of electrolytes,
      from the fluctuations of the dipole of the cell in MD with fixed charges, with
      errors from independent runs or blocks, the running estimate and the
      correlation time of the dipole.
    * Runs its own MD or analyzes the MD earlier in the flowchart, and checks before
      the job runs that the MD samples the positions often enough.
