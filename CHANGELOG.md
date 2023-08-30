
# Changelog

## Release 2.2 (work in progress)

### Package
- Migrated from Python 2.7 to Python 3.9+ (Python modules & Python/C API).
- Migrated from NUPACK 3 to NUPACK 4, while still using thermodynamic models
  compatible with NUPACK 3.
- Updated the Python package definition, simplified the installation and adapted
  the instructions.
- Dropped the dependence on the environment variable `$NUPACKHOME`.
- Created an [Apptainer](https://apptainer.org/) container for fully
  reproducible development and deployment. For reproducible debugging, there is
  also a container variant including debug builds of Python and Multistrand and
  a pre-configured GDB.
- Updated the Jupyter notebooks (`tutorials/under_the_hood_notebooks`).
- Updated, refactored and improved the test suite. Some of the small tutorials
  are now executed as part of the test suite.

### Functionality
- Updated the NUPACK thermodynamic parameters to `dna04`/`rna99`.
- Introduced the `utils.thermo` module, which wraps the updated NUPACK utility
  functions (e.g., Boltzmann sampling, disabling coaxial stacking in
  thermodynamic model, adjusting concentration units in ensemble free energy).
- Defined repeatedly used physical constants (e.g., Boltzmann, Celsius to
  Kelvin) in the `options` and `utils.thermo` modules.
- Updated the default kinetic parameters in `EnergyOptions` (C++) to match the
  parameter preset `Options.JSDefault()` (Python).
- Tied the formerly static `EnergyModel` (C++) instance as a dynamic attribute
  to an `Options` (Python) instance. Henceforth, different `SimSystem` (Python)
  instances can reuse the *same* `EnergyModel` instance *sequentially* (i.e., to
  avoid re-parsing the same NUPACK parameter file), as well as use *different*
  `EnergyModel` instances *concurrently* (e.g., for varying environment
  conditions or kinetic parameters).
- Enabled toggling simulator debug traces from the Python runtime via
  `Options.verbosity`, i.e., without recompiling the C++ extension.
- Improved the reliability of `MergeSim` (Python) by switching from the standard
  library module `multiprocessing` to the `multiprocess` package, and by using
  `spawn` for concurrency in order to properly isolate NUPACK utilities.
- Added an exception in `MergeSim` when the user-defined `OptionsFactory` is not
  deterministic (up to Boltzmann sampling of initial states), or if it is not
  consistent with the `MergeSimSettings`.

### Internals
- Consolidated the C++ extension module based on the new treatment of `Options`.
- Relocated code which was previously inside assertions, in order to enable
  compilation with `NDEBUG`.
- Resolved a number of compiler warnings.
- Introduced `const` expressions for better compile time optimization.
- Introduced `BaseType`, which replaces the previous `char` encoding of primary
  structure, in order to improve maintainability.
- Stopped computing transition types when they are not used in the kinetic model
  (`rate_method != 'Arrhenius'`).
- Prevented repeated intermediate prints when using `multiprocess` in
  `MergeSim`.
- Cleaned up print statements.
- Removed several obsolete code sections.

### Bug fixes
#### Python
- Corrected swapped dimensions in the parameter preset
  `Options.DNA23Arrhenius()`.
- Added explicit type casts in the `Options` interface, in order to guard
  against type errors that would crash the Python/C API.

#### C++
- Fixed some missing return statements which caused segmentation faults.
- Fixed a reference counting bug which previously mangled `Complex` objects,
  leading to a crash when accessing a `Complex` post simulation.
- Removed some unnecessary manual memory management between the Python and C++
  layers.
- Corrected undefined behavior related to loop indexing vs. array indexing.
- Corrected undefined behavior related to printing adjacent loops.
- Caught illegal initial structures that would cause undefined behavior.


## Release 2.1 (May 2018)

- Moved to c++11 standard.
- Reworked c++ internals: new files, classes, structs, enums, consts and various
  tostrings added. For example: `struct BaseCount` instead of `int[]` for
  storing exposed bases in a complex or loop. Settings imported from Python are
  now cached.
- Many files that were no longer relevant or not working are now removed.
  Multistrand no longer supports the ViennaRNA thermodynamic model.
- Tests of long-run equilibrium and equivalence of the energy model (w.r.t.
  NUPACK) are now working.
- `SimSystem` object now supports `.initialInfo()`, which prints the internal
  representation of the initial state (also see `tutorials/misc/inspection.py`).
- Fixed a bug (pointer comparison) that prevented deterministic execution.
- Fixed a bug that required the user to rename parameter files when simulating
  RNA.
- Multistrand now ships with convenience classes for multithreading
  (`interface/concurrent.py` -- various demos, for example
  `tutorials/misc/computeAnnealRate.py`).
- Multistrand will now complain when the user does not explicitly set the uni-
  and bi-molecular rate constants, before setting them automatically (as
  before). Default parameterization is provided in `multistrand.utils`.
- Multistrand now prints a warning when no initial moves are available in first
  step mode.
- Multistrand now prints a warning when the simulation time is exceeded, but no
  stopping condition is met (provided they are set).
- Multistrand now creates a logfile (`multistrandRun.log`) containing some of
  the simulation details. This file is overwritten each time a simulation is
  started.
- Now supports buffer conditions (see f.a.q.).
- Updated documentation.
- Updated the installation guidelines.
- Removed unused code from the installation files.
- Added case study files for hybridization (`tutorials/hybridization_casestudy`).
- Added a case study on leak rates (`tutorials/leak_casestudy`).
- Added a function that generates a unique (hashable) value for each visited
  state (see `tutorials/misc/uniqueID.py`).
- Added a commandline utility (and tutorial file) to compute hybridization rates
  (`tutorials/misc/computeAnnealRate.py`).
- Convenience functions for standard experimental setups are available in
  `multistrand.experiment` (for example
  `multistrand.experiment.hybridization(multistrand.Options(), string:sequence)`
  -- also see `tutorials/misc/computeAnnealRate.py`).

## (Aug 2017)

- Expanding the leak case study (Mirnank Sharma), added FirstStepLeak object to
  handle leak simulations in concurrent.py (FD and MS).
- A macOS bugfix.
- Added error messages for object initialization (Mirnank Sharma).
- Added compute folder for (commandline) utility in `tutorials/`.
- Removed unused code in Make files and initialization routine.
- Added files for simulating Machinek-Turberfield mismatch paper.

## (Jul 2017)

- Added support for nupack 3.2.1


# Known issues

- Because multithreading is handled outside the C++ core, error messages print
  once for each thread.
