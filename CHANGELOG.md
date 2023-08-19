
# Changelog

## Release 2.2 (work in progress)

### Package
- Moved to Python 3.8+ (Python modules and Python/C API).
- Updated the Python package definition, simplified the installation and adapted
  the instructions.
- Created an [Apptainer](https://apptainer.org/) container for a fully
  reproducible installation.
- Fixed numerous compilation problems.
- Updated, refactored and improved the test suite. Some of the small tutorials
  are now executed as part of the test suite.
- Removed several obsolete code sections.

### Functionality
- Default kinetic parameters in `EnergyOptions` now match the parameter preset
  `Options.JSDefault()`.
- A mistake was fixed in the numerical values for the parameter preset
  `DNA23Arrhenius()`.
- The `Options()` interface now uses explicit type casts to guard against type
  errors that would crash the Python/C API.
- The reliability of `MergeSim()` was improved by switching from the standard
  library module `multiprocessing` to the `multiprocess` package.
- `MergeSim()` now throws an error if the user-defined `OptionsFactory()` is not
  deterministic (up to Boltzmann sampling of initial conditions), or if it is
  not consistent with the `MergeSimSettings()`.
- C++ debug traces can now be toggled in the Python runtime via
  `Options.verbosity`.
- Transition types are now computed only if they are used in the kinetic model
  (`rate_method == 'Arrhenius'`).
- Migrated support from Nupack3 to Nupack4
- Transitioned energy models to `dna04`/`rna06`
- Created `utils.thermo` to wrap nupack util functions (multistrand doesn't account for coaxial stacking)
- Internal Boltzmann sampling uses updated nupack's `sample`
- Updated Jupyter Notebooks (`under_the_hood_notebooks`) to reflect changes in options as well as Python Migration
- Added constants (e.g. Boltzmann/C2K) to `options` or `utils.thermo` to minimize redefinitions for user
- Corrected undefined behavior related to loop indexing vs array indexing
- Corrected undefined behavior related to printing adjacent loops
- Added error message when `$NUPACKHOME` isn't set correctly
- Illegal initial structures that cause undefined behavior are now caught
- Shifted code from assertions in order to allow compilation with `NDEBUG`  
- New `BaseType` which replaces the previous `char` encoding method of passing sequence. Better readability
- Reduced number of compiler warnings
- Missing return statements fixed which caused segfaults
- Removed some unnecessary manual memory management between the Python and C++ layers.
- `MergeSim` concurrency uses `spawn` to properly isolate nupack utilities.
- Prevented repeated intermediate prints when using `MergeSim` due to `multiprocess` 
- Created `reuse_energymodel` option to minimize disk reads
- Fixed reference count bug which "mangles" complex objects. No longer crashes when accessing a complex post simulation
- General cleanup of print statements
- Introduced `const` expressions for better compile time optimization

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
