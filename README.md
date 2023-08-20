# About

        ___  ___      _ _   _     _                       _ 
        |  \/  |     | | | (_)   | |                     | |
        | .  . |_   _| | |_ _ ___| |_ _ __ __ _ _ __   __| |
        | |\/| | | | | | __| / __| __| '__/ _` | '_ \ / _` |
        | |  | | |_| | | |_| \__ \ |_| | | (_| | | | | (_| |
        \_|  |_/\__,_|_|\__|_|___/\__|_|  \__,_|_| |_|\__,_|


Multistrand is a nucleic acids kinetic simulator, and is developed by the
Winfree group at the California Institute of Technology in Pasadena, California
(USA). Until 2013, development was lead by Joseph Schaeffer (now Autodesk).
The project is currently maintained by Jake Kaslewicz (Riedel group, University
of Minnesota, Minneapolis-Saint Paul) and Boyan Beronov (Condon group,
University of British Columbia, Vancouver).

Official website: www.multistrand.org

## Licence

    Multistrand nucleic acid kinetic simulator
    Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
    The Multistrand Team (help@multistrand.org)

Using this software is permitted for academic non-commercial purposes only. All
copyright is retained by Caltech.

**Disclaimer:** This software is provided "as is", without warrenty of any kind,
express or implied, including but not limited to the warrenties of
merchantability, fitness of a particular purpose and noninfringement. In no
event shall the authors or copyright holders be liable for any claim, damages or
other liability, whether in an action of contract, tort or otherwise, arising
from, out of or in connection with the software or the use or other dealings in
the software.

## Contributors

* Erik Winfree (winfree@caltech.edu)
* Chris Thachuk
* Frits Dannenberg (fdann@caltech.edu)
* Chris Berlind
* Joshua Loving
* Justin Bois
* Joseph Berleant
* Joseph Schaeffer
* Jake Kaslewicz (kasle001@umn.edu)
* Boyan Beronov (beronov@cs.ubc.ca)

Questions should be directed to help@multistrand.org.

# Usage

## Requirements

| Dependency                        | Notes              |
|-----------------------------------|--------------------|
| C++11                             | gcc 8+ or clang 8+ |
| Python                            | 3.8+               |
| [NUPACK](https://www.nupack.org/) | 4.0.1              |
 
The `numpy` and `scipy` Python packages are installed automatically as
dependencies, and `matplotlib` is added if the package extra `[tutorials]` is
specified (see `setup.cfg` for details).
 
## Installation

### Linux
 
 - `git clone` this repository into your workspace.
 - Run `pip install .` in the Multistrand directory.

### macOS

 - Install `xcode` commandline tools.
 - Install Python through `homebrew`.
 - Follow the Linux installation steps.
 - In `~/.bash_profile`, edit the `$PYTHONPATH` to include
   `/Library/Python/<python version>/site-packages`.
 
### Windows

 - Follow the instructions for installing the latest version of the [Microsoft
   C++ Build Tools](https://wiki.python.org/moin/WindowsCompilers).
 - Follow the Linux installation steps.
 
### [Apptainer](https://apptainer.org/) container

 - Install the Apptainer platform.
 - Build the container: `$> sudo apptainer build tools/multistrand.sif tools/multistrand.def`
 - Start the container: `$> apptainer shell --cleanenv --contain --pwd /dna/multistrand tools/multistrand.sif`

## Source tree

The Multistrand library is located under `src/`. `test/` is the test suite, and
`tools/` provides Apptainer container definitions and maintenance scripts.

## Testing

To execute the currently maintained portion of the test suite (including some of
the small tutorials):

 - Install the test dependencies: `$> pip install ".[testing]"`
 - Run: `$> pytest`

## Examples

Documentation can be found in `doc/` and `tutorials/`, and tutorial files are
organized as follows. The folder `under_the_hood/` contains in-depth tutorials,
and Jupyter versions are located in `under_the_hood_notebooks/`. The folder
`case_hybridization/` contains a case study into hybridization kinetics
(submission pending). Additional demo files are located in `misc/`.

As a very quick primer, we discuss two small scripts below.

### Hybridization trajectory

A quick test to see if Multistrand is working is to run the following script,
which simulates the hybridization of two complementary strands, ending the
simulation when the two strands either completely hybridize or seperate after an
initial collision:

```sh
$> python tutorials/misc/sample_trace.py
--------------------------------------------------------------------------------
GCGTTTCAC+GTGAAACGC
.(.......+.......).   t=0.000000 ms, dG=  0.16 kcal/mol
((.......+.......))   t=0.000060 ms, dG= -1.26 kcal/mol
GTGAAACGC+GCGTTTCAC
....(..((+))..)....   t=0.000526 ms, dG=  0.05 kcal/mol
....((.((+)).))....   t=0.000577 ms, dG= -1.46 kcal/mol
...(((.((+)).)))...   t=0.000608 ms, dG= -2.73 kcal/mol
...((((((+))))))...   t=0.000858 ms, dG= -7.90 kcal/mol
..(((((((+)))))))..   t=0.001025 ms, dG=-10.74 kcal/mol
.((((((((+)))))))).   t=0.001374 ms, dG= -9.79 kcal/mol
..(((((((+)))))))..   t=0.001421 ms, dG=-10.74 kcal/mol
.((((((((+)))))))).   t=0.002326 ms, dG= -9.79 kcal/mol
..(((((((+)))))))..   t=0.002601 ms, dG=-10.74 kcal/mol
..((((((.+.))))))..   t=0.002988 ms, dG= -9.33 kcal/mol
..(((((((+)))))))..   t=0.003122 ms, dG=-10.74 kcal/mol
..((((((.+.))))))..   t=0.003430 ms, dG= -9.33 kcal/mol
..(((((((+)))))))..   t=0.003570 ms, dG=-10.74 kcal/mol
..((((((.+.))))))..   t=0.003705 ms, dG= -9.33 kcal/mol
..(((((((+)))))))..   t=0.004507 ms, dG=-10.74 kcal/mol
..((((((.+.))))))..   t=0.006064 ms, dG= -9.33 kcal/mol
..(((((((+)))))))..   t=0.006210 ms, dG=-10.74 kcal/mol
..((((((.+.))))))..   t=0.006919 ms, dG= -9.33 kcal/mol
.(((((((.+.))))))).   t=0.007772 ms, dG= -8.37 kcal/mol
((((((((.+.))))))))   t=0.007780 ms, dG=-10.96 kcal/mol
(((((((((+)))))))))   t=0.008021 ms, dG=-12.38 kcal/mol
```

### Hybridization rates

The following script estimates the hybridization rate for a strand and its
complement. The computation relies on "first step mode" and will only work if
NUPACK is correctly installed. Alternatively, dissociation rates can be computed
by using `dissociation` as the first commandline argument. In that case, the
dissociation rate is computed indirectly by estimating the association rate, and
working out the dissociation rate from the partition function (e.g. `k_forward /
k_backward = exp(-dG / RT)` where `dG` is the partition function for the
complex, `R` is the gas constant and `T` is the temperature).

```sh
$> python tutorials/compute/rate.py hybridization AGCTGA -bootstrap
--------------------------------------------------------------------------------
2023-08-05 17:23:29   Starting Multistrand 2.2  (c) 2008-2023 Caltech

Running first step mode simulations for AGCTGA (with Boltzmann sampling)...

Start states:
Complex:
         Name: 'automatic0'
     Sequence: AGCTGA
    Structure: ......
      Strands: ['top']
    Boltzmann: True
  Supersample: 1

Complex:
         Name: 'automatic1'
     Sequence: TCAGCT
    Structure: ......
      Strands: ['top*']
    Boltzmann: True
  Supersample: 1

Stop conditions:
Stop Condition, tag:SUCCESS
  Sequence  0: AGCTGA+TCAGCT
  Structure 0: ((((((+))))))

Stop Condition, tag:FAILURE
  Sequence  0: AGCTGA
  Structure 0: ......

Using Results Type: FirstStepRate
Computing 1000 trials, using 10 threads ..
 .. and rolling 100 trajectories per thread until 500 successful trials occur.

Found 558 successful trials, terminating.
Done.  0.25239 seconds -- now processing results

The hybridization rate of AGCTGA and the reverse complement is 3.00e+06 /M /s
Bootstrapping FirstStepRate, using 1200 samples.    ..finished in 1.39 sec.

Estimated 95% confidence interval: [2.84e+06,3.16e+06]
Computing took 2.7906 s
```

## Log files

Multistrand automatically creates a logfile that contains some information on
the used model, e.g.:

```
$> cat multistrandRun.log
--------------------------------------------------------------------------------
Multistrand 2.2

sodium        :  1 M
magnesium     :  0 M
temperature   :  298.15 K
rate method   :  3  (1: Metropolis, 2: Kawasaki, 3: Arrhenius)
dangles       :  1  (0: none, 1: some, 2: all)
GT pairing    :  1  (0: disabled, 1: enabled)
concentration :  1 M
Nupack params :  /opt/bitnami/python/lib/python3.11/site-packages/nupack/parameters/dna04-nupack3.json

Kinetic parameters:
  type          End        Loop       Stack  StackStack     LoopEnd    StackEnd   StackLoop
  A      1.2965e+01  1.6424e+01  1.4184e+01  8.0457e+00  2.4224e+00  1.7524e+01  5.8106e+00
  E      3.4980e+00  4.4614e+00  5.2869e+00 -6.2712e-01  8.4934e-02  2.6559e+00 -1.1276e+00
  R      1.3584e+06  5.3080e+07  3.7097e+04  8.0873e+07  9.5396e+01  2.1242e+11  5.0138e+06

        dS_A        dH_A     biScale        kUni
 -0.0000e+00 -0.0000e+00  1.6006e-02  1.0000e+00

Rate matrix [ concentration . k_bi . k_uni(l,r) ]:
  2.1742e+04  1.3591e+05  3.5931e+03  1.6777e+05  1.8221e+02  8.5980e+06  4.1772e+04
  1.3591e+05  8.4962e+05  2.2461e+04  1.0487e+06  1.1390e+03  5.3747e+07  2.6112e+05
  3.5931e+03  2.2461e+04  5.9378e+02  2.7724e+04  3.0111e+01  1.4209e+06  6.9031e+03
  1.6777e+05  1.0487e+06  2.7724e+04  1.2945e+06  1.4059e+03  6.6342e+07  3.2231e+05
  1.8221e+02  1.1390e+03  3.0111e+01  1.4059e+03  1.5269e+00  7.2053e+04  3.5006e+02
  8.5980e+06  5.3747e+07  1.4209e+06  6.6342e+07  7.2053e+04  3.4000e+09  1.6518e+07
  4.1772e+04  2.6112e+05  6.9031e+03  3.2231e+05  3.5006e+02  1.6518e+07  8.0252e+04
```

# Frequently asked questions

## Capabilities
 
**Q:** Can I simulate leak reactions using Multistrand?

**A:** Yes. We have now added a preliminary tutorial, see `tutorials/leak_casestudy`.

## Troubleshooting

**Q:** How do I adjust the solvent salt concentrations?

**A:** Like so. (units are M = mol / litre)

```python
from multistrand.options import Options
o1 = Options()
o1.sodium = 0.05
o1.magnesium = 0.0125
```
