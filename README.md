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

Official website: www.multistrand.org

## Licence

Multistrand, kinetic simulator for nucleic acids.
Copyright 2023, California Institute of Technology. All rights reserved.

Using this software is permitted for academic non-commercial purposes only. All copyright is retained by Caltech. 

Disclaimer: This software is provided "as is", without warrenty of any kind, express or implied, including
but not limited to the warrenties of merchantability, fitness of a particular purpose and 
noninfringement. In no event shall the authors or copyright holders be liable for any claim,
damages or other liability, whether in an action of contract, tort or otherwise, arising from,
out of or in connection with the software or the use or other dealings in the software.

## Contributors

* Erik Winfree (winfree@caltech.edu)
* Chris Thachuk
* Frits Dannenberg (fdann@caltech.edu)
* Chris Berlind
* Joshua Loving
* Justin Bois
* Joseph Berleant
* Joseph Schaeffer
* Boyan Beronov

Questions should be directed to help@multistrand.org.

# Usage

## Requirements

| Dependency | Notes  |
| ---------- | ------ |
| C++11      | gcc 7.4 or clang 8.0.0  |
| Python     | 3.8+   |
| [NUPACK](https://www.nupack.org/) | 3.2.2 |
 
The `numpy` and `scipy` Python packages are installed automatically as
dependencies, and `matplotlib` is added if the installation target `tutorials`
is specified (see `setup.cfg` for details).
 
## Installation

### Linux
 
 - `git clone` this repository into your workspace.
 - Set the environment variable `$NUPACKHOME` to point to the NUPACK
   installation directory.
 - Run `pip install .` in the Multistrand directory.

### macOS

 - Install `xcode` commandline tools.
 - Install Python through `homebrew`.
 - Follow the Linux installation steps.
 - In `~/.bash_profile`, edit the `$PYTHONPATH` to include
   `/Library/Python/3.8/site-packages`.
 
### Windows

 - Follow the instructions for installing the latest version of the [`Microsoft
   C++ Build Tools](https://wiki.python.org/moin/WindowsCompilers).
 - Follow the Linux installation steps.
 
### [Apptainer](https://apptainer.org/) container

 - Build: `$> sudo apptainer build tools/multistrand.sif tools/multistrand.def`
 - Start: `$> apptainer shell --cleanenv --contain --pwd /dna/multistrand tools/multistrand.sif`

## Source tree

The Multistrand library is located under `src/`, whereas `nupack/` contains the
Nupack wrapper. `test/` is the test suite, and `tools/` provides Apptainer
container definitions.

## Examples

Documentation can be found in `doc/` and `tutorials/`, and tutorial files are
organized as follows. Folder `under_the_hood/` contains in-depth tutorials, and
Jupyter versions are located in `under_the_hood_notebooks/`. The folder
`case_hybridization/` contains a case study into hybridization kinetics
(submission pending). Additional demo files are located in `misc/`.

As a very quick primer, we discuss two small scripts below.

### Hybridization trajectory

A quick test to see if Multistrand is working is to run `python
tutorials/misc/sample_trace.py`. This script simulates the hybridization of two
complementary strands and ends the simulation when the two strands either
completely hybridize or seperate after an initial collision:
```sh
GTGAAACGC GCGTTTCAC
......... .........   t=0.000000 ms,  dG=0.00 kcal/mol  
GTGAAACGC+GCGTTTCAC
(........+........)   t=0.000041 ms,  dG=0.22 kcal/mol  
GCGTTTCAC+GTGAAACGC
......(.(+).)......   t=0.000046 ms,  dG=1.18 kcal/mol  
......(((+)))......   t=0.000049 ms,  dG=-3.22 kcal/mol  
......((.+.))......   t=0.000059 ms,  dG=-0.63 kcal/mol  
......(((+)))......   t=0.000073 ms,  dG=-3.22 kcal/mol  
.....((((+)))).....   t=0.000190 ms,  dG=-3.21 kcal/mol  
....(((((+)))))....   t=0.000273 ms,  dG=-4.48 kcal/mol  
..(.(((((+))))).)..   t=0.000292 ms,  dG=-4.03 kcal/mol  
..(((((((+)))))))..   t=0.000294 ms,  dG=-7.97 kcal/mol  
.((((((((+)))))))).   t=0.000297 ms,  dG=-10.96 kcal/mol  
(((((((((+)))))))))   t=0.000334 ms,  dG=-12.38 kcal/mol  
```

### Hybridization rates

The following script estimates the hybridization rate for a strand and its
complement. The computation relies on "first step mode" and will only work if
NUPACK is correctly installed. Alternatively, dissociation rates can be computed
by using `dissociation` as the first commandline argument. In that case, the
dissociation rate is computed indirectly by estimating the association rate, and
working out the dissociation rate from the partition function (e.g. `k_forward /
k_backward = exp(-dG / RT)` where `dG` is the partition function for the
complex, `R` is the gas constant and `T` is temperature).

```sh
$> python tutorials/compute/rate.py hybridization 'AGCTGA' -bootstrap
2017-10-03 13:31:25  Starting Multistrand 2.1      (c) 2008-2023 Caltech
Running first step mode simulations for AGCTGA (with Boltzmann sampling)...

Computing 1200 trials, using 6 threads .. 
 .. and rolling 200 trajectories per thread until 500 successful trials occur. 
Found 957 successful trials, terminating.

nForward = 957 
nReverse = 1443 
 
k1       = 5.58e+06  /M /s  

Done.  0.91273 seconds 

The hybridization rate of AGCTGA and the reverse complement is 5.58e+06 /M /s
Estimated 95% confidence interval: [5.32e+06,5.86e+06]
```

## Log files

Multistrand automatically creates a logfile (`./multistrandRun.log`) that
contains some information on the used model, e.g.:

```
Sodium      :  0.5 M 
Magnesium   :  0 M 
Temperature :  298.15 K
Rate method :  1           (1: Metropolis, 2:  Kawasaki)
dangles     :  1           (0: none, 1: some, 2: all)
substrate   :  1           (1: DNA)
GT pairing  :  0           (0: disabled)

 biScale     kUni    
 1.4e+06     5e+06
```
 

# Frequently Aksed Questions

## Capabilities
 
**Q:** Can I simulate leak reactions using Multistrand?

**A:** Yes. We have now added a preliminary tutorial, see `tutorials/leak_casestudy`.

## Troubleshooting

**Q:** When I try to run any multistrand script, the console returns `segfault 11`.

**A:** Please make sure `$NUPACKHOME` is set. When you run `echo $NUPACKHOME` in
bash, it should return the directory of your nupack installation.

**Q:** How do I adjust the solvent salt concentrations?

**A:** Like so. (units are M = mol / litre)

```python
from multistrand.options import Options
o1 = Options()
o1.sodium = 0.05        # units: mol / litre
o1.magnesium = 0.0125   # units: mol / litre
```

**Q:** When I run `tutorials/misc/computeAnnealRate.py`, `nForward` and `nReverse`
are both zero and the program does not terminate.

**A:** This occurs when NUPACK returns void output for 'sample'. If NUPACK is installed correctly, then running 

``` bash
./nupack/bin/sample -multi -T 25 -material dna -samples 100
./nupack/build/bin/sample -multi -T 25 -material dna -samples 100      (v.3.2.1 and up)
```
and supplying the arguments `'test' '1' 'AGTGTGCGTAGA' '1'` will result in a
list of 100 non-trivial dot-paren secondary structures in the `test.sample`
file.

**Q:** When I run a Multistrand script, I get `ImportError: No module named Multistrand`.

**A:** Please check that `$PYTHONPATH` includes a link to your Multistrand
executables, or run `sudo pip install .` in the Multistrand directory for a
system-wide installation of the Multistrand Python module.

**Q:** When I try to run `tutorials/misc/sample_trace.py`, I get the error:

```
['tutorials/misc/sample_trace.py']
terminate called after throwing an instance of 'std::logic_error'
  what():  basic_string::_S_construct null not valid
Aborted
```

**A:** You need to set the `$NUPACKHOME` enviroment variable. You can see the
current value by using `echo $NUPACKHOME`.

**Q:** On macOS, I get the following error when I try to install:

```
In file included from system/utility.cc:16:
In file included from ./include/move.h:29:
In file included from ./include/moveutil.h:14:
./include/sequtil.h:66:22: error: no matching constructor for initialization of
      'vector<int>'
        vector<int> count = { 0, 0, 0, 0, 0 }; // use baseType as access
                            ^~~~~~~~~~~~~~~~~
```

**A:** We've noticed that linking against Anaconda sometimes gives errors during
installation on macOS. The issue seems to be that the linked gcc does to not
support certain c++11 functions. macOS users who rely on brew-installed Python
do not have this problem, because the gcc call is piped to llvm, which then
works. To resolve this issue, users should prepend the enviroment variables

```
CC=clang CXX=clang
```
before the installation command.
