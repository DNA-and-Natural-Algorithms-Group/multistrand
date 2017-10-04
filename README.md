# Multistrand #
        ___  ___      _ _   _     _                       _ 
        |  \/  |     | | | (_)   | |                     | |
        | .  . |_   _| | |_ _ ___| |_ _ __ __ _ _ __   __| |
        | |\/| | | | | | __| / __| __| '__/ _` | '_ \ / _` |
        | |  | | |_| | | |_| \__ \ |_| | | (_| | | | | (_| |
        \_|  |_/\__,_|_|\__|_|___/\__|_|  \__,_|_| |_|\__,_|  

Multistrand is a nucleic acids kinetic simulator, and is developed by the Winfree group at the California Institute of Technology in Pasadena, California (USA). Until 2013, development was lead by Joseph Schaeffer (now Autodesk).


Contributors:

* Erik Winfree			      winfree@caltech.edu
* Chris Thachuk
* Frits Dannenberg    	fdann@caltech.edu
* Chris Berlind
* Joshua Loving
* Justin Bois
* Joseph Berleant
* Joseph Schaeffer

Questions should be directed to help@multistrand.org. Also see www.multistrand.org


Frits Dannenberg, May 26rd, 2017


## Licence ##

Multistrand, kinetic simulator for nucleic acids.
Copyright 2017, California Institute of Technology. All rights reserved.

Using this software is permitted for academic non-commercial purposes only. All copyright is retained by Caltech. 

Disclaimer: This software is provided "as is", without warrenty of any kind, express or implied, including
but not limited to the warrenties of merchantability, fitness of a particular purpose and 
noninfringement. In no event shall the authors or copyright holders be liable for any claim,
damages or other liability, whether in an action of contract, tort or otherwise, arising from,
out of or in connection with the software or the use or other dealings in the software.


## Requirements ##

| Dependency | Notes               | 
| ---------- | ------              |
| c++11      | gcc 4.8.5+ or clang/llvm 8.0.0   |  
|  python2   |  2.7.12+       	   | 
| NUPACK - www.nupack.org    |  **3.2.1**  | 
|  make      |  4.0+         | 
 
As of mid-2017, NUPACK 3.2.1 was released, which restores the sample functionality that Multistrand depends on. Users no longer need to patch nupack in order to use first step mode. 

Some users may need to install 'make' first. You can check that make is installed by simplying calling "make" in the terminal, which should return a message similar to the below. You can similarly run "python -V",  "gcc -v" and "clang -V" to check which python, gcc and clang are installed. Mac users may need to install xcode in order to proceed. 
```sh
$ make: *** No targets specified and no makefile found.  Stop.
```

Tutorial files use the 'numpy', 'matplotlib' and 'scipy' python packages. You can install these using 'pip install numpy', and so on.

 
## Installation ##
 
 For linux:
 
 - Clone the repository into your workspace, using 'git clone *repo-url*'.
 - In your enviroment (eclipse, bash, etc), set NUPACKHOME to point the directory where NUPACK is installed. 
 - Build multistrand by running 'make' in the Multistrand directory.
 - Multistrand can be exported as a python library by calling 'sudo make install'.

In Fedora, add 'export NUPACKHOME=/path/to/nupack3.2.1' to ~./bashrc to make the export permanent.
To verify that NUPACKHOME is set correctly in bash, run 'echo $NUPACKHOME':
```sh
$ echo $NUPACKHOME
/home/user/path/to/nupack3.2.1
```
For OS X, the following steps enabled successful installation on a 2017 macbook pro. 

 - Install xcode commandline tools
 - Install brew (homebrew) and install python through homebrew (command: 'brew install python').
 - The linux installation steps should now work. 
 - In the bash profile file ~/.bash_profile, edit the python path to include /Library/Python/2.7/site-packages. For example, the line could read: export PYTHONPATH=/Library/Python/2.7/site-packages:$PYTHONPATH. 

Update 2017-10-1: We've noticed that linking against anaconda sometimes gives errors during installation on OSX. The issue seems to be that the linked gcc, 4.2.1, does to not support certain c++11 functions. OS X users who rely on brew-installed python do not have this problem, because the gcc call is piped to llvm, which then works. 

To resolve this issue, users should add the lines

```
os.environ["CC"] = "clang"
os.environ["CXX"] = "clang"
```


after the line that reads "config_vars = distutils.sysconfig.get_config_vars()" in the setup.py file.




## Package tree ##

Source dirs are /energymodel, /include, /interface/, /loop, /nupack, /state, /system, /test. Build dirs are /multistrand and /obj and buildfiles are Makefile and setup.py.

Tutorial files are organized as follows. Folder under_the_hood contains in depth tutorials. Jupyter versions are located in under_the_hood_notebooks. The folder case_hybridization contains a case study into hybridization kinetics (submission pending). Additional demo files are located in /misc/.

## Using Multistrand ##

In-depth tutorial files on multistrand are found in /tutorials/under_the_hood. The other files in /tutorials/ are typically set up as case studies. As a very quick primer, we discuss two small scripts below. 

### Hybridization trajectory ###

A quick test to see if Multistrand is working is by running 'python tutorials/misc/sample_trace.py'. This example does not use First Step mode, and should work would with all recent versions of NUPACK (3.0.6, 3.1.0, 3.2.1). This script simulates the hybridization of two complementary strands and ends the simulation when the two strands either completely hybridize or seperate after an initial collision:  

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

### Hybridization rates ###

The following script computes a hybridization rate for a strand and its compement. The computation relies on first step mode and will only work if NUPACK is correctly installed. Alternatively, dissociation rates can be computed by using 'dissociation' as the first commandline argument. In that case, the dissociation rate is computed indirectly by estimating the association rate, and working out the dissociation rate from the partition function (e.g. k_forward / k_backward = exp( -dG / RT) where dG is the partition function for the complex, R is the gas constant and T is temperature). 
Note: the released code does not yet include the calibrated parameters from the DNA23 conference. If you would like to use the calibrated model, please contact me directly (FD, 2017-10-03).

```sh
[iris@dhcp-135-182 Multistrand]$ python tutorials/compute/rate.py hybridization 'AGCTGA' -bootstrap
2017-10-03 13:31:25  Starting Multistrand 2.1      (c) 2008-2017 Caltech      
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




### Log files ###

Multistrand automatically creates a logfile ("multistrandRun.log") that contains some information on the used model, like so:
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
 
 ### Frequently Aksed Questions ###
 
Q: Can I simulate leak reactions using Multistrand?

A: Yes. We have now added a preliminary tutorial, see /tutorials/leak_casestudy.

Q: How do I adjust the solvent salt concentrations?

A: Like so. (units are M = mol / litre) 

```python
from multistrand.options import Options
o1 = Options()
o1.sodium = 0.05        # units: mol / litre
o1.magnesium = 0.0125   # units: mol / litre
```

Q: When I run tutorials/misc/computeAnnealRate.py, nForward and nReverse are both zero and the program does not terminate.

A:  This occurs when NUPACK returns void output for 'sample'. If NUPACK is installed correctly, then running 

``` bash
./nupack/bin/sample -multi -T 25 -material dna -samples 100	
./nupack/build/bin/sample -multi -T 25 -material dna -samples 100      (v.3.2.1 and up)
```
and supplying the arguments 'test' '1' 'AGTGTGCGTAGA' '1' will result in a list of 100 non-trivial dot-paren secondary structures in the 'test.sample' file. 
NUPACK 3.0.4 only: if you have patched NUPACK, be sure to rebuild (make clean; make) the package. Unpatched NUPACK 3.0.4 will return a void output. Some nupack releases (3.1, 3.2) do not have the sample function included.

Q: When I run a Multistrand script, I get ImportError: No module named Multistrand. 

A: Please check that PYTHONPATH includes a link to your Multistrand executables, or run 'sudo make install' in the Multistrand directory to install Multistrand python module on your system.
