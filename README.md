# Multistrand #

Multistrand is a tool for kinetic simulations of nucleic acids, and is developed by the Winfree group at Caltech in Pasadena, CA, USA. Until 2013, development was lead by Joseph Schaeffer (now Autodesk).

Contributors:

Erik Winfree			winfree@caltech.edu
Chris Thachuk
Frits Dannenberg    	fdann@caltech.edu
Chris Berlind
Joshua Loving
Justin Bois
Joseph Berleant
Joseph Schaeffer

General questions should be directed to help@multistrand.org. Also see www.multistrand.org


Frits Dannenberg, May 26rd, 2017


## Licence ##

Multistrand, kinetic simulator for nucleic acids.
Copyright 2017, California Institute of Technology. All rights reserved.

Using this software is permitted for academic non-commercial purposes only. All copyright is retained by Caltech. 

Disclaimer: This software is provided "as is", without warrenty of any kind, express or implied, including
but not limited to the warrenties of merchantability, fitness of a particular purpose and 
noninfringement. In no event shall the authors or copyright holders be liable for any claim,
damages or otehr liability, whether in an action of contract, tort or otherwise, arising from,
out of or in connection with the software or the use or other dealings in the software.


## Requirements ##

 -  c++11,  (clang, gcc v4.8.5+) 
 -  python2, 	 	2.7.12+
 -  nupack 3.0.4 only*
 -  make,			4.0+
 
The earliest NUPACK release available for download is 3.0.6 (www.nupack.org). Under 3.0.6, first step mode will not work, which relies on the sample subroutine. To make first step mode work in NUPACK 3.0.6, it is sufficient to copy sample.c from 3.0.4 into the 3.0.6 distribution and then to recompile. First step mode will not work with 3.1.0 and 3.2.0.

Some users may need to install 'make' first. You can check that make is installed by simplying calling "make" in the terminal, which should return a message similar to the below. You can similarly run "python -V" and "gcc -v" to check if python and gcc are installed. Mac users may need to install xcode in order to proceed.
```sh
$ make: *** No targets specified and no makefile found.  Stop.
```

Tutorial files use the 'numpy', 'matplotlib' and 'scipy' python packages (you can install these using 'pip install numpy' and so on).

 
## Installation ##
 
 - Clone the repository or download and unzip Multistrand into your workspace.
 - Make sure that in your enviroment (Eclipse, PyDev, bash, etc), NUPACKHOME points to the directory where nupack is installed. In the terminal, you can verify this by running 'echo $NUPACKHOME'.
 - Build multistrand by running 'make' in the Multistrand directory
 - The installation can be exported as a python library to the appropriate /site-packages/ by calling 'sudo make install'.
 
## Package tree ##


.
/energymodel
/include
/interface
/loop
/multistrand
/nupack
/obj
/state
/system
/test
/tutorials
Makefile
README.md
setup.py

Source dirs are /energymodel, /include, /interface/, /loop, /nupack, /state, /system, /test. Build dirs are /multistrand and /obj and buildfiles are Makefile and setup.py.

Tutorial files are organized as follows. Folder /under_the_hood/ contains in depth tutorials. Jupyter versions of the /under_the_hood/ are found in /under_the_hood_notebooks/. /case_hybridization/ contains an upcomming case study into hybridization kinetics (submission pending). Additional demo files are located in /misc/.

## Using Multistrand ##

Several tutorial files on using multistrand are found in /tutorials/under_the_hood. As a very quick primer, we discuss two small scripts below.  Multistrand also comes with documentation (run 'make docs').


### Hybridization trajectory ###

A quick test to see if Multistrand is working is by running 'python tutorials/misc/sample_trace.py'. This example does not use First Step mode, and should work would with all recent versions of NUPACK (3.0.6, 3.1.0, 3.2.0). This script simulates the hybridization of two complementary strands and ends the simulation when the two strands either completely hybridize or seperate after an initial collision:  

```sh
['/home/iris/workspace/caltech/Multistrand/tutorials/misc/sample_trace.py']
GTGAAACGC GCGTTTCAC
......... ......... t=0.0000 ms, dG=0.00 kcal/mol  
GTGAAACGC+GCGTTTCAC
(........+........) t=0.0406 ms, dG=0.22 kcal/mol  
GCGTTTCAC+GTGAAACGC
......(.(+).)...... t=0.0465 ms, dG=1.18 kcal/mol  
......(((+)))...... t=0.0495 ms, dG=-3.22 kcal/mol  
......((.+.))...... t=0.0588 ms, dG=-0.63 kcal/mol  
......(((+)))...... t=0.0725 ms, dG=-3.22 kcal/mol  
.....((((+))))..... t=0.1901 ms, dG=-3.21 kcal/mol  
....(((((+))))).... t=0.2732 ms, dG=-4.48 kcal/mol  
..(.(((((+))))).).. t=0.2921 ms, dG=-4.03 kcal/mol  
..(((((((+))))))).. t=0.2936 ms, dG=-7.97 kcal/mol  
.((((((((+)))))))). t=0.2971 ms, dG=-10.96 kcal/mol  
(((((((((+))))))))) t=0.3337 ms, dG=-12.38 kcal/mol  
```

### Hybridization Rates ###

The following script computes a hybridization rate for a strand and its compement. The computation relies on first step mode and will only work with the appropriate version of NUPACK.

```sh
[iris@dhcp-135-182 Multistrand]$ python tutorials/misc/computeAnnealRate.py 'AGCTGA' -bootstrap
2017-05-30 18:05:34  Starting Multistrand 2.1      (c) 2008-2017 Caltech      
Running first step mode simulations for AGCTGA (with Boltzmann sampling)...
Computing 800 trials, using 2 threads .. 
Done.  0.25586 seconds 

nForward = 419 
nReverse = 381 
k1 = 7332500.0 /M/s   
k_eff = 7332500.0 /M/s   

The hybridization rate of AGCTGA and the reverse complement is 7.33e+06 /M /s
Estimated 95% confidence interval: [6.86e+06,7.79e+06] 
```


### Log files ###

Multistrand automatically creates a logfile ("multistrandRun.log") that contains some information on the used model, like so:
```
Sodium      :  0.5 M 
Magnesium   :  0 M 
Temperature :  298.15 K
Rate method :  2           (1: Kawasaki, 2: Metropolis)
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
o1.sodium = 0.05		# units: mol / litre
o1.magnesium = 0.0125	# units: mol / litre
```


