# Multistrand

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


# Requirements

 -  c++11,  (clang, gcc v4.8.5+) 
 -  python2, 	 	2.7.12+
 -  nupack 3.0.4 only*
 -  make,			4.0+
 
The oldest NUPACK release available for download is 3.0.6 (www.nupack.org). Under 3.0.6, first step mode will not work, which relies on the sample subroutine. To make first step mode work in NUPACK 3.0.6, it is sufficient to copy sample.c from 3.0.4 into the 3.0.6 distribution and then to recompile.

Some users may need to install 'make' first. You can check that make is installed by simplying calling "make" in the terminal, which should return a message similar to the below. You can similarly run "python -V" and "gcc -v" to check if python and gcc are installed. Mac users may need to install xcode in order to proceed.
```sh
$ make: *** No targets specified and no makefile found.  Stop.
```


Tutorial files use the 'numpy', 'matplotlib' and 'scipy' python packages (you can install these using 'pip install numpy' and so on).

 
 # Installation
 
 - Clone the repository or download and unzip Multistrand into your workspace.
 - Make sure that in your enviroment (Eclipse, PyDev, bash, etc), NUPACKHOME points to the directory where nupack is installed. In the terminal, you can verify this by running 'echo $NUPACKHOME'.
 - Build multistrand by running 'make' in the Multistrand directory
 - The installation can be exported as a python library to the appropriate /site-packages/ by calling 'sudo make install'.
 
 Tutorials are found in /tutorials/. For docs, use 'make docs'



