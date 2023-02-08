from __future__ import print_function
# hairpin_energies.py
#
# This example walks you through how to create a single-stranded complex (i.e. a strand sequence and structure), and how to find its energy. 
#
# Invoke as "python hairpin_energies.py" to see a plot.
# Invoke as "python -i hairpin_energies.py" to see a plot and then drop into the python interpreter to investigate further.
# Invoke as "ipython --pylab -i hairpin_energies.py" to see a plot and then drop into the ipython interpreter to investigate further.

if False:  # only needed if you're having trouble with your Multistrand installation
    import multistrand_setup

try:
    from multistrand.objects import *
    from multistrand.options import Options
    from multistrand.system import energy, initialize_energy_model

except ImportError:
    print("Could not import Multistrand.")
    raise

#############


o = Options(temperature=25,dangles="Some")    # prepares for simulation.
initialize_energy_model(o)                    # necessary if you want to use energy() without running a simulation first. 
# see more about the energy model usage and initialization in threewaybm_trajectories.py

# More meaningful names for argument values to the energy() function call, below.
Loop_Energy = 0    # requesting no dG_assoc or dG_volume terms to be added.  So only loop energies remain.
Volume_Energy = 1  # requesting dG_volume but not dG_assoc terms to be added.  No clear interpretation for this.
Complex_Energy = 2 # requesting dG_assoc but not dG_volume terms to be added.  This is the NUPACK complex microstate energy, sans symmetry terms.
Tube_Energy = 3    # requesting both dG_assoc and dG_volume terms to be added.  Summed over complexes, this is the system state energy.

# Sequence is from Schaeffer's PhD thesis, chapter 7, figure 7.1

# Just for illustration, create a hairping strand with just the outermost 4 base pairs of the stem formed:
c = Complex( strands=[Strand(name="hairpin", sequence="GTTCGGGCAAAAGCCCGAAC")], structure= '((((' + 12*'.' + '))))' )
energy( [c], o, Complex_Energy)  # should be -1.1449...
# Note that energy() takes a *list* of complexes, and returns a tuple of energies.  Don't give it just a complex as input, else all hell may break loose.

# Try 'help(energy)' for a little more information.
# Similarly, 'help(initialize_energy_model)'  or  'help(Options)'  or  'help(Complex)'  or  'help(Strand)' .
# But beware that the help docs are not always complete or up-to-date, sorry.

# Using this sequence, find the energy for a particular secondar structure conformation.
def print_hp(s):
    e = energy( [Complex( strands=[Strand(name="hairpin", sequence="GTTCGGGCAAAAGCCCGAAC")], structure=s)], o, Complex_Energy)[0]  
    print(s + '  (%5.2f)' % e)
    return e

# Manually define a set of secondary structures for our hairpin, closing from the outside.
path1=[0]*9
path1[0] = print_hp('....................')
path1[1] = print_hp('(..................)')
path1[2] = print_hp('((................))')
path1[3] = print_hp('(((..............)))')
path1[4] = print_hp('((((............))))')
path1[5] = print_hp('(((((..........)))))')
path1[6] = print_hp('((((((........))))))')
path1[7] = print_hp('(((((((......)))))))')
path1[8] = print_hp('((((((((....))))))))')

print()

# Algorithmically define a set of secondary structures for our hairpin, closing from the inside.
path2=[0]*9
for i in range(9):
    path2[i] = print_hp('.'*(8-i) + '('*i + '....' + ')'*i + '.'*(8-i))

steps    = range(9)

# If matplotlib is available, we can plot the above computed values.
def myplot():
    import numpy as np
    import matplotlib
    import matplotlib.pylab as plt
    
    plt.figure(1)
    plt.plot(steps,path1,'go', label='Outside bases first')
    plt.hold(True)
    plt.plot(steps,path1,'g-', label='_nolabel_')
    plt.plot(steps,path2,'rx', label='Inside bases first')
    plt.plot(steps,path2,'r-', label='_nolabel_')
    plt.hold(False)
    plt.title("Energy landscape for two hairpin folding pathways")
    plt.xlabel("Base pairs formed",fontsize='larger')
    plt.ylabel("Microstate Energy (kcal/mol)",fontsize='larger')
    plt.yticks(fontsize='larger',va='bottom')
    plt.xticks(fontsize='larger')
    plt.legend()
    plt.show()

# Execute the following only if run from the command line, e.g. as "python -i hairpin_energies.py" or "ipython --pylab -i hairpin_energies.py"
if __name__ == '__main__':
    myplot()
else:
    print("Try:\nhairpin_energies.myplot()\n")
