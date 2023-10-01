# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
This example walks you through how to create a single-stranded complex (i.e. a
strand sequence and structure), and how to find its energy.

Invoke as "python hairpin_energies.py" to see a plot.

Invoke as "python -i hairpin_energies.py" to see a plot and then drop into the
python interpreter to investigate further.

Invoke as "ipython --pylab -i hairpin_energies.py" to see a plot and then drop
into the ipython interpreter to investigate further.
"""

from multistrand.objects import *
from multistrand.options import Options, EnergyType
from multistrand.system import calculate_energy


opt = Options(temperature=25, dangles="Some")   # prepares for simulation.
opt.DNA23Metropolis()
# see more about the energy model usage and initialization in
# threewaybm_trajectories.py

# Sequence is from Schaeffer's PhD thesis, chapter 7, figure 7.1

# Just for illustration, create a hairping strand with just the outermost 4 base
# pairs of the stem formed:
seq = "GTTCGGGCAAAAGCCCGAAC"
N = len(seq)
c = Complex(strands=[Strand(name="hairpin", sequence=seq)],
            structure= '((((' + 12*'.' + '))))')
calculate_energy([c], opt, EnergyType.complex)  # should be -1.1449...
# Note that energy() takes a *list* of complexes, and returns a tuple of
# energies. Don't give it just a complex as input, else all hell may break
# loose.

# For more information, try:
#   'help(multistrand.objects)'
#   'help(multistrand.options)'
#   'help(multistrand.system)'
# But beware that the help docs are not always complete or up-to-date, sorry.

# Using this sequence, find the energy for a particular secondar structure
# conformation.
def print_hp(s):
    e = calculate_energy(
        [Complex(strands=[Strand(name="hairpin", sequence=seq)], structure=s)],
        opt, EnergyType.complex)[0]
    print(f'{s}  ({e:6.2f})')
    return e

# Define a set of secondary structures for our hairpin, closing from the
# outside.
steps = range(9)
path1 = [print_hp(f"{'('*i}{'.'*(N-2*i)}{')'*i}")
         for i in steps]
print()

# Define a set of secondary structures for our hairpin, closing from the inside.
path2 = [print_hp(f"{'.'*((N-4)//2-i)}{'('*i}....{')'*i}{'.'*((N-4)//2-i)}")
         for i in steps]
print()

# If matplotlib is available, we can plot the above computed values.
def myplot():
    import matplotlib
    import matplotlib.pylab as plt
    
    plt.figure(1)
    plt.plot(steps,path1,'go', label='Outside bases first')
    plt.plot(steps,path1,'g-', label='_nolabel_')
    plt.plot(steps,path2,'rx', label='Inside bases first')
    plt.plot(steps,path2,'r-', label='_nolabel_')
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
