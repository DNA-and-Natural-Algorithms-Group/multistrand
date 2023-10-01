# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
This is a follow-up on hairpin_energies.py. Here, rather than manually examining
secondary structures of interest, we let a simulation run explore the full
secondary structure energy landscape according to a Metropolis-biased random
walk.

Try it like this:
  python -i hairpin_trajectories.py
  ipython --pylab -i hairpin_trajectories.py
"""

from multistrand.objects import *
from multistrand.options import Options
from multistrand.system import SimSystem


# Results of the simulation are stored in the Options object 'o' that was used
# to set up the simulation. Since this is a "Trajectory Mode" simulation, we
# will extract the sequence of conformations visited, and print them. Assumes a
# single strand is being simulated.
def print_trajectory(o):
    print(o.full_trajectory[0][0][3])   # the strand sequence
    print(o.start_state[0].structure)   # the starting structure
    for i in range(len(o.full_trajectory)):
        time = o.full_trajectory_times[i]
        state = o.full_trajectory[i][0]
        struct = state[4]
        dG = state[5]
        print(f'{struct} t={time:11.9f} seconds, dG={dG:6.2f} kcal/mol')

# Sequence is from Schaeffer's PhD thesis, chapter 7, figure 7.1 -- start with
# no base pairs formed.
c = Complex(strands=[Strand(name="hairpin", sequence="GTTCGGGCAAAAGCCCGAAC")],
            structure= 20*'.')
# WARNING! Unfortunately, Multistrand currently does not test to check that the
# requested structure is valid. Providing an invalid structure is likely to lead
# to a core dump or segmentation fault.

o = Options(temperature=25.0, 
            dangles='Some', 
            start_state=[c],
            simulation_time=5e-6,  # 5 microseconds
            num_simulations=1,  # don't play it again, Sam
            output_interval=1,  # record every single step
            simulation_mode='Trajectory')
o.DNA23Metropolis()
print(f"k_uni = {o.unimolecular_scaling:g} /s, k_bi = {o.bimolecular_scaling:g} /M/s")
# you can also set them to other values if you want

# This actually runs the simulation.
s = SimSystem(o)
s.start()

print_trajectory(o)        
# Note that the simulation proceeds until the time limit has been EXCEEDED.
# That means that, at the exact time specified, the system is in the PENULTIMATE
# state. Just FYI -- but this is important if you are sampling to get an
# "equilibrium" or time-dependent sample. If you were to take the last state,
# and if that state is very short-lived, then you would over-sample it.

def myplot():
    import matplotlib
    import matplotlib.pylab as plt

    times = o.full_trajectory_times
    states = o.full_trajectory
    energies = [s[0][5] for s in states]  # you can examine 'states' to see what other useful information is there

    plt.figure(1)
    plt.plot(times, energies,'go', times,energies,'g-')
    plt.title("Energy landscape for simulated hairpin folding trajectory")
    plt.xlabel("Time (seconds)", fontsize='larger')
    plt.ylabel("Microstate Energy (kcal/mol)", fontsize='larger')
    plt.yticks(fontsize='larger', va='bottom')
    plt.xticks(fontsize='larger')
    plt.show()


# execute the following only if run from the command line, e.g. as
# "python -i hairpin_trajectories.py" or
# "ipython --pylab -i hairpin_trajectories.py"
if __name__ == '__main__':
    myplot()
else:
    print("Try:\nhairpin_trajectories.myplot()\n")
