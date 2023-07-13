# hairpin_trajectories.py
#
# This is a follow-up on hairpin_energies.py.
# Here, rather than manually examining secondary structures of interest, we let a simulation run explore the full secondary structure 
# energy landscape according to a Metropolis-biased random walk.
# 
#
# Try it like this:
#   python -i hairpin_trajectories.py
#   ipython --pylab -i hairpin_trajectories.py

if False:  # only needed if you're having trouble with your Multistrand installation
    import multistrand_setup

try:
    from multistrand.objects import *
    from multistrand.options import Options, Literals
    from multistrand.system import SimSystem, energy

except ImportError:
    print("Could not import Multistrand.")
    raise

#############

# More meaningful names for argument values to the energy() function call, below.
Loop_Energy = 0    # requesting no dG_assoc or dG_volume terms to be added.  So only loop energies remain.
Volume_Energy = 1  # requesting dG_volume but not dG_assoc terms to be added.  No clear interpretation for this.
Complex_Energy = 2 # requesting dG_assoc but not dG_volume terms to be added.  This is the NUPACK complex microstate energy, sans symmetry terms.
Tube_Energy = 3    # requesting both dG_assoc and dG_volume terms to be added.  Summed over complexes, this is the system state energy.


# Results of the simulation are stored in the Options object 'o' that was used to set up the simulation.
# Since this is a "Trajectory Mode" simulation, we will extract the sequence of conformations visited, and print them.
# Assumes a single strand is being simulated.
def print_trajectory(o):
    print(o.full_trajectory[0][0][3])   # the strand sequence
    print(o.start_state[0].structure)   # the starting structure
    for i in range(len(o.full_trajectory)):
        time = o.full_trajectory_times[i]
        state = o.full_trajectory[i][0]
        struct = state[4]
        dG = state[5]
        print(f'{struct} t={time:11.9f} seconds, dG={dG:6.2f} kcal/mol')

# Sequence is from Schaeffer's PhD thesis, chapter 7, figure 7.1 -- start with no base pairs formed.
c = Complex( strands=[Strand(name="hairpin", sequence="GTTCGGGCAAAAGCCCGAAC")], structure= 20*'.' )
# WARNING!  Unfortunately, Multistrand currently does not test to check that the requested structure is valid.
# Providing an invalid structure is likely to lead to a core dump or segmentation fault.

o = Options(temperature=25.0, 
            dangles='Some', 
            start_state = [c], 
            simulation_time = 0.0000001,  # 0.1 microseconds
            num_simulations = 1,  # don't play it again, Sam
            output_interval = 1,  # record every single step
            rate_method = Literals.metropolis, # the default is 'Kawasaki' (numerically, these are 1 and 2 respectively)
            rate_scaling = 'Calibrated', # this is the same as 'Default'.  'Unitary' gives values 1.0 to both.  
            simulation_mode = 'Trajectory')  # numerically 128.  See interface/_options/constants.py for more info about all this.

print(f"k_uni = {o.unimolecular_scaling:g} /s, k_bi = {o.bimolecular_scaling:g} /M/s")  # you can also set them to other values if you want

# This actually runs the simulation.  
s = SimSystem(o)
s.start()
# Important caveat:  SimSystem will initialize the energy model according to information in Options 'o' 
# if the energy model has not yet been initialized.
# But if prior calls have already initialized the energy model -- even if it's at another temperature or join_concentration -- then
# it will not be automatically re-initialized.  You would have to do this manually.

print_trajectory(o)        
# Note that the simulation proceeds until the time limit has been EXCEEDED.
# That means that, at the exact time specified, the system is in the PENULTIMATE state. 
# Just FYI -- but this is important if you are sampling to get an "equilibrium" or time-dependent sample.
# If you were to take the last state, and if that state is very short-lived, then you would over-sample it.

def myplot():
    import numpy as np
    import matplotlib
    import matplotlib.pylab as plt

    times = o.full_trajectory_times
    states = o.full_trajectory
    energies = [s[0][5] for s in states]  # you can examine 'states' to see what other useful information is there

    plt.figure(1)
    plt.plot(times,energies,'go', times,energies,'g-')
    plt.title("Energy landscape for simulated hairpin folding trajectory")
    plt.xlabel("Time (seconds)",fontsize='larger')
    plt.ylabel("Microstate Energy (kcal/mol)",fontsize='larger')
    plt.yticks(fontsize='larger',va='bottom')
    plt.xticks(fontsize='larger')
    plt.show()


# execute the following only if run from the command line, e.g. as "python -i hairpin_trajectories.py" or "ipython --pylab -i hairpin_trajectories.py"
if __name__ == '__main__':
    myplot()
else:
    print("Try:\nhairpin_trajectories.myplot()\n")
