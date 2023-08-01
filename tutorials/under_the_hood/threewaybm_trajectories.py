# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
This example illustrates how to examine trajectories in multistranded simulations,
where association and dissociation can change the number and type of complexes.
Pay attention to how, when that happens, the strand ordering can change,
and thus dot-paren structures must be displayed differently.

Try it like this, e.g.:
  python -i hairpin_trajectories.py
"""

from multistrand.objects import *
from multistrand.options import Options
from multistrand.system import SimSystem, energy


# More meaningful names for argument values to the energy() function call, below.
Loop_Energy = 0    # requesting no dG_assoc or dG_volume terms to be added.  So only loop energies remain.
Volume_Energy = 1  # requesting dG_volume but not dG_assoc terms to be added.  No clear interpretation for this.
Complex_Energy = 2 # requesting dG_assoc but not dG_volume terms to be added.  This is the NUPACK complex microstate energy, sans symmetry terms.
Tube_Energy = 3    # requesting both dG_assoc and dG_volume terms to be added.  Summed over complexes, this is the system state energy.


def create_setup():

    # build complexes with domain-level information
    toehold_seq = "GTGGGT"
    bm_design_A = "ACCGCACGTCCACGGTGTCG"
    bm_design_B = "ACCGCACGTCACTCACCTCG"

    toehold = Domain(name="toehold",sequence=toehold_seq,length=6)
    branch_migration_A = Domain(name="bm_A", sequence=bm_design_A, seq_length=20)
    branch_migration_B = Domain(name="bm_B", sequence=bm_design_B, seq_length=20)
    
    substrate_A = toehold + branch_migration_A
    substrate_B = toehold + branch_migration_B
    incumbent_A = Strand(name="incumbent",domains=[branch_migration_A.C])
    incumbent_B = Strand(name="incumbent",domains=[branch_migration_B.C])

    incoming_A = substrate_A.C
    incoming_B = substrate_B.C

    # Note that "+" is used to indicate strand breaks.  
    # So the initial structures represent the incoming strand bound by its toehold,
    # and we'll see that either it completes strand displacement, or it dissociates.
    start_complex_A = Complex(strands=[incoming_A, substrate_A, incumbent_A],
                              structure=".(+)(+)")
    start_complex_B = Complex(strands=[incoming_B, substrate_B, incumbent_B],
                              structure=".(+)(+)")
    
    o1 = Options()
    o1.simulation_mode = 0x0080 # trajectory mode
    o1.num_simulations = 1
    o1.simulation_time = 0.00002 # 200 microseconds, about 250 steps
    o1.temperature = 37.0
    o1.dangles = 1
    o1.output_interval = 100   # record every 100 steps (so we'll get around 100 record entries)
    o1.start_state = [start_complex_A]
    o1.rate_scaling="Calibrated"
    o1.join_concentration=1e-6  # 1 uM 
    o1.verbosity=0  # doesn't turn off output during simulation -- but it should.  please wait for multistrand 3.0.
                    # the alternative is to increase the output_interval so something ridiculously high, but this also eliminates the trajectory record

    o2 = Options()
    o2.simulation_mode = 0x0080  # trajectory mode
    o2.num_simulations = 1
    o2.simulation_time = 0.00002 # 200 us, about 250 steps
    o2.temperature = 37.0
    o2.dangles = 1
    o2.start_state = [start_complex_B]
    o2.output_interval = 100   # could do o2.output_time to get trajectory items evenly spaced in time instead of by number of steps
    o2.rate_scaling="Calibrated"
    o1.join_concentration=1e-6  # 1 uM 
    
    return o1,o2


# generalized from "hairpin_trajectories.py" version to allow multistrand complexes and multiple complexes in a tube
def print_trajectory(o):
    seqstring=''
    for i in range(len(o.full_trajectory)): # go through each output microstate of the trajectory
        time = o.full_trajectory_times[i]   # time at which this microstate is entered
        states = o.full_trajectory[i]       # this is a list of the complexes present in this tube microstate
        newseqs = []
        for state in states: newseqs += [ state[3] ]   # extract the strand sequences in each complex (joined by "+" for multistranded complexes)
        newseqstring = ' '.join(newseqs)    # make a space-separated string of complexes, to represent the whole tube system sequence
        if not newseqstring == seqstring : 
            print(newseqstring)
            seqstring=newseqstring          # because strand order can change upon association of dissociation, print it when it changes
        structs = []
        for state in states: structs += [ state[4] ]   # similarly extract the secondary structures for each complex
        tubestruct = ' '.join(structs)      # give the dot-paren secondary structure for the whole test tube
        dG=0
        for state in states: dG += state[5]
        print('%s t=%11.9f seconds, dG=%6.2f kcal/mol' % (tubestruct,time, dG))

        # Needlessly verify that the reported trajectory energies are the Tube_Energy values
        dGv=0
        for state in states:
            cs=state[3].split('+')
            st=state[4]
            dGv += energy( [Complex( strands=[Strand(sequence=s) for s in cs], structure=st)], o, Tube_Energy)[0]  
        if not dGv == dG: print("Energy Mismatch")


# Perform the simulations
o1,o2 = create_setup()
s1 = SimSystem(o1)
s1.start()
s2 = SimSystem(o2)
s2.start()
# see below about the energy model

# Show what happened
print()
print("Sequence Design 1 (shown every 100 steps):")
print_trajectory(o1)
print()
print("Sequence Design 2 (shown every 100 steps):")
print_trajectory(o2)
        

if __name__ == '__main__':
    print("Can you tell the difference between design 1 and design 2?")
    print("On most runs, Design 2 finishes within the allotted time.")
    print("That is, the incumbent strand gets displaced.")
    print("You can see that a \"+\" turns into a \" \", indicating that separate complexes are present.")
    print("Or did the incoming strand get rejected, the toehold spontaneously dissociating?")
    print("Can you tell?")
    print("And can you discern why Design 1 almost never completes displacement in the given time?")

# Some notes for the intrepid explorer: 

# The Options object and the SimSystem object are intimately tied, and after s.start() you can't just change energy model parameters or rate parameters.
# If you do, you need to tell Multistrand to update the energy model call, initialize_energy_model(o), and then make the SimSystem and start it.
# In the example, if these two simulations had used a different join_concentration, a different temperature, a different material (RNA vs DNA), 
# a different dangles option, a different parameter set (NUPACK vs Vienna), a different rate method (Kawasaki or Metropolis), or a different set of rate 
# scaling paramters (e.g. Calibrated, Unitary, etc), then we would need to use initialize_energy_model() in between.
# This is because all simulations share the same energy & kinetics model (which we just call the "energy model" for short)
# and it it not automatically re-adjusted after it has been created.  Note that "same energy model" means "the same sequence and secondary structure will get
# the same numerical value for the energy, and the same move will get the same rate" -- so, even though the underlying model may be the same, a change in, 
# e.g. concentration will imply a "different energy model" as we are using the term here.
#
# system.SimSystem(), system.calculate_energy(), system.calculate_rate() will all call initialize_energy_model() if none has been created yet,
# or else use the existing one even if the passed Options object specified different temperature or other parameters.  
# The onus is on the user to make sure this is all correct, if your code ever changes energy / kinetics parameters.  
# Basically, the only case where you *don't* need to update the energy model is if you change simulation mode, simulation time, start states, stop conditions,
# or number of simulations.
#
# This is a common source of mistakes, so why not always do it automatically?  The reason is so that Multistrand can efficiently exploit multicore processors
# and run with multiple threads, as shown in the threewaybm_first_step_mode.py example.  All threads share the energy model, and therefore we must not
# change the energy model when a new simulation starts but existing simulations are already running.  

# Also, once s.start() has run, you cannot run ('start') for this SimSystem again.  You must make a new Options object and a new SimSystem object.
# But if, in doing so, you changed energy/kinetics model parameters, you must also initialize_energy_model(o), as we've said.

# Note that if num_simulations > 1, then all the recorded states get collected together into the full_trajectory.  
# To tell where one trajectory stops and the next simulation starts, you can look at the time stamps.

# A good exercise to the reader, as a test of understanding: 
# you should be able to extract sequences & structures from trajectory, make a Complex of them, evaluate energy, start new simulation, etc.

# Trajectory like this were used in Schaeffer's MS thesis.
# Note that in the thesis, we displayed them in 5'->3' counterclockwise. Oops.
