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
from multistrand.options import Options, EnergyType
from multistrand.system import SimSystem, calculate_energy


def create_setup():
    # build complexes with domain-level information
    toehold_seq = "GTGGGT"
    bm_design_A = "ACCGCACGTCCACGGTGTCG"
    bm_design_B = "ACCGCACGTCACTCACCTCG"

    toehold = Domain(name="toehold", sequence=toehold_seq)
    branch_migration_A = Domain(name="bm_A", sequence=bm_design_A)
    branch_migration_B = Domain(name="bm_B", sequence=bm_design_B)

    substrate_A = toehold + branch_migration_A
    substrate_B = toehold + branch_migration_B
    incumbent_A = Strand(name="incumbent", domains=[branch_migration_A.C])
    incumbent_B = Strand(name="incumbent", domains=[branch_migration_B.C])

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
    o1.simulation_mode = "Trajectory"
    o1.num_simulations = 1
    o1.simulation_time = 5e-5
    o1.temperature = 37.0
    o1.dangles = 1
    # record every 1000 steps (so we'll get around 100 record entries)
    o1.output_interval = 1000
    o1.start_state = [start_complex_A]
    o1.join_concentration = 1e-6  # 1 uM
    o1.JSMetropolis37()
    o1.verbosity = 0

    o2 = Options()
    o2.simulation_mode = "Trajectory"
    o2.num_simulations = 1
    o2.simulation_time = 5e-5
    o2.temperature = 37.0
    o2.dangles = 1
    o2.start_state = [start_complex_B]
    # could do o2.output_time to get trajectory items evenly spaced in time
    # instead of by number of steps
    o2.output_interval = 1000
    o2.join_concentration = 1e-6  # 1 uM
    o2.JSMetropolis37()
    o2.verbosity = 0

    return o1, o2


# generalized from "hairpin trajectories tutorial" version to allow multistrand complexes and multiple complexes in a tube
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
            dGv += calculate_energy( [Complex( strands=[Strand(sequence=s) for s in cs], structure=st)], o, EnergyType.tube)[0]
        if not dGv == dG: print("Energy Mismatch")


# Perform the simulations
o1,o2 = create_setup()
s1 = SimSystem(o1)
s1.start()
s2 = SimSystem(o2)
s2.start()

# Show what happened
print()
print("Sequence Design 1 (shown every 1000 steps):")
print_trajectory(o1)
print()
print("Sequence Design 2 (shown every 1000 steps):")
print_trajectory(o2)


if __name__ == '__main__':
    print()
    print("Can you tell the difference between design 1 and design 2?")
    print("On most runs, Design 2 finishes within the allotted time.")
    print("That is, the incumbent strand gets displaced.")
    print("You can see that a \"+\" turns into a \" \", indicating that separate complexes are present.")
    print("Or did the incoming strand get rejected, the toehold spontaneously dissociating?")
    print("Can you tell?")
    print("And can you discern why Design 1 almost never completes displacement in the given time?")

### Some notes for the intrepid explorer:

# Once `s.start()` has run, you cannot run `start()` for this `SimSystem` again.
# You must make a new `Options` object, or otherwise retrieve the simulation
# results before you reuse the old `Options` object, and then create a new
# `SimSystem` object.

# Note that if `num_simulations > 1`, then all the recorded states get collected
# together into the `full_trajectory`. To tell where one trajectory stops and
# the next simulation starts, you can look at the time stamps.

# A good exercise to the reader, as a test of understanding: You should be able
# to extract sequences & structures from a trajectory, make a `Complex` of them,
# evaluate the energy, start new simulations, etc.

# Trajectory like this were used in Schaeffer's MS thesis. Note that in the
# thesis, we displayed them in 5'->3' counterclockwise. Oops.
