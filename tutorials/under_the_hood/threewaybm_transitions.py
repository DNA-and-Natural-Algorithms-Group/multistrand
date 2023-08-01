# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
This example is similar to hairpin_transistion.py, except that multistranded
complexes are handles. Mainly, we think about what the exact and loose
macrostate differences are -- what do they tell us?

Try it like this, e.g.:
  python -i threewaybm_transitions.py
"""

import numpy as np

from multistrand.objects import *
from multistrand.options import Options
from multistrand.system import SimSystem


# for StopCondition and Macrostate definitions:
Exact_Macrostate = 0   # match a secondary structure exactly (i.e. any system state that has a complex with this exact structure)
Bound_Macrostate = 1   # match any system state in which the given strand is bound to another strand
Dissoc_Macrostate = 2  # match any system state in which there exists a complex with exactly the given strands, in that order
Loose_Macrostate = 3   # match a secondary structure with "don't care"s, allowing a certain number of disagreements
Count_Macrostate = 4   # match a secondary structure, allowing a certain number of disagreements
# see Schaeffer's PhD thesis, chapter 7.2, for more information


def setup_options_threewaybm(toehold_seq = "GTGGGT", bm_design = "ACCGCACGTCACTCACCTCG"):

    # the structures below are hard-coded for these lengths
    assert len(toehold_seq)==6
    assert len(bm_design)==20

    toehold = Domain(name="toehold",sequence=toehold_seq,length=6)
    branch_migration = Domain(name="bm", sequence=bm_design, seq_length=20)
    
    substrate = toehold + branch_migration
    incumbent = Strand(name="incumbent",domains=[branch_migration.C])

    incoming = substrate.C

    # start with 6-base toehold fully bound
    start_complex = Complex(strands=[incoming, substrate, incumbent], structure=".(+)(+)")

    initial_structure_dp   = "....................((((((+))))))((((((((((((((((((((+))))))))))))))))))))"
    six_bases_structure_dp = "..............((((((((((((+))))))))))))((((((((((((((+))))))))))))))......"
    six_bases_loose_dp     = "**************((**********+**********))((************+************))******"
    twelve_bases_struc_dp  = "........((((((((((((((((((+))))))))))))))))))((((((((+))))))))............"
    twelve_bases_loose_dp  = "********((*****************+***************))((******+******))************"
    eighteen_structure_dp  = "..((((((((((((((((((((((((+))))))))))))))))))))))))((+)).................."
    eighteen_loose_dp      = "**((**********************+**********************))((+))******************"

    six_bases_complex           = Complex(strands=[incoming,substrate,incumbent], structure=six_bases_structure_dp)
    twelve_bases_complex        = Complex(strands=[incoming,substrate,incumbent], structure=twelve_bases_struc_dp)
    eighteen_bases_complex      = Complex(strands=[incoming,substrate,incumbent], structure=eighteen_structure_dp)
    six_basesloose_complex      = Complex(strands=[incoming,substrate,incumbent], structure=six_bases_loose_dp)
    twelve_basesloose_complex   = Complex(strands=[incoming,substrate,incumbent], structure=twelve_bases_loose_dp)
    eighteen_basesloose_complex = Complex(strands=[incoming,substrate,incumbent], structure=eighteen_loose_dp)

    disassoc_complex            = Complex(strands=[incumbent], structure=".")   # succesful strand displacement
    failed_complex              = Complex(strands=[incoming], structure="..")   # failed strand displacement attempt

    start_sc          = Macrostate("INITIAL", [(start_complex,Count_Macrostate,2)])                 # Within distance 2 of the start_complex state.
    six_sc_exact      = Macrostate("SIX_EXACT", [(six_bases_complex,Exact_Macrostate,0)])           # the third parameter is ignored; not needed for exact macrostates
    six_sc_loose      = Macrostate("SIX_LOOSE", [(six_basesloose_complex,Loose_Macrostate,2)])      # 8 base pairs defined; must have at least 6 to match.
    twelve_sc_exact   = Macrostate("TWELVE_EXACT", [(twelve_bases_complex,Exact_Macrostate,0)])
    twelve_sc_loose   = Macrostate("TWELVE_LOOSE", [(twelve_basesloose_complex,Loose_Macrostate,2)])
    eighteen_sc_exact = Macrostate("EIGHTEEN_EXACT", [(eighteen_bases_complex,Exact_Macrostate,0)])
    eighteen_sc_loose = Macrostate("EIGHTEEN_LOOSE", [(eighteen_basesloose_complex,Loose_Macrostate,2)])

    # why bother giving a list of just one macrostate-def tuple?  A Macrostate with a list of multiple tuples give the AND (i.e. intersection) of microstates.

    completed_sc      = StopCondition("stop:COMPLETE", [(disassoc_complex,Dissoc_Macrostate,0)])  # incumbent strand fell off  
    rejected_sc       = StopCondition("stop:REJECTED", [(failed_complex,Dissoc_Macrostate,0)])    # incoming strand fell off

    # join_concentration is not defined, because in this simulation we stop before there's any chance for association steps
    o_exact = Options(simulation_mode="Transition",parameter_type="Nupack", dangles="Some",
                substrate_type="DNA", num_simulations = 5, simulation_time=.01, temperature=310.15,
                start_state=[start_complex], rate_scaling='Calibrated', verbosity=0)
    o_exact.stop_conditions = [start_sc, six_sc_exact, twelve_sc_exact, eighteen_sc_exact, completed_sc, rejected_sc]
    o_loose = Options(simulation_mode="Transition",parameter_type="Nupack", dangles="Some",
                substrate_type="DNA", num_simulations = 5, simulation_time=.01, temperature=310.15,
                start_state=[start_complex], rate_scaling='Calibrated', verbosity=0)

    o_loose.stop_conditions = [start_sc, six_sc_loose, twelve_sc_loose, eighteen_sc_loose, completed_sc, rejected_sc]
    return o_exact, o_loose


# mol will be a list of True/False for which transition macrostates the system has entered
# so in_state(mol) returns True if the system is in at least one of the listed macrostates.
def in_state( mol ): return sum(mol) > 0


# mol is a Boolean descriptor of macrostate occupancy, like mol above.
# a short-hand name for this macrostate (based on the order given in stop_conditions) is provided.
def mol_name(mol):
    charindex = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0'
    names = [charindex[j] for i,j in zip(mol,range(len(mol))) if i]
    if names == []:
        names = charindex[26]
    else:
        names = ",".join(names)
    return names


# t0 and t1 are Boolean descriptors of macrostate occupancy, like mol above.
# here, we provide a printable name for the transition between two macrostate occupancy lists.
def trans_name(t0,t1):
    return mol_name(t0) + ' -> ' + mol_name(t1)


def print_transitions( transition_traj ):
    for t in transition_traj:
        print("%12g : %s" % ( t[0], mol_name(t[1]) ))


# for each simulation, the transition trajectory reports the tuple (time_entered, which_macrostates_the_system_is_now_in)
def parse_transition_lists( transition_traj_list ):
    transition_dict = {}

    # the mol1 --> mol2 transition times represent (time of entry into mol1) to (time of entry into mol2)
    for transition_traj in transition_traj_list:
        truncated = [i for i in transition_traj if in_state(i[1])]
        tt = truncated # only keep the entry times and mol states for non-trivial mols

        for i in range(len(tt)-1):
            nm = trans_name(tt[i][1],tt[i+1][1])
            if nm in transition_dict:
                transition_dict[nm].append( tt[i+1][0] - tt[i][0] )
            else:
                transition_dict[nm] = [tt[i+1][0] - tt[i][0]]

    return transition_dict


def parse_transition_list( transition_traj_list ):
    return parse_transition_lists( [transition_traj_list] )

    
def print_transition_dict( transition_dict, options = None ):
    k = list(transition_dict.keys())
    k.sort() 

    for i in k:
        transition_times = np.array( transition_dict[i] )
        print(("{0}: {2:.2e} ({1})".format(i,len(transition_dict[i]),np.mean(transition_times))))
    
    # also print the true names of the macrostates, if an Options object is provided
    charindex = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0'
    if options:
        for i,idx in zip(options.stop_conditions,range(len(options.stop_conditions))):
            print(("{0}: {1}".format( i.tag, charindex[idx])))


#### Stuff to try... automatically or by hand, line-by-line 

# takes a little time to run, but not too long.
print("--- Running Simulations ---")
#o1,o2 = setup_options_threewaybm()
o1,o2 = setup_options_threewaybm(bm_design="ACCGCACGTCCACGGTGTCG")  # try this too.  toehold dissociates?... much slower!
s=SimSystem(o1)
s.start()
s=SimSystem(o2)
s.start()
# same energy model in both simulations, so no need to update the energy model in between.

print()
print("--- Analysis of simulations with exact transition states ---")
# print "  Coarse-grained trajectory of simulation #1:"
# print_transitions(o1.interface.transition_lists[0])
print("  Transitions from simulation #1:")
parsedlist = parse_transition_list(o1.interface.transition_lists[0])
print_transition_dict(parsedlist)
print("  Transitions averaged over all %d simulations:" % o1.num_simulations)
parsedlist = parse_transition_lists(o1.interface.transition_lists)
print_transition_dict(parsedlist,o1) # adds names for macrostates

print()
print("--- Analysis of simulations with loose transition states ---")
# print "  Coarse-grained trajectory of simulation #1:"
# print_transitions(o2.interface.transition_lists[0])
print("  Transitions from simulation #1:")
parsedlist = parse_transition_list(o2.interface.transition_lists[0])
print_transition_dict(parsedlist)
print("  Transitions averaged over all %d simulations:" % o2.num_simulations)
parsedlist = parse_transition_lists(o2.interface.transition_lists)
print_transition_dict(parsedlist,o2)  # adds names for macrostates

# The lesson here is similar to hairpin_transitions.py:  exact macrostates don't track system behavior well, 
# because too much of the time the system is not in any macrostate -- thus it's easy to avoid the "checkpoints"
# and go straight from start to finish.
