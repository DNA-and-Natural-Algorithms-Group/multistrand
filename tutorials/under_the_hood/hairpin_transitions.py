# hairpin_transitions.py
# 
# Having seen how the simulation generates trajectories in hairpin_trajectories.py, here we extract some information from the trajectories.
#
# Using the example of Chapter 7.3 of Schaeffer's PhD thesis, this script shows how to use macrostate definitions and Transition Mode 
# to extract a coarse-grained description of the system's motion through state space.
#
# This mode requires a little more post-processing of the Multistrand simulations... but don't be put off.  It's easy.
#
# Note that the transition times listed here are briefer than those that appear in the tables of Chapter 7.3 of Schaeffer's thesis.
# That is because we are using different unimolecular and bimolecular rate constant values.
#
#
# Try it like this, e.g.:
#   python -i hairpin_transitions.py

if False:  # only needed if you're having trouble with your Multistrand installation
    import multistrand_setup

try:
    from multistrand.objects import *
    from multistrand.options import Options
    from multistrand.system import SimSystem

except ImportError:
    print("Could not import Multistrand.")
    raise

import numpy as np

#############

# for StopCondition and Macrostate definitions:
Exact_Macrostate = 0   # match a secondary structure exactly (i.e. any system state that has a complex with this exact structure)
Bound_Macrostate = 1   # match any system state in which the given strand is bound to another strand
Dissoc_Macrostate = 2  # match any system state in which there exists a complex with exactly the given strands, in that order
Loose_Macrostate = 3   # match a secondary structure with "don't care"s, allowing a certain number of disagreements
Count_Macrostate = 4   # match a secondary structure, allowing a certain number of disagreements
# see Schaeffer's PhD thesis, chapter 7.2, for more information

def setup_options_hairpin():
    """
    setup_options_hairpin( )

    Returns the options object for simple hairpin example of
    transition mode in Chapter 7.3 of Schaeffer's PhD thesis.
    """

    # Once domains are defined, strands can be built from them using "+".
    stem = Domain(name="stem",sequence="GCATGC",length=6)
    hp = Domain(name="hairpin", sequence="AAAA",length=4)
    s = stem + hp + stem.C
    
    # Note that because we defined domains, we can either represent secondary structures at the domain level or at the sequence level.
    start_complex = Complex(strands=[s], structure="...")
    pathway_endside_exact_complex = Complex(strands=[s], structure="(((..........)))")
    pathway_endside_loose_complex = Complex(strands=[s], structure="(((**********)))")
    pathway_hpside_exact_complex = Complex(strands=[s],  structure="...(((....)))...")
    pathway_hpside_loose_complex = Complex(strands=[s],  structure="***(((****)))***")
    full_complex = Complex( strands=[s], structure="(.)")

    # Define macrostates to be tracked, i.e. we look for transitions between them.
    initial_sc           = Macrostate( "INITIAL", [(start_complex,Exact_Macrostate,0)])
    pathway_hp_exact_sc  = Macrostate( "HPSIDE_EXACT", [(pathway_hpside_exact_complex,Exact_Macrostate,0)])
    pathway_hp_loose_sc  = Macrostate( "HPSIDE_LOOSE", [(pathway_hpside_loose_complex,Loose_Macrostate,2)])   # within distance 2 
    pathway_end_exact_sc = Macrostate( "ENDSIDE_EXACT", [(pathway_endside_exact_complex,Exact_Macrostate,0)])
    pathway_end_loose_sc = Macrostate( "ENDSIDE_LOOSE", [(pathway_endside_loose_complex,Loose_Macrostate,2)]) # within distance 2 
    full_sc              = StopCondition( "stop:FULL", [(full_complex,Exact_Macrostate,0)])
    # Multistrand treats Macrostates and StopConditions interchangeably; here we choose the name as a mnemonic for how it will be used.
    # The simulation will stop the first time that 'full_sc' is reached, so we call it a StopCondition.
    # The simulation keeps track of where it is, but keeps going, when it visits the others -- so they are called Macrostates.
    # But the simulation would proceed identically if randomly called some Macrostates and others StopConditions.  
    # What determines Multistrand's behavior is that the name of 'full_sc' begins with "stop:" -- this is what causes the simulation to stop when 'full_sc' is reached.

    # We will set up two Transition Mode simulations, one looking at transitions between exact states...
    o_exact = Options(simulation_mode="Transition",parameter_type="Nupack",substrate_type="DNA", temperature=310.15,
                num_simulations = 1000, simulation_time=.01, start_state=[start_complex], rate_scaling='Calibrated', verbosity=0)
    o_exact.stop_conditions = [initial_sc, pathway_end_exact_sc, pathway_hp_exact_sc,full_sc]

    # ... and one looking at transitions between loosely defined macrostates
    o_loose = Options(simulation_mode="Transition",parameter_type="Nupack",substrate_type="DNA", temperature=310.15,
                num_simulations = 1000, simulation_time=.01, start_state=[start_complex], rate_scaling='Calibrated', verbosity=0)
    o_loose.stop_conditions = [initial_sc, pathway_end_loose_sc, pathway_hp_loose_sc,full_sc]

    # change verbosity to 1, above, and you'll see a print-out during the simulation runs, every time a stop state is reached.

    return o_exact,o_loose

#### Some helper code follows.

# The following naming routines are motivated by the observation that when "loose" macrostates are used,
# it is possible that a given system state may belong to more than one macrostate.
# Thus, transitions are not between macrostates per se, but are between macrostate membership vectors.
# I.e. if the macrostates are A, B, C, and D, then a single simulation step might take the system state
# from being in both A and C, to being in B, C, and D.  And perhaps thence to being in no macrostates.

# mol will be a list of True/False for which transition macrostates the system has entered,
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
        print "%12g : %s" % ( t[0], mol_name(t[1]) )
                  
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
    k = transition_dict.keys()
    k.sort() 

    for i in k:
        transition_times = np.array( transition_dict[i] )
        print("{0}: {2:.2e} ({1})".format(i,len(transition_dict[i]),np.mean(transition_times)))
    
    # also print the true names of the macrostates, if an Options object is provided
    charindex = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0'
    if options:
        for i,idx in zip(options.stop_conditions,range(len(options.stop_conditions))):
            print("{0}: {1}".format( i.tag, charindex[idx]))


#### Back to stuff to try... automatically or by hand, line-by-line 

print "--- Running Simulations ---"
o1,o2 = setup_options_hairpin()
s=SimSystem(o1)
s.start()
s=SimSystem(o2)  
s.start()
# the second simulation is at the same join_concentration and temperature as the previous simulation,
# so there's no need to re-initialize the energy model.

# Try looking at o1.interface.transition_lists[0] to see each transition from the first simulation.

print
print "--- Analysis of simulations with exact transition states ---"
print "  Coarse-grained trajectory of simulation #1:"
print_transitions(o1.interface.transition_lists[0])
print "  Transitions from simulation #1:"
parsedlist = parse_transition_list(o1.interface.transition_lists[0])
print_transition_dict(parsedlist)
print "  Transitions averaged over all %d simulations:" % o1.num_simulations
parsedlist = parse_transition_lists(o1.interface.transition_lists)
print_transition_dict(parsedlist,o1) # adds names for macrostates

print
print "--- Analysis of simulations with loose transition states ---"
print "  Coarse-grained trajectory of simulation #1:"
print_transitions(o2.interface.transition_lists[0])
print "  Transitions from simulation #1:"
parsedlist = parse_transition_list(o2.interface.transition_lists[0])
print_transition_dict(parsedlist)
print "  Transitions averaged over all %d simulations:" % o2.num_simulations
parsedlist = parse_transition_lists(o2.interface.transition_lists)
print_transition_dict(parsedlist,o2)  # adds names for macrostates

# Ponder why, for exact macrostates, the A -> D transition can occur and is not even rare, 
# while with loose macrostates, A transitions exclusively either to B or to C.
