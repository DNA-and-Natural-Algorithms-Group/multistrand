####################################################################
#                                                                  #
#  Copyright (c) 2010-2015 California Institute of Technology.     #
#  Distributed under the MIT License.                              #
#  (See accompanying file LICENSE or copy at                       #
#  http://opensource.org/licenses/MIT)                             #
#                                                                  #
####################################################################
#                                                                  #
#   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)         #
#                                                                  #
####################################################################

import sys,os,os.path
import cPickle
import numpy as np

multihome = None
if 'MULTISTRANDHOME' in os.environ:
    if not os.path.isfile( os.path.join( os.environ['MULTISTRANDHOME'], 'setup.py') ):
        warnings.warn( ImportWarning("Could not find the file 'setup.py' in your MULTISTRANDHOME [{0}]; this environment variable is possibly out of date or not referring to the new Mulistrand distribution."))
        multihome=None
    else:
        if os.environ['MULTISTRANDHOME'] not in sys.path:
            multihome= os.environ['MULTISTRANDHOME']     

if multihome != None:
    sys.path.append(multihome)

try:

    from multistrand.objects import *
    from multistrand.options import Options
    from multistrand.system import SimSystem

except ImportError:
    print("Could not import Multistrand.")
    raise

def setup_options_hp_test():
    """
    setup_options_hp_test( )

    returns the options object for our simple hairpin test for
    transition mode.
    """

    stem = Domain(name="stem",sequence="AGGACT",length=6)
    hp = Domain(name="hairpin", sequence="AAAA",length=4)
    
    s = stem + hp + stem.C
    
    start_complex = Complex(strands=[s], structure="...")

    pathway_complex = Complex(strands=[s], structure="...(((....)))...")
    pathway_loose_complex = Complex(strands=[s], structure="***(((....)))***")
    
    full_complex = Complex( strands=[s], structure="(.)")
    full_loose_complex = Complex( strands=[s], structure="((((((....))))))")

    initial_sc = StopCondition( "INITIAL", [(start_complex,0,0)])
    pathway_sc = StopCondition( "PATHWAY", [(pathway_complex,0,0)])
    pathway_loose_sc = StopCondition( "PATHWAYLOOSE", [(pathway_loose_complex,3,0)])
    full_sc    = StopCondition( "stop:FULL", [(full_complex,0,0)])
    full_loose_sc = StopCondition( "stop:FULL_LOOSE", [(full_loose_complex,3,2)])
    

    o = Options(simulation_mode="Transition",parameter_type="Nupack",substrate_type="DNA",
                num_sims = 1, sim_time=.01, start_state=[start_complex])

    o.stop_conditions = [initial_sc, pathway_sc, pathway_loose_sc,full_loose_sc]
    o.temperature = 310.15
    return o


def setup_options_hp_final():
    """
    setup_options_hp_final( )

    returns the options object for our simple hairpin example for
    transition mode.
    """

    stem = Domain(name="stem",sequence="GCATGC",length=6)
    hp = Domain(name="hairpin", sequence="AAAA",length=4)
    
    s = stem + hp + stem.C
    
    start_complex = Complex(strands=[s], structure="...")

    pathway_endside_complex = Complex(strands=[s], structure=      "(((..........)))")
    pathway_endside_loose_complex = Complex(strands=[s], structure="(((**********)))")

    pathway_hpside_complex = Complex(strands=[s], structure="...(((....)))...")
    pathway_hpside_loose_complex = Complex(strands=[s], structure="***(((****)))***")



    full_complex = Complex( strands=[s], structure="(.)")
#    full_loose_complex = Complex( strands=[s], structure="((((((....))))))")

    initial_sc = StopCondition( "INITIAL", [(start_complex,0,0)])
    pathway_hp_sc = StopCondition( "HPSIDE", [(pathway_hpside_complex,0,0)])
    pathway_hp_loose_sc = StopCondition( "HPSIDE_LOOSE", [(pathway_hpside_loose_complex,3,2)])

    pathway_end_sc = StopCondition( "ENDSIDE", [(pathway_endside_complex,0,0)])
    pathway_end_loose_sc = StopCondition( "ENDSIDE_LOOSE", [(pathway_endside_loose_complex,3,2)])
    
    full_sc    = StopCondition( "stop:FULL", [(full_complex,0,0)])
#    full_loose_sc = StopCondition( "stop:FULL_LOOSE", [(full_loose_complex,3,2)])
    

    o = Options(simulation_mode="Transition",parameter_type="Nupack",substrate_type="DNA",
                num_sims = 1, sim_time=.01, start_state=[start_complex])

    o.stop_conditions = [initial_sc, pathway_end_sc, pathway_hp_sc,full_sc]
    o.temperature = 310.15

    oloose = Options(simulation_mode="Transition",parameter_type="Nupack",substrate_type="DNA",
                num_sims = 1, sim_time=.01, start_state=[start_complex])

    oloose.stop_conditions = [initial_sc, pathway_end_loose_sc, pathway_hp_loose_sc,full_sc]
    oloose.temperature = 310.15
    return o,oloose

def create_setup():
    toehold_seq = "GTGGGT"
    bm_design = "ACCGCACGTCACTCACCTCG"


    toehold = Domain(name="toehold",sequence=toehold_seq,length=6)
    branch_migration = Domain(name="bm", sequence=bm_design,
                                seq_length=20)
    
    substrate = toehold + branch_migration
    incumbent = Strand(name="incumbent",domains=[branch_migration.C])

    incoming = substrate.C

    start_complex = Complex(strands=[incoming, substrate, incumbent],
                              structure=".(+)(+)")
    stop_conditions = []
    stop_conditions.append( StopCondition("INITIAL", [(start_complex,4,2)]))
                            # Within distance 2 of the start_complex state.

    initial_structure_dp =  "....................((((((+))))))((((((((((((((((((((+))))))))))))))))))))"
    six_bases_structure_dp= "..............((((((((((((+))))))))))))((((((((((((((+))))))))))))))......"
    six_bases_loose_dp=     "**************((**********+**********))((************+************))******"
    twelve_bases_struc_dp  ="........((((((((((((((((((+))))))))))))))))))((((((((+))))))))............"
    twelve_bases_loose_dp=  "********((*****************+***************))((******+******))************"
    eighteen_structure_dp = "..((((((((((((((((((((((((+))))))))))))))))))))))))((+)).................."
    eighteen_loose_dp =     "**((**********************+**********************))((+))******************"

    six_bases_complex = Complex(strands=[incoming,substrate,incumbent],
                                structure=six_bases_structure_dp)
    twelve_bases_complex = Complex(strands=[incoming,substrate,incumbent],
                                structure=twelve_bases_struc_dp)
    eighteen_bases_complex = Complex(strands=[incoming,substrate,incumbent],
                                structure=eighteen_structure_dp)
    six_basesloose_complex = Complex(strands=[incoming,substrate,incumbent],
                                structure=six_bases_loose_dp)
    twelve_basesloose_complex = Complex(strands=[incoming,substrate,incumbent],
                                structure=twelve_bases_loose_dp)
    eighteen_basesloose_complex = Complex(strands=[incoming,substrate,incumbent],
                                structure=eighteen_loose_dp)

    disassoc_complex = Complex(strands=[incumbent], structure=".")

    stop_conditions.append( StopCondition("SIX", [(six_bases_complex,0,0)]))
#    stop_conditions.append( StopCondition("SIX_LOOSE", [(six_basesloose_complex,3,2)]))
    stop_conditions.append( StopCondition("TWELVE", [(twelve_bases_complex,0,0)]))
#    stop_conditions.append( StopCondition("TWELVE_LOOSE", [(twelve_basesloose_complex,3,2)]))

    stop_conditions.append( StopCondition("EIGHTEEN", [(eighteen_bases_complex,0,0)]))
#    stop_conditions.append( StopCondition("EIGHTEEN_LOOSE", [(eighteen_basesloose_complex,3,2)]))
    stop_conditions.append( StopCondition("stop:COMPLETE", [(disassoc_complex,2,0)]))


    o = Options(simulation_mode="Transition",parameter_type="Nupack",
                substrate_type="DNA",num_sims = 1,sim_time=.01,
                start_state=[start_complex])

    o.stop_conditions = stop_conditions
    o.temperature = 310.15
    o.dangles = 1
    
    return o

    # bm_design_A = "ACCGCACGTCCACGGTGTCG"
    # branch_migration_A = Domain(name="bm_A", sequence=bm_design_A,
    #                             seq_length=20)
    # substrate_A = toehold + branch_migration_A
    # incumbent_A = Strand(name="incumbent",domains=[branch_migration_A.C])
    # incoming_A = substrate_A.C
    # start_complex_A = Complex(strands=[incoming_A, substrate_A, incumbent_A],
    #                           structure=".(+)(+)")




def parse_transition_list( transition_traj ):
    transition_dict = {}
    
    def in_state( l ):
        flag = False
        for i in l:
            flag = flag or i
        return flag

    def name(t0,t1):
        charindex = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0'
        names0 = [charindex[j] for i,j in zip(t0,range(len(t0))) if i]
        if names0 == []:
            names0 = charindex[26]
        else:
            names0 = ",".join(names0)
        names1 = [charindex[j] for i,j in zip(t1,range(len(t1))) if i]
        if names1 == []:
            names1 = charindex[26]
        else:
            names1 = ",".join(names1)
        return names0 + " -> " + names1
        
        # try:
        #     c0 = t0.index(True)
        # except ValueError:
        #     c0 = 26
        # try:
        #     c1 = t1.index(True)
        # except ValueError:
        #     c1 = 26
            
        # return charindex[c0] + " -> " + charindex[c1]

    truncated = [i for i in transition_traj if in_state(i[1])]
    tt = truncated # or transition_traj, but that's boring...

    for i in range(len(tt)-1):
        nm = name(tt[i][1],tt[i+1][1])
        if nm in transition_dict:
            transition_dict[nm].append( tt[i+1][0] - tt[i][0] )
        else:
            transition_dict[nm] = [tt[i+1][0] - tt[i][0]]


    return truncated, transition_dict

    
def print_transition_dict( transition_dict, options = None ):
    k = transition_dict.keys();
    k.sort() # Used to be k.sort(reverse=True)  due to T/F strings

    for i in k:
        transition_times = np.array( transition_dict[i] )
        
        print("{0}: {2:.2e} ({1})".format(i,len(transition_dict[i]),np.mean(transition_times)))
    
    charindex = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0'
    if options:
        for i,idx in zip(options.stop_conditions,range(len(options.stop_conditions))):
            print("{0}: {1}".format( i.tag, charindex[idx]))
                  
            


def parse_transition_lists( transition_traj_list ):
    transition_dict = {}
    
    def in_state( l ):
        flag = False
        for i in l:
            flag = flag or i
        return flag

    def name(t0,t1):
        charindex = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0'
        names0 = [charindex[j] for i,j in zip(t0,range(len(t0))) if i]
        if names0 == []:
            names0 = charindex[26]
        else:
            names0 = ",".join(names0)
        names1 = [charindex[j] for i,j in zip(t1,range(len(t1))) if i]
        if names1 == []:
            names1 = charindex[26]
        else:
            names1 = ",".join(names1)
        return names0 + " -> " + names1
    for transition_traj in transition_traj_list:
        truncated = [i for i in transition_traj if in_state(i[1])]
        tt = truncated # or transition_traj, but that's boring...

        for i in range(len(tt)-1):
            nm = name(tt[i][1],tt[i+1][1])
            if nm in transition_dict:
                transition_dict[nm].append( tt[i+1][0] - tt[i][0] )
            else:
                transition_dict[nm] = [tt[i+1][0] - tt[i][0]]


    return truncated, transition_dict
