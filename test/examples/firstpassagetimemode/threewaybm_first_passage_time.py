import sys,os,os.path
import cPickle

import numpy as np
import matplotlib
matplotlib.use('macosx')
import matplotlib.pylab as pyp


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

    # import multiprocessing
    # from multiprocessing import Pool
except ImportError:
    print("Could not import Multistrand.")
    raise

def create_setup():

    toehold_seq = "GTGGGT"

    bm_design_A = "ACCGCACGTCCACGGTGTCG"

    bm_design_B = "ACCGCACGTCACTCACCTCG"


    toehold = Domain(name="toehold",sequence=toehold_seq,length=6)
    branch_migration_A = Domain(name="bm_A", sequence=bm_design_A,
                                seq_length=20)

    branch_migration_B = Domain(name="bm_B", sequence=bm_design_B,
                                seq_length=20)
    
    substrate_A = toehold + branch_migration_A
    substrate_B = toehold + branch_migration_B
    incumbent_A = Strand(name="incumbent",domains=[branch_migration_A.C])
    incumbent_B = Strand(name="incumbent",domains=[branch_migration_B.C])

    incoming_A = substrate_A.C
    incoming_B = substrate_B.C

    start_complex_A = Complex(strands=[incoming_A, substrate_A, incumbent_A],
                              structure=".(+)(+)")
    start_complex_B = Complex(strands=[incoming_B, substrate_B, incumbent_B],
                              structure=".(+)(+)")
    
    complete_complex_A = Complex(strands=[incumbent_A],structure=".")
    complete_complex_B = Complex(strands=[incumbent_B],structure=".")

    failed_complex_A = Complex(strands=[incoming_A],structure="..")
    failed_complex_B = Complex(strands=[incoming_B],structure="..")

    complete_sc_A = StopCondition("complete",[(complete_complex_A,2,0)])

    complete_sc_B = StopCondition("complete",[(complete_complex_B,2,0)])

    failed_sc_A = StopCondition("failed",[(failed_complex_A,2,0)])

    failed_sc_B = StopCondition("failed",[(failed_complex_B,2,0)])


    
    o = Options(simulation_mode="First Passage Time",parameter_type="Nupack",
                substrate_type="DNA",num_sims = 100,sim_time=1.0,
                dangles = "Some", temperature = 310.15,
                start_state=[start_complex_A])

    o.stop_conditions = [complete_sc_A,failed_sc_A]

    o2 = Options(simulation_mode="First Passage Time",parameter_type="Nupack",
                substrate_type="DNA",num_sims = 100,sim_time=1.0,
                dangles = "Some", temperature = 310.15,
                start_state=[start_complex_B])

    o2.stop_conditions = [complete_sc_B,failed_sc_B]
    
    return o,o2

def create_setup_long_toehold():

    toehold_seq = "GTGGGTAGGT"

    bm_design_A = "ACCGCACGTCCACGGTGTCG"

    bm_design_B = "ACCGCACGTCACTCACCTCG"


    toehold = Domain(name="toehold",sequence=toehold_seq,length=10)
    branch_migration_A = Domain(name="bm_A", sequence=bm_design_A,
                                seq_length=20)

    branch_migration_B = Domain(name="bm_B", sequence=bm_design_B,
                                seq_length=20)
    
    substrate_A = toehold + branch_migration_A
    substrate_B = toehold + branch_migration_B
    incumbent_A = Strand(name="incumbent",domains=[branch_migration_A.C])
    incumbent_B = Strand(name="incumbent",domains=[branch_migration_B.C])

    incoming_A = substrate_A.C
    incoming_B = substrate_B.C

    start_complex_A = Complex(strands=[incoming_A, substrate_A, incumbent_A],
                              structure=".(+)(+)")
    start_complex_B = Complex(strands=[incoming_B, substrate_B, incumbent_B],
                              structure=".(+)(+)")
    
    complete_complex_A = Complex(strands=[incumbent_A],structure=".")
    complete_complex_B = Complex(strands=[incumbent_B],structure=".")

    failed_complex_A = Complex(strands=[incoming_A],structure="..")
    failed_complex_B = Complex(strands=[incoming_B],structure="..")

    complete_sc_A = StopCondition("complete",[(complete_complex_A,2,0)])

    complete_sc_B = StopCondition("complete",[(complete_complex_B,2,0)])

    failed_sc_A = StopCondition("failed",[(failed_complex_A,2,0)])

    failed_sc_B = StopCondition("failed",[(failed_complex_B,2,0)])


    
    o = Options(simulation_mode="First Passage Time",parameter_type="Nupack",
                substrate_type="DNA",num_sims = 100,sim_time=1.0,
                dangles = "Some", temperature = 310.15,
                start_state=[start_complex_A])

    o.stop_conditions = [complete_sc_A,failed_sc_A]

    o2 = Options(simulation_mode="First Passage Time",parameter_type="Nupack",
                substrate_type="DNA",num_sims = 100,sim_time=1.0,
                dangles = "Some", temperature = 310.15,
                start_state=[start_complex_B])

    o2.stop_conditions = [complete_sc_B,failed_sc_B]
    
    return o,o2


def plot_histograms_complete_vs_failed( result_list, colors=['b','m'], figure=1 ):
    times_complete = np.array([i.time for i in result_list if i.tag == 'complete'])
    times_failed = np.array([i.time for i in result_list if i.tag == 'failed'])

    min_time = np.min( [np.min(times_complete),np.min(times_failed)] )
    max_time = np.max( [np.max(times_complete),np.max(times_failed)] )

    pyp.figure( figure )

    pyp.hist( times_complete, 50, range=(min_time,max_time), color=colors[0], label="complete" )

    pyp.hold(True)
    
    pyp.hist( times_failed, 50, range=(min_time,max_time), color=colors[1],rwidth=.5, label="failed")

    pyp.xlabel("First Passage Time (s)",fontsize='larger')
    pyp.ylabel("# of Trajectories",fontsize='larger')
    pyp.yticks(fontsize='larger',va='bottom')
    pyp.xticks(fontsize='larger')
    pyp.legend(loc=0)

def plot_histograms_two_designs( result_lists, colors=['magenta','b'], figure=1 ):
    times = []
    times.append(np.array([i.time for i in result_lists[0] if i.tag == 'complete']))
    times.append(np.array([i.time for i in result_lists[1] if i.tag == 'complete']))

    min_time = np.min( [np.min(times[0]),np.min(times[1])] )
    max_time = np.max( [np.max(times[0]),np.max(times[1])] )
    max_time = .02

    pyp.figure( figure )
    pyp.hold(False)

    pyp.hist( times[1], 50, range=(min_time,max_time), color=colors[1],label="Design B")
    # pyp.hist( times[1], 20, range=(min_time,max_time), color=colors[1],rwidth=.5, label="Design B") 


    pyp.hold(True)
    
    pyp.hist( times[0], 50, range=(min_time,max_time), color=colors[0], label="Design A", rwidth=.5 )

    pyp.xlabel("First Passage Time (s)",fontsize='larger')
    pyp.ylabel("# of Trajectories",fontsize='larger')
    pyp.yticks(fontsize='larger',va='bottom')
    pyp.xticks(fontsize='larger')
    pyp.legend(loc=0)



def plot_histograms_single( result_list, colors=['r','g'], figure=1, label="" ):
    times = []
    times.append(np.array([i.time for i in result_list if i.tag == 'complete']))

    pyp.figure( figure )

    if label != "":
        pyp.hist( times[0], 50, color=colors[0], label=label )
    else:
        pyp.hist( times[0], 50, color=colors[0] )
        

    # pyp.hold(True)
    
    # pyp.hist( times[1], 200, range=(min_time,max_time), color=colors[1],rwidth=.5, label="Design B")

    pyp.xlabel("First Passage Time (s)", fontsize='larger')
    pyp.ylabel("# of Trajectories",fontsize='larger')
    pyp.yticks(fontsize='larger',va='bottom')
    pyp.xticks(fontsize='larger')

    pyp.legend(loc=0)


def plot_completion_graph( result_lists, colors=['b','r'], figure=1, labels=[] ):
    times = []
    percents = []
    for rl in result_lists:
        t = [i.time for i in rl if i.tag == 'complete']
        t.sort()
        t.append(1.0)
        times.append(np.array(t))

        p = np.array(range(1,len(t)+1))
        p = p / 10.0
        p[-1] = p[-2]
        percents.append( p )

    pyp.figure( figure )

    for t,p,c in zip(times,percents,colors):
        pyp.plot( t, p, color = c, linewidth=2.0 )
        pyp.hold(True)

    maxtime = np.max([i[-1] for i in times])
    
    pyp.xlabel("Simulation Time (s)",fontsize='larger')
    pyp.ylabel("% of Trajectories Complete",fontsize='larger')
    pyp.yticks([0,20,40,60,80,100],("0%","20%","40%","60%","80%","100%"),fontsize='larger',va='bottom')
    pyp.xticks(fontsize='larger')

#    pyp.title( "Percentage of Total Trajectories Complete by 
    pyp.legend(loc=0)


def plot_completion_graph_complete_vs_failed( result_list, colors=['b','m'], figure=1, labels=[] ):
    times = []
    percents = []
    for rl in ['complete','failed']:
        t = [i.time for i in result_list if i.tag == rl]
        t.sort()
        t.append(1.0)
        times.append(np.array(t))

        p = np.array(range(1,len(t)+1))
        p = p / 10.0
        p[-1] = p[-2]
        percents.append( p )

    pyp.figure( figure )

    for t,p,c,label in zip(times,percents,colors,['complete','failed']):
        pyp.plot( t, p, color = c, linewidth=2.0, label=label )
        pyp.hold(True)

    maxtime = np.max([i[-1] for i in times])
    
    pyp.xlabel("Simulation Time (s)",fontsize='larger')
    pyp.ylabel("% of Trajectories Complete",fontsize='larger')
    pyp.yticks([0,20,40,60,80,100],("0%","20%","40%","60%","80%","100%"),fontsize='larger',va='bottom')
    pyp.xticks(fontsize='larger')

#    pyp.title( "Percentage of Total Trajectories Complete by 
    pyp.legend(loc=0)

        
