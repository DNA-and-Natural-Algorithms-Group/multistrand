# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
First Passage Time can be used to get statistics on how long it takes to reach a
certain state, or one of a set of defined states. Here we compare two sequences
for toehold-mediated three-way strand displacement. Starting in with the
incoming strand binding by the toehold only, we run simulations until either the
incoming strand falls off again (a "failure") or the incumbent top strand is
displaced and dissociates (a "success"). This is specified via Literals.Literals.dissoc_macrostate
StopConditions that look for the incoming strand being in a complex by itself,
or the incumbent strand being in a complex by itself, respectively. We perform
only 100 trials for each sequence, because that already takes around 10 minutes
(most of which is design A). Once that's done, for each sequence we have 100
trial results, each tagged with either "failure" or "success". We can look at
histograms of how long such trajectories took, or alternatively tabulate what
fraction of trajectories completed by a certain time.

Try it like this, e.g.:
  ipython --pylab -i threewaybm_first_passage_time.py
"""

from multistrand.objects import *
from multistrand.options import Options, Literals
from multistrand.system import SimSystem

import numpy as np
import matplotlib
import matplotlib.pylab as plt

######### This is first passage time

# here, design A has a hairpin and is slow / fails often.  design B is fast
# in figure 7.4, design A has the same sequence, but...                                                     design B is "ACCGCACCACGTGGGTGTCG"
def setup_threewaybm_comparison(trials=100, toehold_seq="GTGGGT", bm_design_A="ACCGCACGTCCACGGTGTCG",
                                bm_design_B="ACCGCACGTCACTCACCTCG"):
    toehold = Domain(name="toehold", sequence=toehold_seq, length=len(toehold_seq))
    branch_migration_A = Domain(name="bm_A", sequence=bm_design_A, seq_length=len(bm_design_A))
    branch_migration_B = Domain(name="bm_B", sequence=bm_design_B, seq_length=len(bm_design_B))

    substrate_A = toehold + branch_migration_A
    substrate_B = toehold + branch_migration_B
    incumbent_A = Strand(name="incumbent", domains=[branch_migration_A.C])
    incumbent_B = Strand(name="incumbent", domains=[branch_migration_B.C])

    incoming_A = substrate_A.C
    incoming_B = substrate_B.C

    # start with the full toehold bound
    start_complex_A = Complex(strands=[incoming_A, substrate_A, incumbent_A], structure=".(+)(+)")
    start_complex_B = Complex(strands=[incoming_B, substrate_B, incumbent_B], structure=".(+)(+)")

    complete_complex_A = Complex(strands=[incumbent_A], structure=".")
    complete_complex_B = Complex(strands=[incumbent_B], structure=".")

    failed_complex_A = Complex(strands=[incoming_A], structure="..")
    failed_complex_B = Complex(strands=[incoming_B], structure="..")

    # the "success" state
    complete_sc_A = StopCondition("complete", [(complete_complex_A, Literals.dissoc_macrostate, 0)])
    complete_sc_B = StopCondition("complete", [(complete_complex_B, Literals.dissoc_macrostate, 0)])

    # the "failure" state
    failed_sc_A = StopCondition("failed", [(failed_complex_A, Literals.dissoc_macrostate, 0)])
    failed_sc_B = StopCondition("failed", [(failed_complex_B, Literals.dissoc_macrostate, 0)])

    o1 = Options(simulation_mode="First Passage Time", num_simulations=trials, simulation_time=1.0,
                 dangles="Some", temperature=310.15, join_concentration=1e-6,  # 1 uM concentration
                 start_state=[start_complex_A], verbose=1)
    o1.stop_conditions = [complete_sc_A, failed_sc_A]
    o1.DNA23Metropolis()

    o2 = Options(simulation_mode="First Passage Time", num_simulations=trials, simulation_time=1.0,
                 dangles="Some", temperature=310.15, join_concentration=1e-6,  # 1 uM concentration
                 start_state=[start_complex_B], verbose=1)
    o2.stop_conditions = [complete_sc_B, failed_sc_B]
    o2.DNA23Metropolis()

    return o1, o2


def plot_histograms_complete_vs_failed(result_list, colors=['b', 'm'], figure=1):
    # separate based on which stop condition ended the simulation   ((( what if simulation ran out of time? )))
    times_complete = np.array([i.time for i in result_list if i.tag == 'complete'])
    times_failed = np.array([i.time for i in result_list if i.tag == 'failed'])
    neither = [i for i in result_list if not i.tag == 'complete' and not i.tag == 'failed']
    if len(neither) > 0:
        print("some trajectories did not finish, one way nor the other...")
        for i in neither:
            assert (i.type_name == 'Time')
            assert (i.tag == None)

    if len(times_complete) > 0 and len(times_failed) > 0:
        min_time = np.min([np.min(times_complete), np.min(times_failed)])
        max_time = np.max([np.max(times_complete), np.max(times_failed)])
    else:
        if len(times_complete) > 0:
            min_time = np.min(times_complete)
            max_time = np.max(times_complete)
        if len(times_failed) > 0:
            min_time = np.min(times_failed)
            max_time = np.max(times_failed)

    plt.figure(figure)

    if len(times_complete) > 0:
        plt.hist(times_complete, 50, range=(min_time, max_time), color=colors[0], label="complete")

    if len(times_failed) > 0:
        plt.hist(times_failed, 50, range=(min_time, max_time), color=colors[1], rwidth=.5, label="failed")

    plt.title("Completion times for successful and failed trajectories, Design A")
    plt.xlabel("First Passage Time (s)", fontsize='larger')
    plt.ylabel("# of Trajectories", fontsize='larger')
    plt.yticks(fontsize='larger', va='bottom')
    plt.xticks(fontsize='larger')
    plt.legend(loc=0)
    plt.show()


def plot_histograms_two_designs( result_lists, colors=['magenta','b'], figure=1 ):
    times = []
    times.append(np.array([i.time for i in result_lists[0] if i.tag == 'complete']))
    times.append(np.array([i.time for i in result_lists[1] if i.tag == 'complete']))

    min_time = np.min( [np.min(times[0]),np.min(times[1])] )
    max_time = np.max( [np.max(times[0]),np.max(times[1])] )
    # the above might fail if any list is empty; min, max of empty is undefined

    plt.figure( figure )
    plt.hist( times[1], 50, range=(min_time,max_time), color=colors[1], label="Design B")
    plt.hist( times[0], 50, range=(min_time,max_time), color=colors[0], label="Design A", rwidth=.5 )

    plt.title("Successful trajectories in two designs")
    plt.xlabel("First Passage Time (s)",fontsize='larger')
    plt.ylabel("# of Trajectories",fontsize='larger')
    plt.yticks(fontsize='larger',va='bottom')
    plt.xticks(fontsize='larger')
    plt.legend(loc=0)
    plt.show()


def plot_completion_graph( result_lists, colors=['b','r'], figure=1, labels=['Design A', 'Design B'] ):
    times = []
    percents = []
    for rl in result_lists:
        n = len(rl)
        t = [i.time for i in rl if i.tag == 'complete']
        sorted(t)
        times.append(np.array(t))

        p = np.array(range(1,len(t)+1))
        p = 100 * p / n  # percentage of all trials, including ones that don't complete
        percents.append( p )

    plt.figure( figure )

    for t,p,c,label in zip(times,percents,colors,labels):
        plt.plot( t, p, color = c, linewidth=2.0, label=label )

    plt.xlabel("Simulation Time (s)",fontsize='larger')
    plt.ylabel("% of Trajectories Complete",fontsize='larger')
    plt.yticks([0,20,40,60,80,100],("0%","20%","40%","60%","80%","100%"),fontsize='larger',va='bottom')
    plt.xticks(fontsize='larger')
    plt.title( "Percentage of Total Trajectories Complete by a Given Time" )
    plt.legend(loc=0)
    plt.show()


def plot_completion_graph_complete_vs_failed( result_list, colors=['b','m'], figure=1, labels=['complete','failed'] ):
    n = len(result_list) # number of trials
    times = []
    percents = []
    for rl in ['complete','failed']:
        t = [i.time for i in result_list if i.tag == rl]
        t.sort()
        times.append(np.array(t))
        p = np.array(range(1,len(t)+1))
        p = 100 * p / n
        percents.append( p )
    # the last data points of each type should sum to 1.0 if all trajectories either fail or complete (rather than time out)

    plt.figure( figure )

    for t,p,c,label in zip(times,percents,colors,labels):
        plt.plot( t, p, color = c, linewidth=2.0, label=label )

    plt.xlabel("Simulation Time (s)",fontsize='larger')
    plt.ylabel("% of Trajectories Complete",fontsize='larger')
    plt.yticks([0,20,40,60,80,100],("0%","20%","40%","60%","80%","100%"),fontsize='larger',va='bottom')
    plt.xticks(fontsize='larger')
    plt.title( "Percentage of Total Trajectories Complete by a Given Time, Design A" )
    plt.legend(loc=0)
    plt.show()


# Takes about 10 min for 100 trials.
if __name__ == '__main__':

    # o1,o2 = setup_threewaybm_comparison()  # equivalent to below, since they are defaults
    o1, o2 = setup_threewaybm_comparison(trials=100, toehold_seq="GTGGGT", bm_design_A="ACCGCACGTCCACGGTGTCG",
                                         bm_design_B="ACCGCACGTCACTCACCTCG")
    # o1,o2 = setup_threewaybm_comparison(trials=100, toehold_seq = "GTGGGTAGGT", bm_design_A = "ACCGCACGTCCACGGTGTCG", bm_design_B = "ACCGCACGTCACTCACCTCG")

    print()
    print("Running Design A")
    s = SimSystem(o1)
    s.start()
    print()
    print("Running Design B")
    s = SimSystem(o2)
    s.start()
    # no need to update the energy model in between runs, because conditions are the same.

    result_list1 = o1.interface.results
    result_list2 = o2.interface.results
    # look at a few things by hand, just to check
    print(result_list1[0])

    times_complete1 = np.array([i.time for i in result_list1 if i.tag == 'complete'])
    times_failed1 = np.array([i.time for i in result_list1 if i.tag == 'failed'])
    print("Design A: %d trajectories total, %d completed, %d failed." % (
    len(result_list1), len(times_complete1), len(times_failed1)))

    times_complete2 = np.array([i.time for i in result_list2 if i.tag == 'complete'])
    times_failed2 = np.array([i.time for i in result_list2 if i.tag == 'failed'])
    print("Design B: %d trajectories total, %d completed, %d failed." % (
    len(result_list2), len(times_complete2), len(times_failed2)))

    plot_histograms_complete_vs_failed(o1.interface.results, figure=1)
    plot_completion_graph_complete_vs_failed(o1.interface.results, figure=2)
    plot_histograms_two_designs([o1.interface.results, o2.interface.results], figure=3)
    plot_completion_graph([o1.interface.results, o2.interface.results], figure=4)
