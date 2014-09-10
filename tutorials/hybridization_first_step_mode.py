# hybridization_first_step_mode.py
#
# adapted from three_way_sims_example.py by Niranjan Srinivas and Joseph Schaeffer.
#
# Usage:
# python -i hybridization_first_step_mode.py


# Import things you need
import cPickle
import random
import numpy as np

from multistrand_setup import *

#########

# for StopCondition macrostate definitions:
Exact_Macrostate = 0
Bound_Macrostate = 1
Dissoc_Macrostate = 2
Loose_Macrostate = 3
Count_Macrostate = 4

def create_setup(strand_seq, num_traj, T=25, rate_method_k_or_m="Metropolis"):
    # essentially, creates the options object and prepares to simulate

    onedomain = Domain(name="itall",sequence=strand_seq)
    top = Strand(name="top",domains=[onedomain])
    bot = top.C

    start_complex_top = Complex(strands=[top],structure=".")
    start_complex_bot = Complex(strands=[bot],structure=".")
    start_complex_top.boltzmann_count = num_traj
    start_complex_bot.boltzmann_count = num_traj
    start_complex_top.boltzmann_sample = True
    start_complex_bot.boltzmann_sample = True
    # turns Boltzmann sampling on for this complex and also does sampling more efficiently by sampling #num_traj states

    success_complex = Complex(strands=[top, bot],structure="(+)")
    success_stop_condition = StopCondition("SUCCESS",[(success_complex,Exact_Macrostate,0)])

    failed_complex = Complex(strands = [top], structure=".")
    failed_stop_condition = StopCondition("FAILURE",[(failed_complex,Dissoc_Macrostate,0)])

    o = Options(simulation_mode="First Step",parameter_type="Nupack", substrate_type="DNA",
                rate_method = rate_method_k_or_m, num_simulations = num_traj, simulation_time=1.0,
                dangles = "Some", temperature = T, rate_scaling = "Calibrated", verbosity = 0)

    o.start_state = [start_complex_top, start_complex_bot]
    o.stop_conditions = [success_stop_condition,failed_stop_condition]
    return o

def compute_rate_constants(forward_times, reverse_times, collision_rates, concentration):
    N_forward = len(forward_times)
    N_reverse = len(reverse_times)
    N = N_forward+N_reverse
    dTsuccess_uni = np.mean(forward_times)
    k2 = 1.0/dTsuccess_uni
    dTfail_uni   = np.mean(reverse_times)
    k2prime = 1.0/dTfail_uni
    kcollision = np.mean(collision_rates)
    prob = float(N_forward)/N
    k1 = prob * kcollision
    k1prime = (1-prob) * kcollision

    z = concentration
    dTcoll = 1/((k1+k1prime)*z)
    dTfail = dTcoll + dTfail_uni
    dTforward = dTcoll + dTsuccess_uni
    dTcorrect = dTfail*k1prime/k1 + dTforward
    keff = (1/dTcorrect)*(1/z)

    return N_forward, N_reverse, k1, k2, k1prime, k2prime, keff


def first_step_simulation(strand_seq, num_traj, T=25, rate_method_k_or_m="Metropolis", concentration=50e-9):

    # run the simulations

    o = create_setup(strand_seq, num_traj, T, rate_method_k_or_m)
    s = SimSystem(o)
    s.start()
    dataset = o.interface.results

    # extract the timing information for successful and failed runs

    forward = [i for i in dataset if i.tag == "SUCCESS"]
    forward_times = np.zeros( len(forward))
    forward_times[:] = [i.time for i in forward]

    reverse = [i for i in dataset if i.tag == "FAILURE" or i.tag == None]
    reverse_times = np.zeros( len(reverse))
    reverse_times[:] = [i.time for i in reverse]

    collision_rates = np.zeros( len(dataset))
    collision_rates[:] = [i.collision_rate for i in dataset]

    # basic calculations from Joseph's PhD thesis.  error bar analysis is new here, not checked whether it's the same. (SHOULD CHECK.)
    # note: will fail if there are either no successful trials or no failed trials

    N_forward, N_reverse, k1, k2, k1prime, k2prime, keff = compute_rate_constants(forward_times, reverse_times, collision_rates, concentration)

    # inset bootstrapping calculations here, to replace the stuff below

    N_forward = len(forward_times)
    N_reverse = len(reverse_times)
    N = N_forward+N_reverse
    dTsuccess_uni = np.mean(forward_times)
    k2 = 1.0/dTsuccess_uni
    std_k2 = k2 * np.std(forward_times)/np.sqrt(N_forward)/np.mean(forward_times) # linear approx: same % error in times as in rates
    dTfail_uni   = np.mean(reverse_times)
    k2prime = 1.0/dTfail_uni
    std_k2prime = k2prime * np.std(reverse_times)/np.sqrt(N_reverse)/np.mean(reverse_times) # linear approx: same % error in times as in rates
    kcollision = np.mean(collision_rates)
    std_kcollision = np.std(collision_rates) / np.sqrt(N)
    prob = float(N_forward)/N
    k1 = prob * kcollision
    std_k1 = np.sqrt(prob*(1-prob)) * kcollision / np.sqrt(N)  # does not account for kcollision variability
    k1prime = (1-prob) * kcollision
    std_k1prime = np.sqrt(prob*(1-prob)) * kcollision / np.sqrt(N) # does not account for kcollision variability

    z = concentration

    dTcoll = 1/((k1+k1prime)*z)
    dTfail = dTcoll + dTfail_uni
    dTforward = dTcoll + dTsuccess_uni
    dTcorrect = dTfail*k1prime/k1 + dTforward

    keff = (1/dTcorrect)*(1/z)
    # would be better to provide bootstrapped error bars by dividing the data set into 4, computing the std of the 4 resulting estimates, and dividing by 2.

    # print out the results

    print "N_forward =", N_forward
    print "N_reverse =", N_reverse
    print "k_eff =", keff, "/M/s at 50 nM", "(still needs error bars)"
    print "k1 =", k1, "+/-", std_k1, "/M/s", "(i.e. +/-", 100*std_k1/k1, "%)"
    print "k2 =", k2, "+/-", std_k2, "/s", "(i.e. +/-", 100*std_k2/k2, "%)"
    print "k1prime =", k1prime, "+/-", std_k1prime, "/M/s", "(i.e. +/-", 100*std_k1prime/k1prime, "%)"
    print "k2prime =", k2prime, "+/-", std_k2prime, "/s", "(i.e. +/-", 100*std_k2prime/k2prime, "%)"
    print "k_collision =", kcollision, "+/-", std_kcollision, "/M/s", "(i.e. +/-", 100*std_kcollision/kcollision, "%)"

    return [N_forward, N_reverse, keff, k1, k1prime, k2, k2prime, kcollision, o]    # add error bars here.  change to dict.


if __name__ == '__main__':
    data=first_step_simulation("ACTGGCGCGTATTATCTACTG", 100, concentration=50e-9)
