# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
Perform simulations that directly examine interactions between complexes,
deducing first-order and second-order rate constants for both productive
and non-productive interactions.  Here we examine the association of a
relatively unstructured DNA strand and its perfect complement.  (We also
simulate the structured RNA hybridization from Gao et al 2006, and structured
DNA hybridization as well.)

To extract rate constants, one could simulate two complexes in a volume of
the appropriate size for the concentration of interest, and simulate until the
reaction is complete.  We will actually do that in hybridization_comparison.py,
but it is very slow for low concentrations because most of the simulation time
is spent simulating intramolecular conformational changes prior to the
bimolecular interaction.  Multistrand's "first step mode" skips over those
pre-collision simulation steps by starting each simulation with a bimolecular
base-pair formation move -- i.e. a "collision" between the complexes.  To
approximate the pre-collision conformational changes, first step mode can
use NUPACK to Boltzmann sample the initial conformation -- i.e. assuming that
there is enough time between collisions for both complexes to equilibrate.
(Without Boltzmann sampling, as illustrated in threewaybm_first_step_mode.py,
every collision involves the exact same conformation of each complex, but
the first base pair formed between them will still be random.)

First step mode simulations proceed as usual after the first step; they are
typically set up with at least two StopConditions:  falling appart again into
the original complexes, or reaching some target state (e.g. a full duplex as in
in this example, or perhaps successful displacement of a strand as in the
three-way branch migration example).  From statistics on a large sample of
collision trajectories, including how long they take to reach some StopCondition,
the relevant rate constants can be computed.

In Schaeffer's PhD thesis, and here, we extract parameters for a model that
distinguishes between productive and unproductive reactions at the very onset.
For example, for hybridization, we have:
   A + B -->{k1} I -->{k2} C
and
   A + B -->{k1'} I' -->{k2'} A + B
where
   k1 is the second-order rate constant for successful collisions,
   k2 is the first-order rate constant for how long it takes to proceed
        from collision to reaching the final state,
   k1' is the second-order rate constant for unproductive collisions,
   k2' is the first-order rate constant for how long it remains in the
        unproductive intermediate configurations before dissociating.

Thus, k1 is the probability that a collision will be successful, times
the bimolecular collision rate averaged over all successful collisions.
In contrast, k2 is the inverse of the average completion time for
the successful collisions.  k1' and k2' are computed analogously.

In the limit of low concentrations, k1 is the only thing you care about.
However, at higher concentrations, the time spend "doing" the reaction
(i.e. k2) could be the rate limiting step, and thus must be taken into
account.  Similarly, at higher concentrations the time spent in
unproductive reactions could also be considerable, slowing down the overall
rate of completion of a reaction. k_eff is calculated from all four rate
constants *AND* the concentration to give a single "effective" second-order
rate constant that approximates the rate of production of product at this
concentration.  (Because k_eff is concentration-dependent, it must be
used with caution.)

More details are explained in the code, as well as methods for computing
error bars.

Usage:
    python -i hybridization_first_step_mode.py
Takes about 5 min.
"""

import random
import numpy as np

from multistrand.objects import *
from multistrand.options import Options, Literals
from multistrand.system import SimSystem


# for StopCondition macrostate definitions:
Exact_Macrostate = 0
Bound_Macrostate = 1
Dissoc_Macrostate = 2
Loose_Macrostate = 3
Count_Macrostate = 4

# from nupack import *
# print(str(get_nupack_exec_path('commandhere')))

def concentration_string(concentration):
    if concentration < 1e-12: 
        return "{} fM".format(1e15*concentration)
    if concentration < 1e-9: 
        return "{} pM".format(1e12*concentration)
    if concentration < 1e-6: 
        return "{} nM".format(1e9*concentration)
    if concentration < 1e-3: 
        return "{} uM".format(1e6*concentration)
    if concentration < 1: 
        return "{} mM".format(1e3*concentration)
    return "{} M".format(concentration)


def create_setup(strand_seq, num_traj, T=25, rate_method_k_or_m="Metropolis", material="DNA"):

    # Essentially, creates the options object and prepares to simulate the hybridization of the strand and its complement.
    onedomain = Domain(name="itall",sequence=strand_seq)
    top = Strand(name="top",domains=[onedomain])
    bot = top.C

    # Note that the structure is specified to be single stranded, but this will be over-ridden when Boltzmann sampling is turned on.
    start_complex_top = Complex(strands=[top],structure=".")
    start_complex_bot = Complex(strands=[bot],structure=".")
    start_complex_top.boltzmann_count = num_traj
    start_complex_bot.boltzmann_count = num_traj
    start_complex_top.boltzmann_sample = True
    start_complex_bot.boltzmann_sample = True
    # Turns Boltzmann sampling on for this complex and also does sampling more efficiently by sampling 'num_traj' states.

    # Stop when the exact full duplex is achieved. (No breathing!)
    success_complex = Complex(strands=[top, bot],structure="(+)")
    success_stop_condition = StopCondition("SUCCESS",[(success_complex,Exact_Macrostate,0)])

    # Declare the simulation unproductive if the strands become single-stranded again.
    failed_complex = Complex(strands = [top], structure=".")
    failed_stop_condition = StopCondition("FAILURE",[(failed_complex,Dissoc_Macrostate,0)])

    o = Options(simulation_mode="First Step",parameter_type="Nupack", substrate_type=material,
                rate_method = rate_method_k_or_m, num_simulations = num_traj, simulation_time=1.0,
                dangles = "Some", temperature = T, rate_scaling = "Calibrated", verbosity = 0)

    o.start_state = [start_complex_top, start_complex_bot]
    o.stop_conditions = [success_stop_condition,failed_stop_condition]
    return o


def compute_rate_constants(dataset, concentration, printit=True):

    # Basic calculations from Joseph's PhD thesis.  (See also threewaybm_first_step_mode.py.)
    # The equations there were derived for runs stating with exact microstates.
    # Here, as we are using Boltzmann sampling, which requires that the mean collision rate 
    # be calculated separately for forward and reverse cases, since conditioning on success
    # might alter the distribution of secondary structures, and thus the collision rates.

    # Also, the error bar analysis is new here.
    # Note: will fail if there are either no successful trials or no failed trials

    # Finally, note that the calculations below are fairly redundant; if you just want
    # k1, k2, k1prime, k2prime, then a bunch of stuff can be omitted and/or simplified.
    # hybridization_comparison.py gives more streamlined calculations (but no error bars).

    collision_rates = np.zeros( len(dataset))
    collision_rates[:] = [i.collision_rate for i in dataset]

    # Pull out the duration of successful reactions and 
    # the bimolecular rate constants for collisions between the particuular Boltzmann-sampled complexes for each trial.
    forward = [i for i in dataset if i.tag == "SUCCESS"]
    forward_times = np.zeros( len(forward))
    forward_times[:] = [i.time for i in forward]
    forward_collision_rates = np.zeros( len(forward))
    forward_collision_rates[:] = [i.collision_rate for i in forward]

    # When Boltzmann sampling, it is possible that the two complexes have no possible first base-pairs to form,
    # for example, if they are both blunt-ended hairpins.  In this case, i.collision_rate will be 0, i.tag will be 'noinitial'.  
    # Since the other way that i.tag can be 'None' is when a simulation doesn't
    # reach any StopCondition before timing out by exceeding o.simulation_time, and since both of those cases should be considered "failures"
    reverse = [i for i in dataset if i.tag == "FAILURE" or i.tag == Literals.no_initial_moves or i.tag == Literals.time_out]
    reverse_times = np.zeros( len(reverse))
    reverse_times[:] = [i.time for i in reverse]
    reverse_collision_rates = np.zeros( len(reverse))
    reverse_collision_rates[:] = [i.collision_rate for i in reverse]

    # How many simulation trials were successful, and how many were unproductive?
    N_forward = len(forward_times)
    N_reverse = len(reverse_times)
    N = N_forward + N_reverse

    # Calculate first-order rate constants for the duration of the reactions (both productive and unproductive).
    # The error bar formulas here are estimates, and they may not be accurate if the distributions of completion times are unusual or if there aren't enough trials.
    dTsuccess_uni = np.mean(forward_times)
    k2 = 1.0/dTsuccess_uni
    std_k2 = k2 * np.std(forward_times)/np.sqrt(N_forward)/np.mean(forward_times) # linear approx: same % error in times as in rates
    dTfail_uni   = np.mean(reverse_times)
    k2prime = 1.0/dTfail_uni
    std_k2prime = k2prime * np.std(reverse_times)/np.sqrt(N_reverse)/np.mean(reverse_times) # linear approx: same % error in times as in rates

    # Calculate second-order rate constants, and their error bars.
    kcollision = np.mean(collision_rates)
    reverse_kcoll = np.mean(reverse_collision_rates)
    forward_kcoll = np.mean(forward_collision_rates)
    std_kcollision = np.std(collision_rates) / np.sqrt(N)
    std_forward_kcoll = np.std(forward_collision_rates) / np.sqrt(N_forward)
    std_reverse_kcoll = np.std(reverse_collision_rates) / np.sqrt(N_reverse)
    prob = float(N_forward)/N
    k1 = prob * forward_kcoll    # this is mathematically equivalent to np.mean( collision_rates * was_success )  where * is pointwise, like Matlab .*
    std_k1 = np.std( np.concatenate([forward_collision_rates,np.zeros(N_reverse)]) ) / np.sqrt(N)
    # print "%g =?= %g" % ( k1 , np.mean( np.concatenate([forward_collision_rates,np.zeros(N_reverse)]) ) )   # prove the above claim
    k1prime = (1-prob) * reverse_kcoll
    std_k1prime = np.std( np.concatenate([reverse_collision_rates,np.zeros(N_forward)]) ) / np.sqrt(N)
    # print "%g =?= %g" % ( k1prime, np.mean( np.concatenate([reverse_collision_rates,np.zeros(N_forward)]) ) )

    # keff accounts both for potentially time-consuming unimolecular reactions and for potentially time-consuming "failed" interactions.
    z = concentration
    dTcoll = 1/((k1+k1prime)*z)                  # this is the average time for two single-stranded complexes to collide
    dTfail = dTcoll + dTfail_uni                 # conditioned on failure, the average time to collide and unproductively dissociate
    dTforward = dTcoll + dTsuccess_uni           # conditioned on success, the average time to collide and reach the duplex state
    dTcorrect = dTfail*k1prime/k1 + dTforward    # this is the average time for two single-stranded complexes to reach the duplex state after some failed attempts
    keff = (1/dTcorrect)*(1/z)
    # this is mathematically equivalent to  1/(1/k1 + z/k2 + (k1prime/k1)*(z/k2prime))
    # print "%g =?= %g" % ( keff, 1/(1/k1 + z/k2 + (k1prime/k1)*(z/k2prime)) )  # want me to prove it?
    zcrit = k2*k2prime/(k1*k2prime + k1prime*k2) # this is the critical concentration at which k_eff = k1/2

    # print out the results
    if printit:
        print("N_forward =", N_forward)
        print("N_reverse =", N_reverse)
        # print "k_collision = %g +/- %g /M/s (i.e. +/- %g %%)" % (kcollision,std_kcollision,100*std_kcollision/kcollision)
        print("k_collision_forward = %g +/- %g /M/s (i.e. +/- %g %%)" % (forward_kcoll,std_forward_kcoll,100*std_forward_kcoll/forward_kcoll))
        print("k_collision_reverse = %g +/- %g /M/s (i.e. +/- %g %%)" % (reverse_kcoll,std_reverse_kcoll,100*std_reverse_kcoll/reverse_kcoll))
        print("k1                  = %g +/- %g /M/s (i.e. +/- %g %%)" % (k1,std_k1,100*std_k1/k1))
        print("k2                  = %g +/- %g /s   (i.e. +/- %g %%)" % (k2,std_k2,100*std_k2/k2))
        print("k1prime             = %g +/- %g /M/s (i.e. +/- %g %%)" % (k1prime,std_k1prime,100*std_k1prime/k1prime))
        print("k2prime             = %g +/- %g /s   (i.e. +/- %g %%)" % (k2prime,std_k2prime,100*std_k2prime/k2prime))
        print("k_eff               = %g /M/s at %s (still needs error bars)" % (keff,concentration_string(concentration))) 
        print("z_crit              = %s (still needs error bars)" % (concentration_string(zcrit))) 

    return N_forward, N_reverse, kcollision, forward_kcoll, reverse_kcoll, k1, k2, k1prime, k2prime, keff, zcrit


def resample_with_replacement(mylist,num_samples):
    return [random.choice(mylist) for i in range(num_samples)]


def first_step_simulation(strand_seq, num_traj, T=25, rate_method_k_or_m="Metropolis", concentration=50e-9, material="DNA"):

    # Run the simulations

    print(f"Running first step mode simulations for {strand_seq} (with Boltzmann sampling)...")
    o = create_setup(strand_seq, num_traj, T, rate_method_k_or_m, material)
    s = SimSystem(o)
    s.start()
    dataset = o.interface.results

    # You might be interested in examining the data manually when num_traj < 10
    # for i in dataset:
    #    print i
    
    # Extract the timing information for successful and failed runs
    print()
    print("Inferred rate constants with analytical error bars:")
    N_forward, N_reverse, kcoll, forward_kcoll, reverse_kcoll, k1, k2, k1prime, k2prime, keff, zcrit = compute_rate_constants(dataset,concentration)

    # Bootstrapping is a technique that estimates statistical properties by assuming that the given samples adequately represent the true distribution,
    # and then resampling from that distribution to create as many mock data sets as you want.  The variation of statistical quantities 
    # in the mock data sets are often a good estimate of the true values.
    # We rely on bootstrapping to get error bars for k_eff, and to validate our estimated error bars for k2 and k2prime.

    Nfs, Nrs, kcfs, kcrs, k1s, k2s, k1primes, k2primes, keffs, zcrits = ([],[],[],[],[],[],[],[],[],[])
    for i in range(1000):
        t_dataset = resample_with_replacement(dataset,len(dataset))
        t_N_forward, t_N_reverse, t_kcoll, t_forward_kcoll, t_reverse_kcoll, t_k1, t_k2, t_k1prime, t_k2prime, t_keff, t_zcrit = \
            compute_rate_constants(t_dataset, concentration, printit=False)
        Nfs.append(t_N_forward)
        Nrs.append(t_N_reverse)
        kcfs.append(t_forward_kcoll)
        kcrs.append(t_reverse_kcoll)
        k1s.append(t_k1)
        k2s.append(t_k2)
        k1primes.append(t_k1prime)
        k2primes.append(t_k2prime)
        keffs.append(t_keff)
        zcrits.append(t_zcrit)

    std_Nfs = np.std(Nfs)
    std_Nrs = np.std(Nrs)
    std_kcfs = np.std(kcfs)
    std_kcrs = np.std(kcrs)
    std_k1 = np.std(k1s)
    std_k2 = np.std(k2s)
    std_k1prime = np.std(k1primes)
    std_k2prime = np.std(k2primes)
    std_keff = np.std(keffs)        
    std_zcrit = np.std(zcrits)        

    print()
    print("Re-sampled rate constants with bootstrapped error bars:")
    if True:
        print("N_forward = %d +/- %g" % (t_N_forward, std_Nfs))
        print("N_reverse = %d +/- %g" % (t_N_reverse, std_Nrs))
        print("k_collision_forward = %g +/- %g /M/s (i.e. +/- %g %%)" % (t_forward_kcoll, std_kcfs, 100*std_kcfs/forward_kcoll))
        print("k_collision_reverse = %g +/- %g /M/s (i.e. +/- %g %%)" % (t_reverse_kcoll, std_kcrs, 100*std_kcrs/reverse_kcoll))
        print("k1                  = %g +/- %g /M/s (i.e. +/- %g %%)" % (t_k1,std_k1,100*std_k1/k1))
        print("k2                  = %g +/- %g /s   (i.e. +/- %g %%)" % (t_k2,std_k2,100*std_k2/k2))
        print("k1prime             = %g +/- %g /M/s (i.e. +/- %g %%)" % (t_k1prime,std_k1prime,100*std_k1prime/k1prime))
        print("k2prime             = %g +/- %g /s   (i.e. +/- %g %%)" % (t_k2prime,std_k2prime,100*std_k2prime/k2prime))
        print("k_eff               = %g +/- %g /M/s (i.e. +/- %g %%) at %s" % (t_keff,std_keff,100*std_keff/keff,concentration_string(concentration))) 
        print("z_crit              = %s +/- %s (i.e. +/- %g %%)" % (concentration_string(t_zcrit),concentration_string(std_zcrit),100*std_zcrit/zcrit))
    print()

    return [N_forward, N_reverse, k1, k1prime, k2, k2prime, keff, zcrit, o]


if __name__ == '__main__':
    trials=1000 
    # Note that the "analytic" formulas for error bars agree well for N=1000, but are likely to disagree for N=100.

    print()
#     print "Simulating unstructured DNA strand, at 25 C."
#     data=first_step_simulation("ACTGGCGCGTATTATCTACTG", 1000, concentration=50e-9, T=25, material="DNA") 
#     print "Simulating DNA strand with 4-bp blunt-end hairpin, at 25 C."
#     data=first_step_simulation("ACTGGCGCGTATTATCTCAGT", 10000, concentration=50e-9, T=25, material="DNA") 
#     # We can do a lot of simulations in the above case, because the vast majority are very very fast.
#     # (In fact, most Boltzmann-sampled structures are hairpins, so no first base pair reactions are possible.)
#     # Note that k_collision_forward is very different from k_collision_reverse, because they sample markedly different initial structures.
#     print "Simulating DNA strand with 7-bp dual 5-nt tailed hairpin, at 25 C."
#     data=first_step_simulation("CTAACTGGCGCGTATTCGCGCCTTCAC", 1000, concentration=50e-9, T=25, material="DNA") 

    print("Simulating unstructured DNA strand hybridization from Gao et al 2006, at 20 C and 1 uM.")
    data0=first_step_simulation("GTTGTCAAGATGCTACCGTTCAGAG", trials, concentration=1e-6, T=20, material="DNA")
    print("Simulating tailed 3-bp hairpin DNA strand hybridization from Gao et al 2006, at 20 C and 1 uM.")
    data3=first_step_simulation("AGATCAGTGCGTCTGTACTAGCAGT", trials, concentration=1e-6, T=20, material="DNA")
    print("Simulating tailed 4-bp hairpin DNA strand hybridization from Gao et al 2006, at 20 C and 1 uM.")
    data4=first_step_simulation("AGATCAGTGCGTCTGTACTAGCACA", trials, concentration=1e-6, T=20, material="DNA")

    print("Gao et al report 12.0, 7.2, and 2.0 x 10^5 /M/s respectively for their 0-, 3-, and 4-bp hairpin DNA strands.")
    print("Multistrand observed %g, %g, and %g x 10^6 /M/s respectively, which is in qualitative agreement for the relative trends." % \
        (data0[6]/1e6, data3[6]/1e6, data4[6]/1e6))
    print()
    # Also note that for these concentrations, k_eff and k1 agree pretty darn well, for these molecules.

    # Reference for experimental rates:
    #    "Secondary structure effects on DNA hybridization kinetics: a solution versus surface comparison"
    #    Y Gao, LK Wolf, RM Georgiadis
    #    Nucleic acids research, vol. 34, pp. 3370-3377 (2006)
    #    http://nar.oxfordjournals.org/content/34/11/3370.short
    #
    # Reference for coarse-grained molecular dynamics study of the same:
    #    "DNA hairpins primarily promote duplex melting rather than inhibiting hybridization"
    #    John S. Schreck, Thomas E. Ouldridge, Flavio Romano, Petr Sulc, Liam Shaw, Ard A. Louis, Jonathan P. K. Doye
    #    arXiv:1408.4401 [cond-mat.soft]  (2014)
    #    http://arxiv.org/abs/1408.4401





## FD: My output:


# Simulating unstructured DNA strand hybridization from Gao et al 2006, at 20 C and 1 uM.
# Running first step mode simulations for GTTGTCAAGATGCTACCGTTCAGAG (with Boltzmann sampling)...
# simulation_mode = 48 
# simulation_count = 1000 
# o_interval = -1 
# o_time = -1 
# stop_options = 1 
# stop_count = 2 
# max_sim_time = 1 
# seed = 382255406 
# myComplexes = { { GTTGTCAAGATGCTACCGTTCAGAG, ..............(....)..... }{ CTCTGAACGGTAGCATCTTGACAAC, ..(((.....)))............ }} 
#  myStopComplexes = { } 
# temperature = 293.15 
# dangles = 1 
# logml = 0 
# gtenable = 0 
# kinetic_rate_method = 1 
# joinConcentration = 1 
# biScale = 1.4e+06 
# uniScale = 7.3e+08 
#  substrate_type = 0 
# 
# Inferred rate constants with analytical error bars:
# N_forward = 707
# N_reverse = 293
# k_collision_forward = 6.1101e+07 +/- 1.20677e+06 /M/s (i.e. +/- 1.97504 %)
# k_collision_reverse = 5.52785e+07 +/- 1.68592e+06 /M/s (i.e. +/- 3.04986 %)
# k1                  = 4.31984e+07 +/- 1.22527e+06 /M/s (i.e. +/- 2.83638 %)
# k2                  = 1.224e+06 +/- 71196.4 /s   (i.e. +/- 5.81672 %)
# k1prime             = 1.61966e+07 +/- 936485 /M/s (i.e. +/- 5.78199 %)
# k2prime             = 2.21389e+06 +/- 218127 /s   (i.e. +/- 9.85265 %)
# k_eff               = 4.31966e+07 /M/s at 1.0 uM (still needs error bars)
# z_crit              = 23.4693015469 mM (still needs error bars)
# 
# Re-sampled rate constants with bootstrapped error bars:
# N_forward = 718 +/- 14.5006
# N_reverse = 282 +/- 14.5006
# k_collision_forward = 6.19822e+07 +/- 1.19648e+06 /M/s (i.e. +/- 1.95821 %)
# k_collision_reverse = 5.71617e+07 +/- 1.6516e+06 /M/s (i.e. +/- 2.98777 %)
# k1                  = 4.45032e+07 +/- 1.23164e+06 /M/s (i.e. +/- 2.85112 %)
# k2                  = 1.24367e+06 +/- 72116.4 /s   (i.e. +/- 5.89189 %)
# k1prime             = 1.61196e+07 +/- 921377 /M/s (i.e. +/- 5.68871 %)
# k2prime             = 2.47933e+06 +/- 220825 /s   (i.e. +/- 9.97454 %)
# k_eff               = 4.45013e+07 +/- 1.23154e+06 /M/s (i.e. +/- 2.85102 %) at 1.0 uM
# z_crit              = 23.6489126432 mM +/- 1.34052276489 mM (i.e. +/- 5.71181 %)
# 
# Simulating tailed 3-bp hairpin DNA strand hybridization from Gao et al 2006, at 20 C and 1 uM.
# Running first step mode simulations for AGATCAGTGCGTCTGTACTAGCAGT (with Boltzmann sampling)...
# simulation_mode = 48 
# simulation_count = 1000 
# o_interval = -1 
# o_time = -1 
# stop_options = 1 
# stop_count = 2 
# max_sim_time = 1 
# seed = 1507406588 
# myComplexes = { { AGATCAGTGCGTCTGTACTAGCAGT, ((..((.......))..))...... }{ ACTGCTAGTACAGACGCACTGATCT, ((.....)).(((.....))).... }} 
#  myStopComplexes = { } 
# temperature = 293.15 
# dangles = 1 
# logml = 0 
# gtenable = 0 
# kinetic_rate_method = 1 
# joinConcentration = 1 
# biScale = 1.4e+06 
# uniScale = 7.3e+08 
#  substrate_type = 0 
# 
# Inferred rate constants with analytical error bars:
# N_forward = 407
# N_reverse = 593
# k_collision_forward = 2.74909e+07 +/- 826213 /M/s (i.e. +/- 3.00541 %)
# k_collision_reverse = 2.3545e+07 +/- 503948 /M/s (i.e. +/- 2.14036 %)
# k1                  = 1.11888e+07 +/- 543579 /M/s (i.e. +/- 4.85824 %)
# k2                  = 111262 +/- 10035.5 /s   (i.e. +/- 9.01972 %)
# k1prime             = 1.39622e+07 +/- 472338 /M/s (i.e. +/- 3.38298 %)
# k2prime             = 587752 +/- 120604 /s   (i.e. +/- 20.5196 %)
# k_eff               = 1.11874e+07 /M/s at 1.0 uM (still needs error bars)
# z_crit              = 8.04389016313 mM (still needs error bars)
# 
# Re-sampled rate constants with bootstrapped error bars:
# N_forward = 434 +/- 16.3619
# N_reverse = 566 +/- 16.3619
# k_collision_forward = 2.67452e+07 +/- 822781 /M/s (i.e. +/- 2.99292 %)
# k_collision_reverse = 2.35576e+07 +/- 498290 /M/s (i.e. +/- 2.11633 %)
# k1                  = 1.16074e+07 +/- 560847 /M/s (i.e. +/- 5.01257 %)
# k2                  = 108529 +/- 10331.3 /s   (i.e. +/- 9.28554 %)
# k1prime             = 1.33336e+07 +/- 495807 /M/s (i.e. +/- 3.55107 %)
# k2prime             = 598863 +/- 129415 /s   (i.e. +/- 22.0186 %)
# k_eff               = 1.16059e+07 +/- 560720 /M/s (i.e. +/- 5.01206 %) at 1.0 uM
# z_crit              = 7.73892417558 mM +/- 754.985543712 uM (i.e. +/- 9.38583 %)
# 
# Simulating tailed 4-bp hairpin DNA strand hybridization from Gao et al 2006, at 20 C and 1 uM.
# Running first step mode simulations for AGATCAGTGCGTCTGTACTAGCACA (with Boltzmann sampling)...
# simulation_mode = 48 
# simulation_count = 1000 
# o_interval = -1 
# o_time = -1 
# stop_options = 1 
# stop_count = 2 
# max_sim_time = 1 
# seed = 1973482372 
# myComplexes = { { AGATCAGTGCGTCTGTACTAGCACA, ......((((((....))..)))). }{ TGTGCTAGTACAGACGCACTGATCT, .((((..........))))...... }} 
#  myStopComplexes = { } 
# temperature = 293.15 
# dangles = 1 
# logml = 0 
# gtenable = 0 
# kinetic_rate_method = 1 
# joinConcentration = 1 
# biScale = 1.4e+06 
# uniScale = 7.3e+08 
#  substrate_type = 0 
# 
# Inferred rate constants with analytical error bars:
# N_forward = 381
# N_reverse = 619
# k_collision_forward = 2.8327e+07 +/- 364877 /M/s (i.e. +/- 1.28809 %)
# k_collision_reverse = 2.88346e+07 +/- 298932 /M/s (i.e. +/- 1.03671 %)
# k1                  = 1.07926e+07 +/- 456693 /M/s (i.e. +/- 4.23153 %)
# k2                  = 49917.1 +/- 2644.12 /s   (i.e. +/- 5.29703 %)
# k1prime             = 1.78486e+07 +/- 479920 /M/s (i.e. +/- 2.68884 %)
# k2prime             = 1.91464e+06 +/- 423070 /s   (i.e. +/- 22.0965 %)
# k_eff               = 1.07902e+07 /M/s at 1.0 uM (still needs error bars)
# z_crit              = 4.43394486229 mM (still needs error bars)
# 
# Re-sampled rate constants with bootstrapped error bars:
# N_forward = 392 +/- 15.0265
# N_reverse = 608 +/- 15.0265
# k_collision_forward = 2.80464e+07 +/- 362484 /M/s (i.e. +/- 1.27964 %)
# k_collision_reverse = 2.83984e+07 +/- 301484 /M/s (i.e. +/- 1.04556 %)
# k1                  = 1.09942e+07 +/- 447421 /M/s (i.e. +/- 4.14563 %)
# k2                  = 48988.7 +/- 2689.9 /s   (i.e. +/- 5.38874 %)
# k1prime             = 1.72662e+07 +/- 475332 /M/s (i.e. +/- 2.66313 %)
# k2prime             = 1.72001e+06 +/- 455941 /s   (i.e. +/- 23.8133 %)
# k_eff               = 1.09916e+07 +/- 447231 /M/s (i.e. +/- 4.1448 %) at 1.0 uM
# z_crit              = 4.26509130373 mM +/- 284.957346709 uM (i.e. +/- 6.42672 %)
# 
# Gao et al report 12.0, 7.2, and 2.0 x 10^5 /M/s respectively for their 0-, 3-, and 4-bp hairpin DNA strands.
# Multistrand observed 43.1966, 11.1874, and 10.7902 x 10^6 /M/s respectively, which is in qualitative agreement for the relative trends.
