####################################################################
#                                                                  #
#  Copyright (c) 2010-2015 California Institute of Technology.     #
#  Distributed under the MIT License.                              #
#  (See accompanying file LICENSE or copy at                       #
#  http://opensource.org/licenses/MIT)                             #
#                                                                  #
####################################################################

# hybridization_comparison.py
# 
# Again using hybridization of two short strands as an example, here we show that 
# multiple ways of extracting rate constants from Multistrand simulations give the 
# same result, at least if the system is simple enough that the model 
#   A + B -->{k_eff} C 
#       C -->{K_rev} A + B
# is appropriate.  We emphasize that only *short* strands can be studied this way,
# because simulations of the dissociation of long duplexes are prohibitively slow.
#
# Normally, k_eff should be the concentration-independent second-order bimolecular
# rate constant for association.  However, for sufficiently high concentrations, the
# unimolecular steps or unproductive interactions will be rate-limiting, and k_eff
# (inferred from data to fit the above model) will decrease.  Thus, for the general 
# case we consider k_eff to be a concentration-dependent approximation, in contrast
# to the more fundamental rate constants k1, k2, k1prime, k2prime derived in
# hybridization_first_step_mode.py.  
#
# We consider three ways to obtain k_eff.  
# 1.  First step mode.  Estimate more fundamental rate constants from unimolecular 
#     simulations of collisions, from which k_eff can be computed as a function of 
#     concentration.  This is the fastest.  See hybridization_first_step_mode.py 
#     for more explanation.
# 2.  First passage time simulations starting with (Boltzmann sampled) individual
#     strands within a concentration-dependent volume, wait for them to collide
#     as many times as necessary until a duplex is formed, and compute the rate
#     constant k_eff (for that concentration) from the average time.
# 3.  Transition mode simulations, where two strands are placed in a concentration-
#     dependent volume and simulated for a good long time as they hybridize and
#     fall apart again and again. Rate constants for both k_eff and k_rev are
#     estimated from the average times to transition from single-stranded to 
#     double-stranded macrostates (for k_eff) and back (for k_rev).
#
# There are also three ways to obtain k_rev.
# 1.  Trust the first step mode k_eff value, trust that a bimolecular model
#     is sufficient, so k_rev/k_eff = exp(dG/RT), and compute dG from the
#     NUPACK-predicted partition functions of the single strands and the duplex.
#     This is the fastest.
# 2.  First passage time simulations starting with the duplex state, stopping
#     when dissociation occurs.
# 3.  The same transition mode simulation described above.
#
# Although the calculations for these three methods are quite different, for a 
# suffient number of simulation trials they agree quite well.  Method 1 is faster
# than method 2 is faster than method 3.  (A proper comparison would also evaluate
# error bars, so that the amount of compute time required for a given level of accuracy
# is fairly compared.  We don't do that here.)
#
# A significant issue when using methods 2 and 3 is that, a priori, you won't know what 
# the critical concentration is, below which k_eff is constant.  (This value is usually 
# what you want.)  On the other hand, it is desirable to do first passage time and 
# transition mode simulations at a high concentration, so as not to waste time doing 
# simulations of single-stranded molecules just waiting for a collision to happen.  Practically,
# one must do a simulation at a high concentration, then decrease the concentration by a factor of
# two and perform another simulation, verifying that the estimated k_eff values don't change
# significantly.  First step mode simulations are not only faster, but they also obviate this
# concern by separately extracting the bimolecular (k1) and unimolecular (k2) rate constants.
#
# Try it like this, e.g.:
#   python -i hybridization_comparison.py
# The default 5-mer hybridization trials take about 10 minutes to complete.  
# Other options can be uncommented at the bottom of the file.

import numpy as np
import matplotlib
import matplotlib.pylab as plt
import time

if False:  # only needed if you're having trouble with your Multistrand installation
    import multistrand_setup

try:
    from multistrand.objects import *
    from multistrand.options import Options
    from multistrand.system import SimSystem, initialize_energy_model

except ImportError:
    print("Could not import Multistrand.")
    raise

######### This is first passage time

# for StopCondition and Macrostate definitions:
Exact_Macrostate = 0   # match a secondary structure exactly (i.e. any system state that has a complex with this exact structure)
Bound_Macrostate = 1   # match any system state in which the given strand is bound to another strand
Dissoc_Macrostate = 2  # match any system state in which there exists a complex with exactly the given strands, in that order
Loose_Macrostate = 3   # match a secondary structure with "don't care"s, allowing a certain number of disagreements
Count_Macrostate = 4   # match a secondary structure, allowing a certain number of disagreements
# see Schaeffer's PhD thesis, chapter 7.2, for more information

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


def first_step_simulation(strand_seq, trials, T=25, material="DNA"):

   print "Running first step mode simulations for %s (with Boltzmann sampling)..." % (strand_seq)

   # Using domain representation makes it easier to write secondary structures.
   onedomain = Domain(name="itall",sequence=strand_seq)
   top = Strand(name="top",domains=[onedomain])
   bot = top.C

   # Note that the structure is specified to be single stranded, but this will be over-ridden when Boltzmann sampling is turned on.
   start_complex_top = Complex(strands=[top],structure=".")
   start_complex_bot = Complex(strands=[bot],structure=".")
   start_complex_top.boltzmann_count = trials
   start_complex_bot.boltzmann_count = trials
   start_complex_top.boltzmann_sample = True
   start_complex_bot.boltzmann_sample = True
   # Turns Boltzmann sampling on for this complex and also does sampling more efficiently by sampling 'trials' states.

   # Stop when the exact full duplex is achieved. (No breathing!)
   success_complex = Complex(strands=[top, bot],structure="(+)")
   success_stop_condition = StopCondition("SUCCESS",[(success_complex,Exact_Macrostate,0)])

   # Declare the simulation unproductive if the strands become single-stranded again.
   failed_complex = Complex(strands = [top], structure=".")
   failed_stop_condition = StopCondition("FAILURE",[(failed_complex,Dissoc_Macrostate,0)])

   o = Options(simulation_mode="First Step",parameter_type="Nupack", substrate_type=material,
               rate_method = "Metropolis", num_simulations = trials, simulation_time=1.0,
               dangles = "Some", temperature = T, rate_scaling = "Calibrated", verbosity = 0)

   o.start_state = [start_complex_top, start_complex_bot]
   o.stop_conditions = [success_stop_condition,failed_stop_condition]

   # Now go ahead and run the simulations.
   initialize_energy_model(o)  # concentration changes, so we must make sure energies are right
   s = SimSystem(o)
   s.start()
   dataset = o.interface.results

   # Now determine the reaction model parameters from the simulation results.  (Simplified from hybridization_first_step_mode.py.)
   collision_rates = np.array( [i.collision_rate for i in dataset] )
   was_success = np.array([1 if i.tag=="SUCCESS" else 0 for i in dataset])
   was_failure = np.array([0 if i.tag=="SUCCESS" else 1 for i in dataset])
   forward_times = np.array( [i.time for i in dataset if i.tag == "SUCCESS"] )
   reverse_times = np.array( [i.time for i in dataset if i.tag == "FAILURE" or i.tag == None] )

   # Calculate first-order rate constants for the duration of the reactions (both productive and unproductive).
   k2 = 1.0/np.mean(forward_times)
   k2prime = 1.0/np.mean(reverse_times)

   # Calculate second-order rate constants for productive and unproductive reactions.
   k1 = np.mean( collision_rates * was_success )
   k1prime = np.mean( collision_rates * was_failure )

   return k1, k2, k1prime, k2prime


def first_passage_dissociation(strand_seq, trials, T=25, material="DNA"):

   print "Running first passage time simulations for dissociation of %s..." % (strand_seq)

   # Using domain representation makes it easier to write secondary structures.
   onedomain = Domain(name="itall",sequence=strand_seq)
   top = Strand(name="top",domains=[onedomain])
   bot = top.C
   single_strand_top = Complex(strands=[top],structure=".")
   single_strand_bot = Complex(strands=[bot],structure=".")
   duplex_complex = Complex(strands=[top, bot],structure="(+)")

   # Declare the simulation complete if the strands become single-stranded again.
   success_stop_condition = StopCondition("SUCCESS",[(single_strand_top,Dissoc_Macrostate,0)])

   o = Options(simulation_mode="First Passage Time",parameter_type="Nupack", substrate_type=material,
               rate_method = "Metropolis", num_simulations = trials, simulation_time=10.0, 
               join_concentration=1e-6, # 1 uM concentration, but doesn't matter for dissociation
               dangles = "Some", temperature = T, rate_scaling = "Calibrated", verbosity = 0)
   o.start_state = [duplex_complex]
   o.stop_conditions = [success_stop_condition]

   # Now go ahead and run the simulations.
   initialize_energy_model(o)  # concentration changes, so we must make sure energies are right
   s = SimSystem(o)
   s.start()
   dataset = o.interface.results

   times = np.array([i.time for i in dataset])
   timeouts = [i for i in dataset if not i.tag == 'SUCCESS']
   if len(timeouts)>0 :
        print "some dissociation trajectories did not finish..."
        for i in timeouts :
            assert (i.type_name=='Time')
            assert (i.tag == None )
            assert (i.time >= 10.0)
   
   krev = 1.0/np.mean( times )

   return krev

def first_passage_association(strand_seq, trials, concentration, T=25, material="DNA"):

   print "Running first passage time simulations for association of %s at %s..." % (strand_seq, concentration_string(concentration))

   # Using domain representation makes it easier to write secondary structures.
   onedomain = Domain(name="itall",sequence=strand_seq)
   top = Strand(name="top",domains=[onedomain])
   bot = top.C
   duplex_complex = Complex(strands=[top, bot],structure="(+)")
   single_strand_top = Complex(strands=[top],structure=".")
   single_strand_bot = Complex(strands=[bot],structure=".")
   # Start with Boltzmann-sampled single-strands... it only seems fair.
   single_strand_top.boltzmann_count = trials
   single_strand_bot.boltzmann_count = trials
   single_strand_top.boltzmann_sample = True
   single_strand_bot.boltzmann_sample = True

   # Declare the simulation complete if the strands become a perfect duplex.
   success_stop_condition = StopCondition("SUCCESS",[(duplex_complex,Exact_Macrostate,0)])

   o = Options(simulation_mode="First Passage Time",parameter_type="Nupack", substrate_type=material,
               rate_method = "Metropolis", num_simulations = trials, simulation_time=10.0, 
               join_concentration=concentration,
               dangles = "Some", temperature = T, rate_scaling = "Calibrated", verbosity = 0)
   o.start_state = [single_strand_top, single_strand_bot]
   o.stop_conditions = [success_stop_condition]

   # Now go ahead and run the simulations.
   initialize_energy_model(o)  # concentration changes, so we must make sure energies are right
   s = SimSystem(o)
   s.start()
   dataset = o.interface.results

   times = np.array([i.time for i in dataset])
   timeouts = [i for i in dataset if not i.tag == 'SUCCESS']
   if len(timeouts)>0 :
        print "some association trajectories did not finish..."
        for i in timeouts :
            assert (i.type_name=='Time')
            assert (i.tag == None )
            assert (i.time >= 10.0)
   
   print "average completion time = %g seconds at %s" % (np.mean(times),concentration_string(concentration))

   keff = 1.0/np.mean( times )/concentration

   return keff

######## stuff for transition mode

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

def print_transition_dict( transition_dict, options = None ):
    k = transition_dict.keys()
    k.sort() 

    for i in k:
        transition_times = np.array( transition_dict[i] )
        print("{0}: {2:.2e} s ({1} events)".format(i,len(transition_dict[i]),np.mean(transition_times)))
    
    # also print the true names of the macrostates, if an Options object is provided
    charindex = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0'
    if options:
        for i,idx in zip(options.stop_conditions,range(len(options.stop_conditions))):
            print("{0}: {1}".format( i.tag, charindex[idx]))


def transition_mode_simulation(strand_seq, duration, concentration, T=25, material="DNA"):

   print "Running transition mode simulations of %s at %s..." % (strand_seq, concentration_string(concentration))

   # Using domain representation makes it easier to write secondary structures.
   onedomain = Domain(name="itall",sequence=strand_seq)
   top = Strand(name="top",domains=[onedomain])
   bot = top.C
   duplex_complex = Complex(strands=[top, bot],structure="(+)")
   single_strand_top = Complex(strands=[top],structure=".")
   single_strand_bot = Complex(strands=[bot],structure=".")

   # Declare macrostates 
   single_stranded_macrostate = Macrostate("SINGLE",[(single_strand_top,Dissoc_Macrostate,0)])
   duplex_macrostate = Macrostate("DUPLEX",[(duplex_complex,Loose_Macrostate,4)])

   o = Options(simulation_mode="Transition",parameter_type="Nupack", substrate_type=material,
               rate_method = "Metropolis", num_simulations = 1, simulation_time=float(duration),   # time must be passed as float, not int
               join_concentration=concentration,
               dangles = "Some", temperature = T, rate_scaling = "Calibrated", verbosity = 0)
   o.start_state = [single_strand_top, single_strand_bot]
   o.stop_conditions = [single_stranded_macrostate, duplex_macrostate] # not actually stopping, just tracking

   # Now go ahead and run the simulations until time-out.
   initialize_energy_model(o)  # concentration changes, so we must make sure energies are right
   s = SimSystem(o)
   s.start()
   
   # Now make sense of the results.
   transition_dict = parse_transition_lists(o.interface.transition_lists)
   print_transition_dict( transition_dict, o )

   # A is SINGLE, B is DUPLEX
   N_AtoA = float(len( transition_dict['A -> A'] )) if 'A -> A' in transition_dict else 0
   dT_AtoA = np.mean( transition_dict['A -> A'] ) if N_AtoA > 0 else 1 # will be mult by zero in that case
   N_AtoB = float(len( transition_dict['A -> B'] )) if 'A -> B' in transition_dict else 0
   dT_AtoB = np.mean( transition_dict['A -> B'] ) if N_AtoB > 0 else 1
   N_BtoB = float(len( transition_dict['B -> B'] )) if 'B -> B' in transition_dict else 0
   dT_BtoB = np.mean( transition_dict['B -> B'] ) if N_BtoB > 0 else 1
   N_BtoA = float(len( transition_dict['B -> A'] )) if 'B -> A' in transition_dict else 0
   dT_BtoA = np.mean( transition_dict['B -> A'] ) if N_BtoA > 0 else 1

   keff = 1.0/(dT_AtoB + (N_AtoA/N_AtoB)*dT_AtoA)/concentration if N_AtoB > 0 else None
   krev = 1.0/(dT_BtoA + (N_BtoB/N_BtoA)*dT_BtoB) if N_BtoA > 0 else None
   
   return keff, krev

def compare_hybridization(seq, concentrations, T=25, material="DNA"):
   # track time for each kind of simulation, using time.time(), which has units of seconds
   time1=time.time()

   # do one "first step mode" run, get k1, k2, etc, from which z_crit and k_eff(z) can be computed
   k1, k2, k1prime, k2prime = first_step_simulation(seq, 10000, T=T, material=material) 
   time2=time.time()
   time1step_for = time2-time1
   print "k1 = %g /M/s, k2 = %g /s, k1prime = %g /M/s, k2prime = %g /s  (%g seconds)" % (k1, k2, k1prime, k2prime, time1step_for)
   zcrit = k2*k2prime/(k1*k2prime + k1prime*k2) # this is the critical concentration at which k_eff = k1/2
   print "zcrit = %s" % (concentration_string(zcrit))
   for z in concentrations:
       keff_1s = 1/(1/k1 + z/k2 + (k1prime/k1)*(z/k2prime))  # first-step mode predictions
       print "keff = %g /M/s at %s" % (keff_1s, concentration_string(z))
   print

   # call NUPACK for pfunc dG of the reaction, calculate krev based on keff
   print "Calculating dissociate rate constant based on NUPACK partition function energies and first step mode k_eff..."
   import nupack
   dG_top = nupack.pfunc([seq], T=T)
   dG_bot = nupack.pfunc([ Strand(sequence=seq).C.sequence ], T=T)
   dG_duplex = nupack.pfunc([ seq, Strand(sequence=seq).C.sequence ], T=T)
   RT = 1.987e-3 * (273.15+T)
   time3=time.time()
   time_nupack = time3-time2
   krev_nupack = keff_1s * np.exp( (dG_duplex - dG_top - dG_bot)/RT )
   print "krev = %g /s (%g seconds)" % (krev_nupack, time_nupack)
   print

   # do one "first passage time" run for dissociation, and get k_rev
   krev_1p = first_passage_dissociation(seq, 500, T=T, material=material) 
   time4=time.time()
   time1passage_rev = time4-time3
   print "krev = %g /s (%g seconds)" % (krev_1p, time1passage_rev)
   print

   # for each concentration z, do one "first passage time" run for association, and get k_eff(z)
   keffs_1p=[]
   for concentration in concentrations:
       keff = first_passage_association(seq, 10000, concentration=concentration, T=T, material=material)
       keffs_1p.append((keff,concentration))
   time5=time.time()
   time1passage_for = time5-time4
   for (keff,z) in keffs_1p:
       print "keff = %g /M/s at %s" % (keff, concentration_string(z))
   print "(took %g seconds total)" % (time1passage_for)
   print

   # for each concentration z, do one long "transition mode" run, and compute k_eff(z) and k_rev for each run.
   keffs_tm = []
   krevs_tm = []
   for concentration in concentrations:
       keff, krev = transition_mode_simulation(seq, 0.5, concentration, T=T, material=material)
       keffs_tm.append((keff,concentration))
       if not keff is None:
           print "keff = %g /M/s at %s" % (keff, concentration_string(concentration))
       krevs_tm.append((krev,concentration))
       if not krev is None:
           print "krev = %g /s at %s" % (krev, concentration_string(concentration))
   time6=time.time()
   time_transitions = time6-time5
   print "(took %g seconds total)" % (time_transitions)
   print

   # One could plot k_eff vs z, and k_rev vs z, comparing the methods.

if __name__ == '__main__':

    # dissociation and transision simulations are slow -- let's consider short strands here!
    compare_hybridization(seq='TCGAT', concentrations=[1e-2,1e-3,1e-4,1e-5,1e-6,1e-7]) # takes about 10 minutes

    ### The following is super-slow for some modes, but it will run.  Expect 1 day or more.  
    # compare_hybridization(seq='CGTTTCG', concentrations=[1e-1,1e-2,1e-3,1e-4,1e-5]) 

    ### These three are super super super slow, and did not finish all modes for me:
    # compare_hybridization(seq='TCGATTTTGTA', concentrations=[1e-3,1e-4,1e-5])  
    # compare_hybridization(seq='TCGATTTTTCGA', concentrations=[1e-3,1e-4,1e-5])
    # compare_hybridization(seq='ACTGGCGCGTATTATCTACTG', concentrations=[1e-3,1e-4,1e-5])
