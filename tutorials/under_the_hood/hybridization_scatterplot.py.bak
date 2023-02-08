# hybridization_scatterplot.py
#
# Here we use the function calls defined in hybridization_first_step_mode.py
# to examine the influence of secondary structure on hybridization rates,
# according to the Multistrand model.  The inspiration is a pair of papers,
# Gao et al 2006 and Schreck et al 2014 (full citations at the bottom), the
# latter of which argues that secondary structure accelerates dissociation
# rather than inhibiting association -- at least in the case of the experimental
# data provided by Gao et al -- and performs coarse-grained molecular dynamics
# simulations as evidence. Our goal is to refine that insight by looking at 
# a larger number of sequences and understanding the role of "toeholds" and
# hairpin stem length in determining whether/when secondary structure accelerates 
# dissociation or inhibits association.
# 
# Usage:
#     python -i hybridization_scatterplot generate random <len> <num seqs> <data file name>
#     python -i hybridization_scatterplot generate iso-random <len> <num seqs> <data file name>
#     python -i hybridization_scatterplot generate structured <len> <max toe> <max stem> <num seqs> <data file name>
#     python -i hybridization_scatterplot generate iso-structured <len> <max toe> <max stem> <num seqs> <data file name>
#     python -i hybridization_scatterplot generate fixed-stem <len> <max toe> <stem size> <num seqs> <data file name>
#     python -i hybridization_scatterplot generate iso-fixed-stem <len> <max toe> <stem size> <num seqs> <data file name>
#     python -i hybridization_scatterplot plot <data file names>
#
# Usage examples:
#     python -i hybridization_scatterplot.py generate random 15 1000 r15-1000
#     python -i hybridization_scatterplot.py generate iso-random 25 1000 ir25-1000
#     python -i hybridization_scatterplot.py generate structured 15 5 5 100 s15-100
#     python -i hybridization_scatterplot.py generate iso-structured 15 5 5 300 is15-300
#     python -i hybridization_scatterplot.py generate iso-structured 25 8 8 100 is25-100
#     python -i hybridization_scatterplot.py generate fixed-stem 15 5 4 100 f15-100
#     python -i hybridization_scatterplot.py generate iso-fixed-stem 25 8 5 100 if25-100
#     python -i hybridization_scatterplot.py plot r15-1000 s15-100 f15-100
#     open scatterdata/r15-1000+s15-100+f15-100.pdf
#

# Reference for experimental rates:
#    "Secondary structure effects on DNA hybridization kinetics: a solution versus surface comparison"
#    Y Gao, LK Wolf, RM Georgiadis
#    Nucleic acids research, vol. 34, pp. 3370-3377 (2006)
#    http://nar.oxfordjournals.org/content/34/11/3370.short
#
# Reference for coarse-grained molecular dynamics study of the same:
#    "DNA hairpins primarily promote duplex melting rather than inhibiting hybridization"
#    John S. Schreck, Thomas E. Ouldridge, Flavio Romano, Petr Sulc, Liam Shaw, Ard A. Louis, Jonathan P. K. Doye
#    Nucleic Acids Research, 2015, Vol. 43, No. 13 6181-6190
#    arXiv:1408.4401 [cond-mat.soft]  (2014)
#    http://arxiv.org/abs/1408.4401




# Import things you need
import random, string, pickle, sys
import numpy as np
import nupack

if False:  # only needed if you're having trouble with your Multistrand installation
    import multistrand_setup

try:
    from multistrand.objects import *
    from multistrand.options import Options, Literals
    from multistrand.system import SimSystem, initialize_energy_model

except ImportError:
    print("Could not import Multistrand.")
    raise

#########

# for StopCondition macrostate definitions:
Exact_Macrostate = 0
Bound_Macrostate = 1
Dissoc_Macrostate = 2
Loose_Macrostate = 3
Count_Macrostate = 4

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
    # for example, if they are both blunt-ended hairpins.  In this case, i.collision_rate will be 0, i.tag "noinitial"
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

    # problems arise below if either N_forward == 0 or N_reverse == 0... should clean up code to handle it more smoothly

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
        print "N_forward =", N_forward
        print "N_reverse =", N_reverse
        # print "k_collision = %g +/- %g /M/s (i.e. +/- %g %%)" % (kcollision,std_kcollision,100*std_kcollision/kcollision)
        print "k_collision_forward = %g +/- %g /M/s (i.e. +/- %g %%)" % (forward_kcoll,std_forward_kcoll,100*std_forward_kcoll/forward_kcoll)
        print "k_collision_reverse = %g +/- %g /M/s (i.e. +/- %g %%)" % (reverse_kcoll,std_reverse_kcoll,100*std_reverse_kcoll/reverse_kcoll)
        print "k1                  = %g +/- %g /M/s (i.e. +/- %g %%)" % (k1,std_k1,100*std_k1/k1)
        print "k2                  = %g +/- %g /s   (i.e. +/- %g %%)" % (k2,std_k2,100*std_k2/k2)
        print "k1prime             = %g +/- %g /M/s (i.e. +/- %g %%)" % (k1prime,std_k1prime,100*std_k1prime/k1prime)
        print "k2prime             = %g +/- %g /s   (i.e. +/- %g %%)" % (k2prime,std_k2prime,100*std_k2prime/k2prime)
        print "k_eff               = %g /M/s at %s (still needs error bars)" % (keff,concentration_string(concentration)) 
        print "z_crit              = %s (still needs error bars)" % (concentration_string(zcrit)) 

    return N_forward, N_reverse, kcollision, forward_kcoll, reverse_kcoll, k1, k2, k1prime, k2prime, keff, zcrit

def resample_with_replacement(mylist,num_samples):
    return [random.choice(mylist) for i in range(num_samples)]

def first_step_simulation(strand_seq, num_traj, T=25, rate_method_k_or_m="Metropolis", concentration=50e-9, material="DNA"):

    # Run the simulations

    print "Running %d first step mode simulations for %s (with Boltzmann sampling)..." % (num_traj,strand_seq)
    o = create_setup(strand_seq, num_traj, T, rate_method_k_or_m, material)
    initialize_energy_model(o)  # Prior simulations could have been for different temperature, material, etc.
                                # But Multistrand "optimizes" by sharing the energy model parameters from sim to sim.
                                # So if in the same python session you have changed parameters, you must re-initialize.
    s = SimSystem(o)
    s.start()
    return o
    

def WC(seq):
    """Computes the Watson-Crick complement for a DNA sequence."""
    return seq.translate(string.maketrans('ACTG','TGAC'))[::-1]

def randomseq(length,bases='ACTG'):
    """Chooses a random DNA sequence."""
    return ''.join(random.choice(bases) for i in range(length))

def stemsize(seq, T=25, material='dna'):
    result = nupack.mfe([seq], T=T, material=material)
    struct = result[0][0]
    n = len(struct)
    try: 
        stemstart5 = struct.index('(')
        stemlen5   = struct[stemstart5:].index('.')
        revstruct  = struct[::-1]
        stemstart3 = revstruct.index(')')
        stemlen3   = revstruct[stemstart3:].index('.')
        stemlen    = min(stemlen5,stemlen3)
    except ValueError:
        stemlen = 0
    return stemlen

def toeholds(seq, T=25, material="dna"):
    result = nupack.mfe([seq], T=T, material=material)
    struct = result[0][0]
    n = len(struct)
    try: 
        toe5len = struct.index('(')
        toe3len = struct[::-1].index(')')
    except ValueError:
        toe5len = n
        toe3len = n
    # print "    %s = %s , toes %d and %d" % (seq, struct, toe5len, toe3len)
    return (toe5len, toe3len)
    
def binding_dG(seq, T=25, material='dna'):
    dG_top = nupack.pfunc([seq], T=T, material=material)
    dG_bot = nupack.pfunc([ WC(seq) ], T=T, material=material)
    dG_duplex = nupack.pfunc([ seq, WC(seq) ], T=T, material=material)
    return (dG_duplex - dG_top - dG_bot)

def duplex_dG(seq, T=25, material='dna'):
    return nupack.pfunc([ seq, WC(seq) ], T=T, material=material)

def strand_dG(seq, T=25, material='dna'):
    return nupack.pfunc([ seq ], T=T, material=material)

def reverse_rate(seq, kf, T=25, material='dna'):
    dG   = binding_dG(seq, T=25, material='dna')
    RT   = 0.001987 * (273.15+25)  # kcal/mol at 25 C
    kr   = kf*np.exp( dG/RT )
    return kr

# this just looks at struct of sense strand.  it would be better to look at both, even if they have different structures. TO DO.
def toebinding(seq, T=25, material='dna'):
    n = len(seq)
    (toe5len, toe3len) = toeholds(seq, T=T, material=material)   
    if toe5len < n:
        toe5seq = seq[:toe5len] if toe5len > 0 else ''
        toe3seq = seq[-toe3len:] if toe3len > 0 else ''
        base = seq[toe5len]
        fake_seq1 = toe5seq + base + 'GGGGTTTTCCCC' + WC(base) + toe3seq
        fake_struct1 = '.'*toe5len + '(((((....)))))' + '.'*toe3len
        fake_seq2 = WC(fake_seq1)
        fake_struct2 = '.'*toe3len + '(((((....)))))' + '.'*toe5len
        fake_struct_toe5 = '('*toe5len + '(((((....)))))' + '.'*toe3len + '+' + '.'*toe3len + '(((((....)))))' + ')'*toe5len
        fake_struct_toe3 = '.'*toe5len + '(((((....)))))' + '('*toe3len + '+' + ')'*toe3len + '(((((....)))))' + '.'*toe5len
        fake_struct_toes = '('*toe5len + '(((((....)))))' + '('*toe3len + '+' + ')'*toe3len + '(((((....)))))' + ')'*toe5len
        # print "Sequence %s:  %s %s and %s %s" % (seq,fake_seq1,fake_struct1,fake_seq2,fake_struct2)
        dG1 = nupack.energy([fake_seq1],fake_struct1,T=T,material=material)
        dG2 = nupack.energy([fake_seq2],fake_struct2,T=T,material=material)
        if toe5len > 0:
            dG_toe5 = nupack.energy([fake_seq1,fake_seq2],fake_struct_toe5,T=T,material=material) - dG1 - dG2
        else:
            dG_toe5 = 0
        if toe3len > 0:
            dG_toe3 = nupack.energy([fake_seq1,fake_seq2],fake_struct_toe3,T=T,material=material) - dG1 - dG2
        else:
            dG_toe3 = 0
        if toe5len > 0 and toe3len > 0:
            dG_toes = nupack.energy([fake_seq1,fake_seq2],fake_struct_toes,T=T,material=material) - dG1 - dG2
        else:
            dG_toes = 0
        return max(-dG_toe5,-dG_toe3,-dG_toes)
    else:
        # print "Sequence %s: no structure!" % seq
        return -duplex_dG(seq, T=T, material=material)

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print """Usage:
              python -i hybridization_scatterplot generate random <len> <num seqs> <data file name>
              python -i hybridization_scatterplot generate iso-random <len> <num seqs> <data file name>
              python -i hybridization_scatterplot generate structured <len> <max toe> <max stem> <num seqs> <data file name>
              python -i hybridization_scatterplot generate iso-structured <len> <max toe> <max stem> <num seqs> <data file name>
              python -i hybridization_scatterplot generate fixed-stem <len> <max toe> <stem size> <num seqs> <data file name>
              python -i hybridization_scatterplot generate iso-fixed-stem <len> <max toe> <stem size> <num seqs> <data file name>
              python -i hybridization_scatterplot plot <data file names>
              """
        sys.exit()

    if sys.argv[1] == 'generate':
        if sys.argv[2] in ['random','iso-random'] and len(sys.argv) == 6:
            L = int(sys.argv[3])
            N = int(sys.argv[4])
            filename = sys.argv[5]
        elif sys.argv[2] in ['structured','iso-structured','fixed-stem','iso-fixed-stem'] and len(sys.argv) == 8:
            L = int(sys.argv[3])
            toemax = int(sys.argv[4])
            stemmax = int(sys.argv[5])
            N = int(sys.argv[6])
            filename = sys.argv[7]
        else:
            print """Usage:
              python -i hybridization_scatterplot generate random <len> <num seqs> <data file name>
              python -i hybridization_scatterplot generate iso-random <len> <num seqs> <data file name>
              python -i hybridization_scatterplot generate structured <len> <max toe> <max stem> <num seqs> <data file name>
              python -i hybridization_scatterplot generate iso-structured <len> <max toe> <max stem> <num seqs> <data file name>
              python -i hybridization_scatterplot generate fixed-stem <len> <max toe> <stem size> <num seqs> <data file name>
              python -i hybridization_scatterplot generate iso-fixed-stem <len> <max toe> <stem size> <num seqs> <data file name>
              python -i hybridization_scatterplot plot <data file names>
              """
            sys.exit()
        if not set(filename) <= set('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-'): # don't allow crazy symbols
            print "Data file postfix names must be alphanumeric only."
            sys.exit()

    if sys.argv[1] == 'plot':
        filenames = []
        for fn in sys.argv[2:]:
            if not set(fn) <= set('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-'): # don't allow crazy symbols
                print "Data file postfix names must be alphanumeric only."
                sys.exit()
            filenames.append(fn)

    print
    if sys.argv[1] == 'generate' and (sys.argv[2] in ['random','iso-random','structured','iso-structured','fixed-stem','iso-fixed-stem']):   
        
        if 'iso' in sys.argv[2]:
            Ngenerate = 100*N  # we'll generate many sequences and choose the ones with the closest-to-average duplex energy
            print "Generating %d sequences each of length %d, with similar duplex energies..." % (N, L)
        else:
            Ngenerate = N

        # first choose the sequences
        if 'random' in sys.argv[2]:
            candidates = [ randomseq(L,'ACTG') for i in range(Ngenerate) ]  
        elif 'structured' in sys.argv[2] or 'fixed-stem' in sys.argv[2]:    
            # generate sequences according to domain-level specifications with a range of hairpin stem and toehold lengths
            # (L,am,bm,cm) = (30,6,8,6)   # up to size 6+2*8+6 = 28 nt in the toeholds and stems.
            # (L,am,bm,cm) = (21,5,6,3)    # up to size 5+2*6+3 = 20 nt in the toeholds and stems.
            (am,bm,cm) = (toemax,stemmax,toemax)
            candidates = []
            for i in range(Ngenerate):
                d = 0
                while d <= 0:
                    a = random.randint(1,am) if sys.argv[2] == 'structured' else random.randint(min(2,am),am)
                    b = random.randint(1,bm) if sys.argv[2] == 'structured' else bm
                    c = random.randint(1,cm) if sys.argv[2] == 'structured' else random.randint(min(2,cm),cm)
                    d = L-a-2*b-c
                stem = randomseq(b)
                #seq  = randomseq(a,'ATC')+stem+randomseq(d,'T')+WC(stem)+randomseq(c,'ATC')   # tends to have different duplex_dG, which confounds analysis
                seq  = randomseq(a,'ATCG')+stem+randomseq(d,'ACTG')+WC(stem)+randomseq(c,'ATCG')
                candidates.append(seq)

        if 'iso' in sys.argv[2]:  
            candidates = list(set(candidates))  # get rid of duplicates
            dGs = [ duplex_dG(seq) for seq in candidates ]
            mean_dG = np.mean(dGs)
            dist_dG = [ abs(dG-mean_dG) for dG in dGs ]
            ranked_candidates = sorted( zip(candidates, dist_dG), key = lambda z: z[1])
            candidates = [ s for (s,d) in ranked_candidates[0:N] ]


        # then evaluate the sequences...
        results=[]
        for i in range(N):
            print
            seq  = candidates[i]
            print "%d: Simulating %s DNA strand at 25 C:  %s with toeholds %d and %d, stem %d, ss dG = %g and %g, ds dG = %g " % \
                    ( i,sys.argv[2],seq,toeholds(seq)[0],toeholds(seq)[1],stemsize(seq),strand_dG(seq),strand_dG(WC(seq)),duplex_dG(seq) )

            # Accumulate statistics from all the runs, until enough succesful simulations are collected.
            N_forward_total=0
            N_reverse_total=0
            total_time_forward=0.0
            total_time_reverse=0.0
            k1_net = 0.0
            k1prime_net = 0.0
            trials=30
            while N_forward_total < 25:   # this is a bit wasteful, but "only" by a factor of about 1.5.  Aims for 20% error bars on k1.
                print
                o = first_step_simulation(seq, trials, concentration=50e-9, T=25, material="DNA") 
                print "Inferred rate constants with analytical error bars (this batch only):"
                N_forward, N_reverse, kcoll, forward_kcoll, reverse_kcoll, k1, k2, k1prime, k2prime, keff, zcrit = compute_rate_constants(o.interface.results,50e-9)

                # update cumulative averages
                total_time_forward += float(N_forward)/k2 if N_forward > 0 else 0
                total_time_reverse += float(N_reverse)/k2prime if N_reverse > 0 else 0
                if N_forward_total > 0:
                    k1_net = (N_forward_total * k1_net + N_forward * k1) / (N_forward_total + N_forward) if N_forward > 0 else k1_net
                else:
                    k1_net = (N_forward * k1) / N_forward if N_forward > 0 else k1_net
                if N_reverse_total > 0:
                    k1prime_net = (N_reverse_total * k1prime_net + N_reverse * k1prime) / (N_reverse_total + N_reverse) if N_reverse > 0 else k1prime_net
                else:
                    k1prime_net = (N_reverse * k1prime) / N_reverse if N_reverse > 0 else k1prime_net
                N_forward_total += N_forward
                N_reverse_total += N_reverse

                if trials < 500000:  # memory problems running huge numbers of trials.  maybe we need to do better:  compress dataset incrementally...?
                    trials*=2

            # use cumulative averages to determine overall statistics
            k2_net = N_forward_total / total_time_forward
            k2prime_net = N_reverse_total / total_time_reverse
            print " %s ==> N_forward = %d, N_reverse = %d, k1 = %g, k1_prime = %g, k2 = %g, k2prime = %g" %  \
                (seq, N_forward_total, N_reverse_total, k1_net, k1prime_net, k2_net, k2prime_net)

            results.append([seq, N_forward_total, N_reverse_total, k1_net, k1prime_net, k2_net, k2prime_net])

        # now save the results for later...
        with open("scatterdata/" + filename + ".pkl", "wb") as f:
            pickle.dump(results, f)
        pdfname = filename

    elif sys.argv[1] == 'plot': 

        results = []
        for fn in filenames:
            print "Loading data from %s ..." % ("scatterdata/" + fn + ".pkl")
            with open("scatterdata/" + fn + ".pkl", "rb") as f:
                these_results = pickle.load(f)
            results += these_results
        pdfname = "+".join(filenames)
            
    else:
        print "Didn't understand args."
        sys.exit()

    # some old data files contain only (seq,kf); newer ones contains (seq, Nf, Nr, k1, k1p, k2, k2p); but here we only need kf == k1
    results = [ ( (data[0],data[3]) if len(data)==7 else (data[0],data[1]) ) for data in results ]  


    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages


    print
    print "Calculating stats on database of simulated strands..."
    toe_adj     = 2
    toes        = [ min(max(max(toeholds(seq)), sum(toeholds(seq))-toe_adj),17) for (seq,kf) in results ]  # effective length of toeholds
    stems       = [ stemsize(seq) for (seq,kf) in results ]                                                # length of duplex stem
    toe_dGs     = [ toebinding(seq) for (seq,kf) in results ]                               # (positive number, 1 or 2) for best energy of binding by 1 or 2 toeholds
    stem_dGs    = [ max( -strand_dG(seq), -strand_dG(WC(seq)) ) for (seq,kf) in results ]   # positive number, whichever is stronger
    binding_dGs = [ binding_dG(seq) for (seq,kf) in results ]
    duplex_dGs  = [ duplex_dG(seq) for (seq,kf) in results ]
    strand_dGs  = [ strand_dG(seq) + strand_dG(WC(seq)) for (seq,kf) in results ]   
    kfs         = [ kf for (seq,kf) in results ]
    krs         = [ reverse_rate(seq,kf) for (seq,kf) in results ]
    log_kfs     = [ np.log10(kf) for (seq,kf) in results]  # forward rates
    log_krs     = [ np.log10(kr) for kr in krs]            # reverse rates

    extremes = ''
    i = np.argmax(kfs)
    extremes += "\nFastest association: %s at %g /M/s" % (results[i][0],kfs[i])
    i = np.argmin(kfs)
    extremes += "\nSlowest association: %s at %g /M/s" % (results[i][0],kfs[i])
    i = np.argmax(krs)
    extremes += "\nFastest dissociation: %s at %g /s" % (results[i][0],krs[i])
    i = np.argmin(krs)
    extremes += "\nSlowest dissociation: %s at %g /s" % (results[i][0],krs[i])

    with PdfPages("scatterdata/"+pdfname+".pdf") as pdf:
        print "Drawing plots, saving to 'scatterdata/" + pdfname + ".pdf'..." 

        # Two subplots, the axes array is 1-d
        plt.figure(1)
        plt.subplots_adjust( hspace=1.0)
        plt.subplot(211)
        n, bins, patches = plt.hist(log_kfs, 20, normed=1, facecolor='green', alpha=0.75)
        plt.title("Association and dissociation rate distributions")
        plt.ylabel("frequency of\nassociation rates",fontsize='large')
        plt.yticks(fontsize='larger',va='bottom')
        plt.xlabel("Log10 rate constant kf (/M/s)" + extremes,fontsize='large')
        plt.xticks(fontsize='larger')
        plt.subplot(212)
        n, bins, patches = plt.hist(log_krs, 20, normed=1, facecolor='green', alpha=0.75)
        plt.ylabel("frequency of\ndissociation rates",fontsize='large')
        plt.yticks(fontsize='larger',va='bottom')
        plt.xlabel("Log10 rate constant kr (/s)",fontsize='large')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        # Do the rates depend upon how much secondary structure the single strands have?
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(strand_dGs, log_kfs, s=[10*s for s in stems], alpha=0.5)
        plt.title("Association and dissociation rates")
        plt.ylabel("Log10 rate constant kf (/M/s) \n",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.subplot(212)
        plt.scatter(strand_dGs, log_krs, s=[10*s for s in stems], alpha=0.5)
        plt.ylabel("Log10 rate constant kr (/s)",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.xlabel("single-strand secondary structure (kcal/mol)",fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        # Do the rates depend upon the overall dG of the reaction
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(binding_dGs, log_kfs, s=[10*s for s in stems], alpha=0.5)
        plt.title("Association and dissociation rates")
        plt.ylabel("Log10 rate constant kf (/M/s) \n",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.subplot(212)
        plt.scatter(binding_dGs, log_krs, s=[10*s for s in stems], alpha=0.5)
        plt.ylabel("Log10 rate constant kr (/s)",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.xlabel("net binding strength (kcal/mol)",fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        # Do the rates depend upon the strength of the duplex state?
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(duplex_dGs, log_kfs, s=[10*s for s in stems], alpha=0.5)
        plt.title("Association and dissociation rates")
        plt.ylabel("Log10 rate constant kf (/M/s) \n",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.subplot(212)
        plt.scatter(duplex_dGs, log_krs, s=[10*s for s in stems], alpha=0.5)
        plt.ylabel("Log10 rate constant kr (/s)",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.xlabel("duplex strength (kcal/mol)",fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        singletoes = [ max(toeholds(seq)) for (seq,kf) in results ]
        toecolors = ['red' if s==t else 'blue' for (s,t) in zip(singletoes,toes)]

        # Do the rates depend upon toehold lengths?
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(toes, log_kfs, s=[10*s for s in stems], color=toecolors, alpha=0.5)
        plt.title("Association and dissociation rates (double-toeholds are blue)")
        plt.ylabel("Log10 rate constant kf (/M/s) \n",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.subplot(212)
        plt.scatter(toes, log_krs, s=[10*s for s in stems], color=toecolors, alpha=0.5)
        plt.ylabel("Log10 rate constant kr (/s)",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.xlabel("effective toehold length (nt)",fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        toes        = np.array(toes)
        stems       = np.array(stems)
        toe_dGs     = np.array(toe_dGs)       # max of the MFE difference for toehold binding configurations (i.e. positive number)
        stem_dGs    = np.array(stem_dGs)      # max of the two (negated) ss pfuncs (i.e. positive number)
        binding_dGs = np.array(binding_dGs)
        duplex_dGs  = np.array(duplex_dGs)
        strand_dGs  = np.array(strand_dGs)    # sum of both ss pfuncs  (negative number)
        kfs         = np.array(kfs)
        krs         = np.array(krs)
        log_kfs     = np.array(log_kfs)
        log_krs     = np.array(log_krs)
        long = np.array([ s > 3 for s in stems ])   # numpy, when indexed by Boolean array, selects that subset of items
        short = np.array([ s < 3 for s in stems ])

        # Do the rates depend upon toehold lengths, for long-stemmed strands?
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(toes[long], log_kfs[long], s=[10*s for s in stems[long]], color=toecolors, alpha=0.5)
        plt.title("Association and dissociation rates for long-stem strands")
        plt.ylabel("Log10 rate constant kf (/M/s) \n",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.subplot(212)
        plt.scatter(toes[long], log_krs[long], s=[10*s for s in stems[long]], color=toecolors, alpha=0.5)
        plt.ylabel("Log10 rate constant kr (/s)",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.xlabel("effective toehold length (nt)",fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        # Do the rates depend upon toehold lengths, for short-stemmed strands?
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(toes[short], log_kfs[short], s=[10*s for s in stems[short]], color=toecolors, alpha=0.5)
        plt.title("Association and dissociation rates for short-stem strands")
        plt.ylabel("Log10 rate constant kf (/M/s) \n",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.subplot(212)
        plt.scatter(toes[short], log_krs[short], s=[10*s for s in stems[short]], color=toecolors, alpha=0.5)
        plt.ylabel("Log10 rate constant kr (/s)",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.xlabel("effective toehold length (nt)",fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        # Do the rates depend upon toehold strengths?
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(toe_dGs, log_kfs, s=[10*s for s in stems], color=toecolors, alpha=0.5)
        plt.title("Association and dissociation rates (double-toeholds are blue)")
        plt.ylabel("Log10 rate constant kf (/M/s) \n",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.subplot(212)
        plt.scatter(toe_dGs, log_krs, s=[10*s for s in stems], color=toecolors, alpha=0.5)
        plt.ylabel("Log10 rate constant kr (/s)",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.xlabel("toehold binding strength (kcal/mol)",fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()


        r1 = np.array( [ (random.random()-0.5) / 2.0 for i in range(len(toes)) ] )  # jitter to avoid super-positioning
        r2 = np.array( [ (random.random()-0.5) / 2.0 for i in range(len(toes)) ] )  # jitter to avoid super-positioning

        # Not quite a contour plot, but can we illustrate kf as a function of toehold length and stem length
        plt.figure(1)
        plt.subplot(111)
        plt.scatter(toes+r1, stems+r2, s=50, c = log_kfs, alpha=0.5)
        # plt.gray()
        cb=plt.colorbar()
        cb.set_label('log10 of kf')
        plt.title("Association rates in toehold-length vs stem-length space\n")
        plt.ylabel("stem length",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.xlabel("effective toehold length",fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        plt.figure(1)
        plt.subplot(111)
        plt.scatter(toes+r1, stems+r2, s=50, c = log_krs, alpha=0.5)
        # plt.gray()
        cb=plt.colorbar()
        cb.set_label('log10 of kr')
        plt.title("Disassociation rates in toehold-length vs stem-length space\n")
        plt.ylabel("stem length",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.xlabel("effective toehold length",fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()


        # Not quite a contour plot, but can we illustrate kf as a function of toehold strength and stem strength
        plt.figure(1)
        plt.subplot(111)
        plt.scatter(toe_dGs, stem_dGs, s=50, c = log_kfs, alpha=0.5)
        # plt.gray()
        cb=plt.colorbar()
        cb.set_label('log10 of kf')
        plt.title("Association rates in toehold-strength vs stem-strength space\n")
        plt.ylabel("stem strength",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.xlabel("toehold strength",fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        plt.figure(1)
        plt.subplot(111)
        plt.scatter(toe_dGs, stem_dGs, s=50, c = log_krs, alpha=0.5)
        # plt.gray()
        cb=plt.colorbar()
        cb.set_label('log10 of kr')
        plt.title("Disassociation rates in toehold-strength vs stem-strength space\n")
        plt.ylabel("stem strength",fontsize='larger')
        plt.yticks(fontsize='larger',va='bottom')
        plt.xlabel("toehold strength",fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        number_distinct_lengths = len(set([len(seq) for (seq,kf) in results]))
        print "duplex dG standard deviation = %g" % np.std(duplex_dGs)
        if number_distinct_lengths == 1 and np.std(duplex_dGs) < .2:
            log_fastest_assoc = np.max(log_kfs)
            log_slowest_dissoc = np.min(log_krs)
            dissoc_speed_ups = log_krs - log_slowest_dissoc
            assoc_slow_downs = log_fastest_assoc - log_kfs
            updown_range = round(max(max(dissoc_speed_ups),max(assoc_slow_downs)))

            # Does secondary structure primarily speed up dissociation rather than slow down association?
            plt.figure(1)
            plt.subplot(111)
            plt.scatter(dissoc_speed_ups, assoc_slow_downs, s=20, alpha=0.5)
            plt.plot([0, updown_range],[0,updown_range],'r--')
            plt.title("Secondary structure slows down association \n or speeds up dissociation?")
            plt.ylabel("Log10 slow-down of kf (relative to fastest)",fontsize='larger')
            plt.yticks(fontsize='larger',va='bottom')
            plt.xlabel("Log10 speed-up of kr (relative to slowest)",fontsize='larger')
            plt.xticks(fontsize='larger')
            pdf.savefig()
            plt.close()

            # Does secondary structure primarily speed up dissociation rather than slow down association?
            plt.figure(1)
            plt.subplot(111)
            plt.scatter(dissoc_speed_ups, assoc_slow_downs, s=20, c = stem_dGs, alpha=0.5)
            plt.plot([0, updown_range],[0,updown_range],'r--')
            cb=plt.colorbar()
            cb.set_label('stem strength (kcal/mol)')
            plt.title("Secondary structure slows down association \n or speeds up dissociation?")
            plt.ylabel("Log10 slow-down of kf (relative to fastest)",fontsize='larger')
            plt.yticks(fontsize='larger',va='bottom')
            plt.xlabel("Log10 speed-up of kr (relative to slowest)",fontsize='larger')
            plt.xticks(fontsize='larger')
            pdf.savefig()
            plt.close()

            # Does secondary structure primarily speed up dissociation rather than slow down association?
            plt.figure(1)
            plt.subplot(111)
            plt.scatter(dissoc_speed_ups, assoc_slow_downs, s=20, c = toe_dGs, alpha=0.5)
            plt.plot([0, updown_range],[0,updown_range],'r--')
            cb=plt.colorbar()
            cb.set_label('toehold strength (kcal/mol)')
            plt.title("Secondary structure slows down association \n or speeds up dissociation?")
            plt.ylabel("Log10 slow-down of kf (relative to fastest)",fontsize='larger')
            plt.yticks(fontsize='larger',va='bottom')
            plt.xlabel("Log10 speed-up of kr (relative to slowest)",fontsize='larger')
            plt.xticks(fontsize='larger')
            pdf.savefig()
            plt.close()

            # Does secondary structure primarily speed up dissociation rather than slow down association?
            plt.figure(1)
            plt.subplot(111)
            plt.scatter(dissoc_speed_ups, assoc_slow_downs, s=20, c = stem_dGs - toe_dGs, alpha=0.5)
            plt.plot([0, updown_range],[0,updown_range],'r--')
            cb=plt.colorbar()
            cb.set_label('stem strength - toehold strength (kcal/mol)')
            plt.title("Secondary structure slows down association \n or speeds up dissociation?")
            plt.ylabel("Log10 slow-down of kf (relative to fastest)",fontsize='larger')
            plt.yticks(fontsize='larger',va='bottom')
            plt.xlabel("Log10 speed-up of kr (relative to slowest)",fontsize='larger')
            plt.xticks(fontsize='larger')
            pdf.savefig()
            plt.close()
