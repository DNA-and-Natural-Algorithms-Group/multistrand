# hybridization_transitions.py
# 
# This example is modified from threeway_transitions.py, but provides a deeper analysis approach (for "simple" hybridization).
# The two strands are constrained to have embedded hairpin-loop secondary structure, as in Gao et al 2006.
# Simulation proceeds from completely unfolded strands separate in the given volume, completing when the full hybrid is formed.
# The full set of domain-level secondary structures (all 44 of them) are used as macrostates, analogous to KinDA's approach.
# With enough simulations, rates between macrostates should be sufficient to estimate both the overall k_f and k_r rate constants,
# by numerically computing expected hitting times from initial to final states within the macrostate continuous-time markov chain,
# and from the final to the initial states as well. This is somewhat reminiscent of forward flux sampling.
#
# In addition to obtaining rates, this analysis ought to allow us to make claims about the pathway taken, e.g. 
#   "Did hybridization initiate via hairpins or unfolded strands?"
#   "Was 3-way or 4-way branch migration involved, or neither?"
#
# Try it like this, e.g.:
#   python -i hybridization_transitions.py medium count 100
#   python -i hybridization_transitions.py medium count2 100
#   python -i hybridization_transitions.py medium loose 100
#   python -i hybridization_transitions.py medium exact 100
#   python -i hybridization_transitions.py medium loose2 1000
#   python -i hybridization_transitions.py medium distance 100
#   python -i hybridization_transitions.py medium distance2 300
# where the second argument is the number of trajectories to simulate, and the first argument indicates the coarse-graining method to be used.
# The postfix "2" requests a second-order markov chain approximation, which entails a larger state space and thus requires more simulation trials.
#
# For comparison to first-step mode, for the 'medium' sequence, try:
#   python -i hybridization_first_step_mode.py DNA ACTCCACATTGTCTTTATATTGTGACAATGAAGCGT 25 1e-3 1000
# and the expected time to full duplex, including time wasted on failed collisions, is
#   1/(k_eff * c)
# where c is the concentration 1e-3 (i.e. 1 mM).  
#
# When the reverse pathway hitting time can be computed, we can also compare the dissociation rate inferred from first-step mode + thermodynamics:
#   kr = kf exp(dG/RT)
# where using NUPACK we can calculate
#   dG = pfunc(ds) - pfunc(ss1) - pfunc(ss2)
# and thus the expected time for dissociation is
#   (1/kf) exp(-dG/RT)
# 
# The "small", "medium", and "large" sequences can each be evaluated similarly
#   python -i hybridization_first_step_mode.py DNA TCTACTTCTACTATAC 25 1e-3 1000
#      kf = k1 = 77.7 x 10^6 /M/s     k_eff =  23.6 x 10^6 /M/s      1/(k_eff * c ) =  42.4 usec
#      dG = -18.76 - (-0.22 + -0.10) = -18.44  kcal/mol
#      kr = kf exp(-18.44 / 0.5924) ==> 8.43 days
#   python -i hybridization_first_step_mode.py DNA ACTCCACATTGTCTTTATATTGTGACAATGAAGCGT 25 1e-3 1000
#      kf = k1 = 7.4 x 10^6 /M/s     k_eff = 7.1 x 10^6 /M/s      1/(k_eff * c ) = 141 usec
#      dG = -52.87 - (-8.45 + -8.33) = -36.10  kcal/mol
#      kr = kf exp(-36.10 / 0.5924) ==> 4.55 x 10^14 days
#   python -i hybridization_first_step_mode.py DNA CATTGTAACTGGCGATGCTACCTGTATTTTTACAGGTAGCATCGCCCCATTAACTC 25 1e-3 1000
#      kf = k1 = 54.2 x 10^6 /M/s     k_eff =  8.6 x 10^6 /M/s      1/(k_eff * c ) =  116 usec
#      dG = -86.42 - (-25.08 + -24.72) = -36.6  kcal/mol
#      kr = kf exp(-36.6 / 0.5924) ==>  1.4 x 10^14 days
#


# To-do note: must describe macrostate approaches and first-order vs second-order chains


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
import sys

#############

# for StopCondition and Macrostate definitions:
Exact_Macrostate = 0   # match a secondary structure exactly (i.e. any system state that has a complex with this exact structure)
Bound_Macrostate = 1   # match any system state in which the given strand is bound to another strand
Dissoc_Macrostate = 2  # match any system state in which there exists a complex with exactly the given strands, in that order
Loose_Macrostate = 3   # match a secondary structure with "don't care"s, allowing a certain number of disagreements
Count_Macrostate = 4   # match a secondary structure, allowing a certain number of disagreements
# see Schaeffer's PhD thesis, chapter 7.2, for more information


# Here, in one of two approaches illustrated in this file, we define macrostates for the state of the tube by explicitly describing 
# each perfect-domain-level-matched secondary structure, and allowing a fixed number of incorrectly-paired nucleotides.  
# This requires enumerating all 44 such domain-level secondary structures.  The Count macrostate type is used.
# The resulting macrostates seldom overlap, resulting in a coarse-grained model roughly the same number of states (e.g. around 24).
# However, this doesn't end up providing a detailed picture of the concluding steps of hybridization, 
# from the four-way junction to its resolution as a complete duplex.
# Consequently, while the overall forward rate (hybridization) can be calculated, the reverse rate (dissociation) cannot.
# Regarding estimation of the dissociation rate, note that the step from the exact perfect duplex to the nearest macrostate is 
# (a) a large step, involving the fraying or breathing of a full domain, and 
# (b) not measured, because simulation stops once that state is reached.


def setup_options_hybridization_Count(trials, toe5_seq, stem_seq, loop_seq, toe3_seq, temp=298.15, conc=1e-3, time=0.0002):

    # I am not sure whether 0-length domains will break Multistrand.  Don't take chances for now...
    assert len(toe5_seq)>0
    assert len(stem_seq)>0
    assert len(loop_seq)>0
    assert len(toe3_seq)>0
    # Thus, to approximate a zero-toehold case, use 1 nt for both toeholds, with complementary sequences... 

    L = len(toe5_seq+stem_seq+loop_seq+stem_seq+toe3_seq)

    toe5 = Domain(name="toe5",sequence=toe5_seq,length=len(toe5_seq))
    stem = Domain(name="stem",sequence=stem_seq,length=len(stem_seq))
    loop = Domain(name="loop",sequence=loop_seq,length=len(loop_seq))
    toe3 = Domain(name="toe3",sequence=toe3_seq,length=len(toe3_seq))
    
    sense_strand     = toe5+stem+loop+stem.C+toe3
    antisense_strand = toe3.C+stem+loop.C+stem.C+toe5.C

    # start with individual strands folded as hairpins, end with complete hybrid
    sense_folded       = Complex(strands=[sense_strand], structure=".(.).")
    antisense_folded   = Complex(strands=[antisense_strand], structure=".(.).")
    hybrid_complete    = Complex(strands=[sense_strand,antisense_strand], structure="(((((+)))))")

    # define the intermediates
    sense_unfolded     = Complex(strands=[sense_strand], structure=".....")
    antisense_unfolded = Complex(strands=[antisense_strand], structure=".....")

       # domain notation ABDbC and cBdba   ..... 39 structures plus the above 5
    complex_A       = Complex(strands=[sense_strand,antisense_strand], structure="(....+....)")
    complex_B1      = Complex(strands=[sense_strand,antisense_strand], structure=".(...+...).")
    complex_D       = Complex(strands=[sense_strand,antisense_strand], structure="..(..+..)..")
    complex_B2      = Complex(strands=[sense_strand,antisense_strand], structure="...(.+.)...")
    complex_C       = Complex(strands=[sense_strand,antisense_strand], structure="....(+)....")

    complex_Ahp_s   = Complex(strands=[sense_strand,antisense_strand], structure="((.).+....)")
    complex_Ahp_a   = Complex(strands=[sense_strand,antisense_strand], structure="(....+.(.))")
    complex_Ahp_as  = Complex(strands=[sense_strand,antisense_strand], structure="((.).+.(.))")
    complex_Chp_s   = Complex(strands=[sense_strand,antisense_strand], structure=".(.)(+)....")
    complex_Chp_a   = Complex(strands=[sense_strand,antisense_strand], structure="....(+)(.).")
    complex_Chp_as  = Complex(strands=[sense_strand,antisense_strand], structure=".(.)(+)(.).")

    complex_AChp_s  = Complex(strands=[sense_strand,antisense_strand], structure="((.)(+)...)")
    complex_AChp_a  = Complex(strands=[sense_strand,antisense_strand], structure="(...(+)(.))")
    complex_AChp_as = Complex(strands=[sense_strand,antisense_strand], structure="((.)(+)(.))")

    complex_AB1     = Complex(strands=[sense_strand,antisense_strand], structure="((...+...))")
    complex_AD      = Complex(strands=[sense_strand,antisense_strand], structure="(.(..+..).)")
    complex_AB2     = Complex(strands=[sense_strand,antisense_strand], structure="(..(.+.)..)")
    complex_AC      = Complex(strands=[sense_strand,antisense_strand], structure="(...(+)...)")
    complex_B1D     = Complex(strands=[sense_strand,antisense_strand], structure=".((..+..)).")
    complex_B1B2    = Complex(strands=[sense_strand,antisense_strand], structure=".(.(.+.).).")
    complex_B1C     = Complex(strands=[sense_strand,antisense_strand], structure=".(..(+)..).")
    complex_DB2     = Complex(strands=[sense_strand,antisense_strand], structure="..((.+.))..")
    complex_DC      = Complex(strands=[sense_strand,antisense_strand], structure="..(.(+).)..")
    complex_B2C     = Complex(strands=[sense_strand,antisense_strand], structure="...((+))...")

    complex_AB1D    = Complex(strands=[sense_strand,antisense_strand], structure="(((..+..)))")
    complex_AB1B2   = Complex(strands=[sense_strand,antisense_strand], structure="((.(.+.).))")
    complex_AB1C    = Complex(strands=[sense_strand,antisense_strand], structure="((..(+)..))")
    complex_ADB2    = Complex(strands=[sense_strand,antisense_strand], structure="(.((.+.)).)")
    complex_ADC     = Complex(strands=[sense_strand,antisense_strand], structure="(.(.(+).).)")
    complex_AB2C    = Complex(strands=[sense_strand,antisense_strand], structure="(..((+))..)")
    complex_B1DB2   = Complex(strands=[sense_strand,antisense_strand], structure=".(((.+.))).")
    complex_B1DC    = Complex(strands=[sense_strand,antisense_strand], structure=".((.(+).)).")
    complex_B1B2C   = Complex(strands=[sense_strand,antisense_strand], structure=".(.((+)).).")
    complex_DB2C    = Complex(strands=[sense_strand,antisense_strand], structure="..(((+)))..")

    complex_AB1DB2  = Complex(strands=[sense_strand,antisense_strand], structure="((((.+.))))")
    complex_AB1DC   = Complex(strands=[sense_strand,antisense_strand], structure="(((.(+).)))")
    complex_AB1B2C  = Complex(strands=[sense_strand,antisense_strand], structure="((.((+)).))")
    complex_ADB2C   = Complex(strands=[sense_strand,antisense_strand], structure="(.(((+))).)")
    complex_B1DB2C  = Complex(strands=[sense_strand,antisense_strand], structure=".((((+)))).")
    complex_AB1DB2C = Complex(strands=[sense_strand,antisense_strand], structure="(((((+)))))")

    ## The above could be done equivalently with combinatorial sets of loose macrostates for each domain.  (See "Loose" setup below.)
    ## KinDA defines macrostates in terms of each domain being x% correct; we'd should do that here, but we take the shortcut.
    
    incorrect_frac = 0.1
    d1 = int(L*incorrect_frac)
    d2 = int(2*L*incorrect_frac)

    print("Macrostates: SS structures allowing up to %d nt incorrect, DS structures allowing up to %d nt incorrect." % (d1, d2))

    macro_U1_U2           = Macrostate("U1_and_U2", [(sense_unfolded,Count_Macrostate,d1),(antisense_unfolded,Count_Macrostate,d1)])          
    macro_U1_F2           = Macrostate("U1_and_F2", [(sense_unfolded,Count_Macrostate,d1),(antisense_folded,Count_Macrostate,d1)])
    macro_F1_U2           = Macrostate("F1_and_U2", [(sense_folded,Count_Macrostate,d1),(antisense_unfolded,Count_Macrostate,d1)])          
    macro_F1_F2           = Macrostate("F1_and_F2", [(sense_folded,Count_Macrostate,d1),(antisense_folded,Count_Macrostate,d1)])

    macro_A               = Macrostate("COMPLEX_A", [(complex_A,Count_Macrostate,d2)])
    macro_B1              = Macrostate("COMPLEX_B1", [(complex_B1,Count_Macrostate,d2)])
    macro_D               = Macrostate("COMPLEX_D", [(complex_D,Count_Macrostate,d2)])
    macro_B2              = Macrostate("COMPLEX_B2", [(complex_B2,Count_Macrostate,d2)])
    macro_C               = Macrostate("COMPLEX_C", [(complex_C,Count_Macrostate,d2)])
    macro_Ahp_s           = Macrostate("COMPLEX_Ahp_s", [(complex_Ahp_s,Count_Macrostate,d2)])
    macro_Ahp_a           = Macrostate("COMPLEX_Ahp_a", [(complex_Ahp_a,Count_Macrostate,d2)])
    macro_Ahp_as          = Macrostate("COMPLEX_Ahp_as", [(complex_Ahp_as,Count_Macrostate,d2)])
    macro_Chp_s           = Macrostate("COMPLEX_Chp_s", [(complex_Chp_s,Count_Macrostate,d2)])
    macro_Chp_a           = Macrostate("COMPLEX_Chp_a", [(complex_Chp_a,Count_Macrostate,d2)])
    macro_Chp_as          = Macrostate("COMPLEX_Chp_as", [(complex_Chp_as,Count_Macrostate,d2)])
    macro_AChp_s          = Macrostate("COMPLEX_AChp_s", [(complex_AChp_s,Count_Macrostate,d2)])
    macro_AChp_a          = Macrostate("COMPLEX_AChp_a", [(complex_AChp_a,Count_Macrostate,d2)])
    macro_AChp_as         = Macrostate("COMPLEX_AChp_as", [(complex_AChp_as,Count_Macrostate,d2)])
    macro_AB1             = Macrostate("COMPLEX_AB1", [(complex_AB1,Count_Macrostate,d2)])
    macro_AD              = Macrostate("COMPLEX_AD", [(complex_AD,Count_Macrostate,d2)])
    macro_AB2             = Macrostate("COMPLEX_AB2", [(complex_AB2,Count_Macrostate,d2)])
    macro_AC              = Macrostate("COMPLEX_AC", [(complex_AC,Count_Macrostate,d2)])
    macro_B1D             = Macrostate("COMPLEX_B1D", [(complex_B1D,Count_Macrostate,d2)])
    macro_B1B2            = Macrostate("COMPLEX_B1B2", [(complex_B1B2,Count_Macrostate,d2)])
    macro_B1C             = Macrostate("COMPLEX_B1C", [(complex_B1C,Count_Macrostate,d2)])
    macro_DB2             = Macrostate("COMPLEX_DB2", [(complex_DB2,Count_Macrostate,d2)])
    macro_DC              = Macrostate("COMPLEX_DC", [(complex_DC,Count_Macrostate,d2)])
    macro_B2C             = Macrostate("COMPLEX_B2C", [(complex_B2C,Count_Macrostate,d2)])
    macro_AB1D            = Macrostate("COMPLEX_AB1D", [(complex_AB1D,Count_Macrostate,d2)])
    macro_AB1B2           = Macrostate("COMPLEX_AB1B2", [(complex_AB1B2,Count_Macrostate,d2)])
    macro_AB1C            = Macrostate("COMPLEX_AB1C", [(complex_AB1C,Count_Macrostate,d2)])
    macro_ADB2            = Macrostate("COMPLEX_ADB2", [(complex_ADB2,Count_Macrostate,d2)])
    macro_ADC             = Macrostate("COMPLEX_ADC", [(complex_ADC,Count_Macrostate,d2)])
    macro_AB2C            = Macrostate("COMPLEX_AB2C", [(complex_AB2C,Count_Macrostate,d2)])
    macro_B1DB2           = Macrostate("COMPLEX_B1DB2", [(complex_B1DB2,Count_Macrostate,d2)])
    macro_B1DC            = Macrostate("COMPLEX_B1DC", [(complex_B1DC,Count_Macrostate,d2)])
    macro_B1B2C           = Macrostate("COMPLEX_B1B2C", [(complex_B1B2C,Count_Macrostate,d2)])
    macro_DB2C            = Macrostate("COMPLEX_DB2C", [(complex_DB2C,Count_Macrostate,d2)])
    macro_AB1DB2          = Macrostate("COMPLEX_AB1DB2", [(complex_AB1DB2,Count_Macrostate,d2)])
    macro_AB1DC           = Macrostate("COMPLEX_AB1DC", [(complex_AB1DC,Count_Macrostate,d2)])
    macro_AB1B2C          = Macrostate("COMPLEX_AB1B2C", [(complex_AB1B2C,Count_Macrostate,d2)])
    macro_ADB2C           = Macrostate("COMPLEX_ADB2C", [(complex_ADB2C,Count_Macrostate,d2)])
    macro_B1DB2C          = Macrostate("COMPLEX_B1DB2C", [(complex_B1DB2C,Count_Macrostate,d2)])
    macro_AB1DB2C         = Macrostate("COMPLEX_AB1DB2C", [(complex_AB1DB2C,Count_Macrostate,d2)])

    initial_hairpins      = Macrostate("HAIRPINS", [(sense_folded,Exact_Macrostate,0),(antisense_folded,Exact_Macrostate,0)])   # perfect hairpins achieved
    completed_hybrid      = Macrostate("COMPLETE", [(hybrid_complete,Exact_Macrostate,0)])                                      # perfect hybrid achieved

    o = Options(simulation_mode="Transition", parameter_type="Nupack", dangles="Some",
                rate_method = "Metropolis", num_simulations = trials, simulation_time=time,  # default 0.0002 sec
                substrate_type="DNA", temperature=temp, join_concentration=conc,           # default 1 mM at 25C
                start_state=[sense_folded,antisense_folded], rate_scaling='Calibrated', verbosity=0)

    # for transition mode, these are mostly macrostates not stop conditions, since the simulation won't stop unless the name begins with "stop"
    o.stop_conditions = [initial_hairpins, macro_U1_U2, macro_U1_F2, macro_F1_U2, macro_F1_F2,                         
                macro_A, macro_B1, macro_D, macro_B2, macro_C, 
                macro_Ahp_s, macro_Ahp_a, macro_Ahp_as, macro_Chp_s, macro_Chp_a, macro_Chp_as,
                macro_AChp_s, macro_AChp_a, macro_AChp_as,
                macro_AB1, macro_AD, macro_AB2, macro_AC, macro_B1D, macro_B1B2, macro_B1C, macro_DB2, macro_DC, macro_B2C,
                macro_AB1D, macro_AB1B2, macro_AB1C, macro_ADB2, macro_ADC, macro_AB2C, macro_B1DB2, macro_B1DC, macro_B1B2C, macro_DB2C,
                macro_AB1DB2, macro_AB1DC, macro_AB1B2C, macro_ADB2C, macro_B1DB2C, macro_AB1DB2C, completed_hybrid]
    return o

# Here, in the second of four approaches illustrated in this file, we define macrostates for the state of the tube as combinatorial sets
# of macrostates that each define the presence of a single domain, allowing a fixed fraction of incorrectly-paired nucleotides.  
# Using the Loose macrostate type, this requires only 11 macrostate definitions (in addition to the exact start state and exact final state).
# Because only a fraction of the 2^11 subsets of (possibly intersecting) macrostates are actually reached,
# the overall coarse-grained system state space is only somewhat larger than in the Count macrostate approach (e.g. around 41 states).
# Pleasantly, this provides a more detailed view of the end stage of the hybridization process.
# Further, if the final state is used, rather than the perfect duplex, as the starting point,
# then the dissociate rate can be estimated, since it is possible to collect enough transitions for reverse steps.
# Unfortunately, this estimate appears to be poor when creating a first-order markov chain of macrostate transitions.
# A major problem is that since the macrostate configurations (membership subsets of the defined macrostates) are immediately contiguous,
# a single-base-pair move can transition directly from one to another -- and back!  So reverse transitions are common.
# Second-order chains (which remember the previous macrostate configuration) improve the estimate considerably, but are still suspect.

def setup_options_hybridization_Loose(trials, toe5_seq, stem_seq, loop_seq, toe3_seq, temp=298.15, conc=1e-3, time=0.0002):

    # I am not sure whether 0-length domains will break Multistrand.  Don't take chances for now...
    assert len(toe5_seq)>0
    assert len(stem_seq)>0
    assert len(loop_seq)>0
    assert len(toe3_seq)>0
    # Thus, to approximate a zero-toehold case, use 1 nt for both toeholds, with complementary sequences... 

    L = len(toe5_seq+stem_seq+loop_seq+stem_seq+toe3_seq)
    
    toe5 = Domain(name="toe5",sequence=toe5_seq,length=len(toe5_seq))
    stem = Domain(name="stem",sequence=stem_seq,length=len(stem_seq))
    loop = Domain(name="loop",sequence=loop_seq,length=len(loop_seq))
    toe3 = Domain(name="toe3",sequence=toe3_seq,length=len(toe3_seq))
    
    sense_strand     = toe5+stem+loop+stem.C+toe3
    antisense_strand = toe3.C+stem+loop.C+stem.C+toe5.C

    # start with individual strands folded as hairpins, end with complete hybrid
    sense_folded     = Complex(strands=[sense_strand], structure=".(.).")
    antisense_folded = Complex(strands=[antisense_strand], structure=".(.).")
    hybrid_complete  = Complex(strands=[sense_strand,antisense_strand], structure="(((((+)))))")

    # define the Loose base pairing patterns, with * for "don't care".
    loose_unfolded1 = Complex(strands=[sense_strand], structure="*.*.*")
    loose_unfolded2 = Complex(strands=[antisense_strand], structure="*.*.*")
    loose_folded1   = Complex(strands=[sense_strand], structure="*(*)*")
    loose_folded2   = Complex(strands=[antisense_strand], structure="*(*)*")
         # domain notation ABDbC and cBdba
    loose_A         = Complex(strands=[sense_strand,antisense_strand], structure="(****+****)")
    loose_B1        = Complex(strands=[sense_strand,antisense_strand], structure="*(***+***)*")
    loose_D         = Complex(strands=[sense_strand,antisense_strand], structure="**(**+**)**")
    loose_B2        = Complex(strands=[sense_strand,antisense_strand], structure="***(*+*)***")
    loose_C         = Complex(strands=[sense_strand,antisense_strand], structure="****(+)****")
    loose_H1        = Complex(strands=[sense_strand,antisense_strand], structure="*(*)*+*****")
    loose_H2        = Complex(strands=[sense_strand,antisense_strand], structure="*****+*(*)*")

    incorrect_frac = 0.5
    d_A = int(len(toe5_seq)*incorrect_frac)
    d_B = int(2*len(stem_seq)*incorrect_frac)
    d_C = int(len(toe3_seq)*incorrect_frac)
    d_D = int(len(loop_seq)*incorrect_frac)

    print("Macrostates: each domain allowing up to fraction %4.2g of nucleotides to be incorrect." % (incorrect_frac))

    macro_U1              = Macrostate("UNFOLDED1", [(loose_unfolded1,Loose_Macrostate,d_B)])
    macro_U2              = Macrostate("UNFOLDED2", [(loose_unfolded2,Loose_Macrostate,d_B)])
    macro_F1              = Macrostate("FOLDED1", [(loose_folded1,Loose_Macrostate,d_B)])
    macro_F2              = Macrostate("FOLDED2", [(loose_folded2,Loose_Macrostate,d_B)])

    macro_A               = Macrostate("DOMAIN_A", [(loose_A,Loose_Macrostate,d_A)])
    macro_B1              = Macrostate("DOMAIN_B1", [(loose_B1,Loose_Macrostate,d_B)])
    macro_D               = Macrostate("DOMAIN_D", [(loose_D,Loose_Macrostate,d_D)])
    macro_B2              = Macrostate("DOMAIN_B2", [(loose_B2,Loose_Macrostate,d_B)])
    macro_C               = Macrostate("DOMAIN_C", [(loose_C,Loose_Macrostate,d_C)])
    macro_H1              = Macrostate("DOMAIN_H1", [(loose_H1,Loose_Macrostate,d_B)])
    macro_H2              = Macrostate("DOMAIN_H2", [(loose_H2,Loose_Macrostate,d_B)])

    initial_hairpins      = Macrostate("HAIRPINS", [(sense_folded,Exact_Macrostate,0),(antisense_folded,Exact_Macrostate,0)])   # perfect hairpins achieved
    completed_hybrid      = Macrostate("COMPLETE", [(hybrid_complete,Exact_Macrostate,0)])                                      # perfect hybrid achieved

    o = Options(simulation_mode="Transition", parameter_type="Nupack", dangles="Some",
                rate_method = "Metropolis", num_simulations = trials, simulation_time=time,       # default 0.0002 sec
                substrate_type="DNA", temperature=temp, join_concentration=conc,                  # default 1 mM at 25C
                start_state=[sense_folded,antisense_folded], rate_scaling='Calibrated', verbosity=0)

    # for transition mode, these are macrostates not stop conditions, since the simulation won't stop unless the name begins with "stop"
    o.stop_conditions = [initial_hairpins, macro_U1, macro_U2, macro_F1, macro_F2, macro_A, macro_B1, macro_D, macro_B2, macro_C, macro_H1, macro_H2, completed_hybrid]
    return o

# Here there are two exact macrostates (the initial separate hairpins, and the final perfect duplex) and the rest measure distances from these two.
# All systems states containing separate strands can be measured as some distance from the initial hairpins; 
# those containing a single complex can be measured as some distance from the final duplex.
# Thus, the coarse-graining consists of a (quasi) linear set of states: linear branches each for separate strands and for complexes, 
# with potential transitions between all separate-strand macrostates and "early" (i.e. far-from-final) single-complex macrostates.
# If large enough distances are considered, then no microstates are excluded; transitions are exclusively between macrostate subsets.
# The intervals of distances considered -- here, 3 -- determine the coarseness; naturally, coarser macrostates make it easier to collect statistics.
# It also reduces the stored trajectory length.  With interval 1, every single step corresponds to a coarse-grained transition.
# The "Distance" coarse-graining appears to be best for estimating reverse-pathway rates (as long as second-order chains are used), 
# but it cannot provide much information about what type of pathway (e.g. 3-way vs 4-way branch migration) dominates reactions.

def setup_options_hybridization_Distance(trials, toe5_seq, stem_seq, loop_seq, toe3_seq, temp=298.15, conc=1e-3, time=0.0002):

    # I am not sure whether 0-length domains will break Multistrand.  Don't take chances for now...
    assert len(toe5_seq)>0
    assert len(stem_seq)>0
    assert len(loop_seq)>0
    assert len(toe3_seq)>0
    # Thus, to approximate a zero-toehold case, use 1 nt for both toeholds, with complementary sequences... 

    L = len(toe5_seq+stem_seq+loop_seq+stem_seq+toe3_seq)
    
    toe5 = Domain(name="toe5",sequence=toe5_seq,length=len(toe5_seq))
    stem = Domain(name="stem",sequence=stem_seq,length=len(stem_seq))
    loop = Domain(name="loop",sequence=loop_seq,length=len(loop_seq))
    toe3 = Domain(name="toe3",sequence=toe3_seq,length=len(toe3_seq))
    
    sense_strand     = toe5+stem+loop+stem.C+toe3
    antisense_strand = toe3.C+stem+loop.C+stem.C+toe5.C

    # start with individual strands folded as hairpins, end with complete hybrid
    sense_folded     = Complex(strands=[sense_strand], structure=".(.).")
    antisense_folded = Complex(strands=[antisense_strand], structure=".(.).")
    hybrid_complete  = Complex(strands=[sense_strand,antisense_strand], structure="(((((+)))))")

    stepsize=3
    print("Macrostates: nucleotides-incorrect distances from the initial and final states in steps of %d." % stepsize)

    macros = []
    for d in range(stepsize,2*L+1,stepsize):
        macros.append( Macrostate("I"+str(d), [(sense_folded,Count_Macrostate,d),(antisense_folded,Count_Macrostate,d)]) )
    for d in range(stepsize,2*L+1,stepsize):
        macros.append( Macrostate("F"+str(d), [(hybrid_complete,Count_Macrostate,d)]) )

    initial_hairpins      = Macrostate("HAIRPINS", [(sense_folded,Exact_Macrostate,0),(antisense_folded,Exact_Macrostate,0)])   # perfect hairpins achieved
    completed_hybrid      = Macrostate("COMPLETE", [(hybrid_complete,Exact_Macrostate,0)])                                      # perfect hybrid achieved

    o = Options(simulation_mode="Transition", parameter_type="Nupack", dangles="Some",
                rate_method = "Metropolis", num_simulations = trials, simulation_time=time,  # default 0.0002 sec
                substrate_type="DNA", temperature=temp, join_concentration=conc,             # default 1 mM at 25C
                start_state=[sense_folded,antisense_folded], rate_scaling='Calibrated', verbosity=0)

    # for transition mode, these are mostly macrostates not stop conditions, since the simulation won't stop unless the name begins with "stop"
    o.stop_conditions = [initial_hairpins]+macros+[completed_hybrid]
    return o

# In this fourth approach to macrostate definition, we chose a set of exact microstates along some canonical pathways.
# Being exact, the markov property for the coarse-grained chain is guaranteed, and second-order chains are unnecessary.
# However, because in principle trajectories can go for a long time between hitting any particular exact state,
# the transition matrix is not likely to be sparse -- and accuracy requires a lot of simulations.
# Furthermore, if (nearly) complete idealized trajectories are used, these "strings in state space" consist of a lot of exact macrostates.
# So it requires more elaborate set-up, and the transition matrix is large.  Also, the Multistrand simulations are slower than usual.
# The real danger is that (perhaps) the majority of trajectories almost never go through our "strings in state space".
# If these strings are really the center of a "river" between "steep mountainous slopes", then we should be OK; otherwise, maybe not.

def setup_options_hybridization_Exact(trials, toe5_seq, stem_seq, loop_seq, toe3_seq, temp=298.15, conc=1e-3, time=0.0002):

    # I am not sure whether 0-length domains will break Multistrand.  Don't take chances for now...
    assert len(toe5_seq)>0
    assert len(stem_seq)>0
    assert len(loop_seq)>0
    assert len(toe3_seq)>0
    # Thus, to approximate a zero-toehold case, use 1 nt for both toeholds, with complementary sequences... 

    L = len(toe5_seq+stem_seq+loop_seq+stem_seq+toe3_seq)
    
    toe5 = Domain(name="toe5",sequence=toe5_seq,length=len(toe5_seq))
    stem = Domain(name="stem",sequence=stem_seq,length=len(stem_seq))
    loop = Domain(name="loop",sequence=loop_seq,length=len(loop_seq))
    toe3 = Domain(name="toe3",sequence=toe3_seq,length=len(toe3_seq))
    
    sense_strand     = toe5+stem+loop+stem.C+toe3
    antisense_strand = toe3.C+stem+loop.C+stem.C+toe5.C

    # start with individual strands folded as hairpins, end with complete hybrid
    sense_folded     = Complex(strands=[sense_strand], structure=".(.).")
    antisense_folded = Complex(strands=[antisense_strand], structure=".(.).")
    hybrid_complete  = Complex(strands=[sense_strand,antisense_strand], structure="(((((+)))))")

    print("Macrostates: exact configurations along idealized pathways for hybridization, three-way branch migration, and four-way branch migration.")

    macros = []

    # single-strands opening and re-forming hairpins
    for i in range(len(stem_seq)):
        partial_sense=Complex(strands=[sense_strand], structure = "."*len(toe5_seq)+"."*(len(stem_seq)-i)+"("*i+"."*len(loop_seq)+")"*i+"."*(len(stem_seq)-i)+"."*len(toe3_seq))
        macros.append( Macrostate("S"+str(i), [(partial_sense,Exact_Macrostate,0)]) )
        partial_antisense=Complex(strands=[sense_strand], structure = "."*len(toe3_seq)+"."*(len(stem_seq)-i)+"("*i+"."*len(loop_seq)+")"*i+"."*(len(stem_seq)-i)+"."*len(toe5_seq))
        macros.append( Macrostate("A"+str(i), [(partial_antisense,Exact_Macrostate,0)]) )

    # hybridization / fraying:  formation from the inside out, fraying from the outside in
    for i in range(int(L/2)):
        partial_hybrid=Complex(strands=[sense_strand,antisense_strand], structure = "."*i + "("*(L-2*i) + "."*i + "+" + "."*i + ")"*(L-2*i) + "."*i)
        macros.append( Macrostate("H"+str(L-2*i), [(partial_hybrid,Exact_Macrostate,0)]) )

    # four-way branch migration: first form one toehold, then form another toehold, then branch migrate, then close the bubble
    for i in range(1,len(toe5_seq)+1):  # toe5 first...
        s1= "."*(len(toe5_seq)-i)+"("*i+"("*len(stem_seq)+"."*len(loop_seq)+")"*len(stem_seq)+"."*len(toe3_seq)
        s2= "."*len(toe3_seq)+"("*len(stem_seq)+"."*len(loop_seq)+")"*len(stem_seq)+")"*i+"."*(len(toe5_seq)-i)
        partial_toes=Complex(strands=[sense_strand,antisense_strand], structure = s1+"+"+s2)
        macros.append( Macrostate("F5-"+str(i), [(partial_toes,Exact_Macrostate,0)]) )
    for i in range(1,len(toe3_seq)+1):  # then toe3...
        s1= "("*len(toe5_seq)+"("*len(stem_seq)+"."*len(loop_seq)+")"*len(stem_seq)+"("*i+"."*(len(toe3_seq)-i)
        s2= "."*(len(toe3_seq)-i)+")"*i+"("*len(stem_seq)+"."*len(loop_seq)+")"*len(stem_seq)+")"*len(toe5_seq)
        partial_toes=Complex(strands=[sense_strand,antisense_strand], structure = s1+"+"+s2)
        macros.append( Macrostate("F53-"+str(i), [(partial_toes,Exact_Macrostate,0)]) )
    for i in range(1,len(toe3_seq)+1):  # toe3 first...
        s1= "."*len(toe5_seq)+"("*len(stem_seq)+"."*len(loop_seq)+")"*len(stem_seq)+"("*i+"."*(len(toe3_seq)-i)
        s2= "."*(len(toe3_seq)-i)+")"*i+"("*len(stem_seq)+"."*len(loop_seq)+")"*len(stem_seq)+"."*len(toe5_seq)
        partial_toes=Complex(strands=[sense_strand,antisense_strand], structure = s1+"+"+s2)
        macros.append( Macrostate("F3-"+str(i), [(partial_toes,Exact_Macrostate,0)]) )
    for i in range(1,len(toe5_seq)+1):  # then toe5...
        s1= "."*(len(toe5_seq)-i)+"("*i+"("*len(stem_seq)+"."*len(loop_seq)+")"*len(stem_seq)+"("*len(toe3_seq)
        s2= ")"*len(toe3_seq)+"("*len(stem_seq)+"."*len(loop_seq)+")"*len(stem_seq)+")"*i+"."*(len(toe5_seq)-i)
        partial_toes=Complex(strands=[sense_strand,antisense_strand], structure = s1+"+"+s2)
        macros.append( Macrostate("F35-"+str(i), [(partial_toes,Exact_Macrostate,0)]) )
    for i in range(1,len(stem_seq)+1): # now branch migrate (closed intermediates only)
        s1= "("*len(toe5_seq)+"("*len(stem_seq)+"."*len(loop_seq)+")"*(len(stem_seq)-i)+"("*i+"("*len(toe3_seq)
        s2= ")"*len(toe3_seq)+")"*i+"("*(len(stem_seq)-i)+"."*len(loop_seq)+")"*len(stem_seq)+")"*len(toe5_seq)
        partial_4way=Complex(strands=[sense_strand,antisense_strand], structure = s1+"+"+s2)
        macros.append( Macrostate("Fc"+str(i), [(partial_4way,Exact_Macrostate,0)]) )
    for i in range(1,int(len(loop_seq)/2)): # now close the bubble
        s1= "("*len(toe5_seq)+"("*len(stem_seq)+"("*i+"."*(len(loop_seq)-2*i)+"("*i+"("*len(stem_seq)+"("*len(toe3_seq)
        s2= ")"*len(toe3_seq)+")"*len(stem_seq)+")"*i+"."*(len(loop_seq)-2*i)+")"*i+")"*len(stem_seq)+")"*len(toe5_seq)
        partial_bubble=Complex(strands=[sense_strand,antisense_strand], structure = s1+"+"+s2)
        macros.append( Macrostate("Fb"+str(i), [(partial_bubble,Exact_Macrostate,0)]) )

    # three-way branch migration by open sense strand: first form 5' toehold, then branch migrate, then zipper up
    for i in range(1,len(toe5_seq)+1):  # toe5 first...
        s1= "."*(len(toe5_seq)-i)+"("*i+"."*len(stem_seq)+"."*len(loop_seq)+"."*len(stem_seq)+"."*len(toe3_seq)
        s2= "."*len(toe3_seq)+"("*len(stem_seq)+"."*len(loop_seq)+")"*len(stem_seq)+")"*i+"."*(len(toe5_seq)-i)
        partial_toes=Complex(strands=[sense_strand,antisense_strand], structure = s1+"+"+s2)
        macros.append( Macrostate("TS5t-"+str(i), [(partial_toes,Exact_Macrostate,0)]) )
    for i in range(1,len(stem_seq)+1):  # now branch migrate (closed intermediates only)
        s1= "("*len(toe5_seq)+"("*i+"."*(len(stem_seq)-i)+"."*len(loop_seq)+"."*len(stem_seq)+"."*len(toe3_seq)
        s2= "."*len(toe3_seq)+"."*i+"("*(len(stem_seq)-i)+"."*len(loop_seq)+")"*len(stem_seq)+")"*len(toe5_seq)
        partial_toes=Complex(strands=[sense_strand,antisense_strand], structure = s1+"+"+s2)
        macros.append( Macrostate("TS5o-"+str(i), [(partial_toes,Exact_Macrostate,0)]) )
    for i in range(1,len(loop_seq+stem_seq+toe3_seq)): # now zipper up
        s1= "("*len(toe5_seq)+"("*len(stem_seq)+"("*i+"."*(len(loop_seq+stem_seq+toe3_seq)-i)
        s2= "."*(len(loop_seq+stem_seq+toe3_seq)-i)+")"*i+")"*len(stem_seq)+")"*len(toe5_seq)
        partial_toes=Complex(strands=[sense_strand,antisense_strand], structure = s1+"+"+s2)
        macros.append( Macrostate("T5z-"+str(i), [(partial_toes,Exact_Macrostate,0)]) )

    # three-way branch migration by open antisense strand: first form 5' toehold, then branch migrate, then zipper up
    for i in range(1,len(toe5_seq)+1):  # toe5 first...
        s1= "."*(len(toe5_seq)-i)+"("*i+"("*len(stem_seq)+"."*len(loop_seq)+")"*len(stem_seq)+"."*len(toe3_seq)
        s2= "."*len(toe3_seq)+"."*len(stem_seq)+"."*len(loop_seq)+"."*len(stem_seq)+")"*i+"."*(len(toe5_seq)-i)
        partial_toes=Complex(strands=[sense_strand,antisense_strand], structure = s1+"+"+s2)
        macros.append( Macrostate("TA5t-"+str(i), [(partial_toes,Exact_Macrostate,0)]) )
    for i in range(1,len(stem_seq)+1):  # now branch migrate (closed intermediates only)
        s1= "("*len(toe5_seq)+"("*i+"("*(len(stem_seq)-i)+"."*len(loop_seq)+")"*len(stem_seq)+"."*len(toe3_seq)
        s2= "."*len(toe3_seq)+"."*len(stem_seq)+"."*len(loop_seq)+"."*(len(stem_seq)-i)+")"*i+")"*len(toe5_seq)
        partial_toes=Complex(strands=[sense_strand,antisense_strand], structure = s1+"+"+s2)
        macros.append( Macrostate("TA5o-"+str(i), [(partial_toes,Exact_Macrostate,0)]) )
        # zippering up is the same as for the open sense strand pathway

    # three-way branch migration by open sense strand: first form 3' toehold, then branch migrate, then zipper up
    for i in range(1,len(toe3_seq)+1):  # toe3 first...
        s1= "."*len(toe5_seq)+"."*len(stem_seq)+"."*len(loop_seq)+"."*len(stem_seq)+"("*i+"."*(len(toe3_seq)-i)
        s2= "."*(len(toe3_seq)-i)+")"*i+"("*len(stem_seq)+"."*len(loop_seq)+")"*len(stem_seq)+"."*len(toe5_seq)
        partial_toes=Complex(strands=[sense_strand,antisense_strand], structure = s1+"+"+s2)
        macros.append( Macrostate("TS3t-"+str(i), [(partial_toes,Exact_Macrostate,0)]) )
    for i in range(1,len(stem_seq)+1):  # now branch migrate (closed intermediates only)
        s1= "."*len(toe5_seq)+"."*len(stem_seq)+"."*len(loop_seq)+"."*(len(stem_seq)-i)+"("*i+"("*len(toe3_seq)
        s2= ")"*len(toe3_seq)+")"*i+"("*(len(stem_seq)-i)+"."*len(loop_seq)+")"*(len(stem_seq)-i)+"."*i+"."*len(toe5_seq)
        partial_toes=Complex(strands=[sense_strand,antisense_strand], structure = s1+"+"+s2)
        macros.append( Macrostate("TS3o-"+str(i), [(partial_toes,Exact_Macrostate,0)]) )
    for i in range(1,len(loop_seq+stem_seq+toe5_seq)): # now zipper up
        s1= "."*(len(toe5_seq+stem_seq+loop_seq)-i)+"("*i+"("*len(stem_seq)+"("*len(toe3_seq)
        s2= ")"*len(toe3_seq)+")"*len(stem_seq)+")"*i+"."*(len(toe5_seq+stem_seq+loop_seq)-i)
        partial_toes=Complex(strands=[sense_strand,antisense_strand], structure = s1+"+"+s2)
        macros.append( Macrostate("T3z-"+str(i), [(partial_toes,Exact_Macrostate,0)]) )

    # three-way branch migration by open antisense strand: first form 3' toehold, then branch migrate, then zipper up
    for i in range(1,len(toe5_seq)+1):  # toe3 first...
        s1= "."*len(toe5_seq)+"("*len(stem_seq)+"."*len(loop_seq)+")"*len(stem_seq)+"("*i+"."*(len(toe3_seq)-i)
        s2= "."*(len(toe3_seq)-i)+")"*i+"."*len(stem_seq)+"."*len(loop_seq)+"."*len(stem_seq)+"."*len(toe5_seq)
        partial_toes=Complex(strands=[sense_strand,antisense_strand], structure = s1+"+"+s2)
        macros.append( Macrostate("TA3t-"+str(i), [(partial_toes,Exact_Macrostate,0)]) )
    for i in range(1,len(stem_seq)+1):  # now branch migrate (closed intermediates only)
        s1= "."*len(toe5_seq)+"("*len(stem_seq)+"."*len(loop_seq)+")"*(len(stem_seq)-i)+"("*i+"("*len(toe3_seq)
        s2= ")"*len(toe3_seq)+")"*i+"."*(len(stem_seq)-i)+"."*len(loop_seq)+"."*(len(stem_seq)-i)+"."*i+"."*len(toe5_seq)
        partial_toes=Complex(strands=[sense_strand,antisense_strand], structure = s1+"+"+s2)
        macros.append( Macrostate("TA3o-"+str(i), [(partial_toes,Exact_Macrostate,0)]) )
        # zippering up is the same as for the open sense strand pathway

    initial_hairpins      = Macrostate("HAIRPINS", [(sense_folded,Exact_Macrostate,0),(antisense_folded,Exact_Macrostate,0)])   # perfect hairpins achieved
    completed_hybrid      = Macrostate("COMPLETE", [(hybrid_complete,Exact_Macrostate,0)])                                      # perfect hybrid achieved

    o = Options(simulation_mode="Transition", parameter_type="Nupack", dangles="Some",
                rate_method = "Metropolis", num_simulations = trials, simulation_time=time,  # default 0.0002 sec
                substrate_type="DNA", temperature=temp, join_concentration=conc,             # default 1 mM at 25C
                start_state=[sense_folded,antisense_folded], rate_scaling='Calibrated', verbosity=0)

    # for transition mode, these are mostly macrostates not stop conditions, since the simulation won't stop unless the name begins with "stop"
    o.stop_conditions = [initial_hairpins]+macros+[completed_hybrid]
    return o


# mol will be a list of True/False for which transition macrostates the system has entered
# so in_state(mol) returns True if the system is in at least one of the listed macrostates.
def in_state( mol ): return sum(mol) > 0


# Gives shorthand name number n.
# The shorthand is a single character if no more than 52 macrostates are defined; goes into multiple characters otherwise.
def shorthand(n): 
    charindex = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
    m = len(charindex)
    nm = charindex[n%m]
    while n>=m:
        n=n/m  # integer arithmetic, truncates
        nm=charindex[n%m]+nm
    return nm


# mol is a Boolean descriptor of macrostate occupancy, like mol above.
# A shorthand name for this macrostate (based on the order given in stop_conditions) is provided.
def mol_name(mol):
    names = [shorthand(j) for i,j in zip(mol,range(len(mol))) if i]
    if names == []:
        names = '0'  
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
# Reminder: mol = macrostate = subset of defined macrostates that we are currently in.  
# Trivial mols are the empty subset, i.e. the macrostate of not being in a defined macrostate.
# An error is generated if either the initial or final macrostate is not visited.
def parse_transition_lists( transition_traj_list ):
    transition_dict = {}
    first_name = None   # find the full name, indicating all defined macrostates that we are in, if we visit the initial exact state
    last_name = None    # find the full name, indicating all defined macrostates that we are in, if we visit the final exact state
    hitting_times = []

    # the mol1 --> mol2 transition times represent (time of entry into mol1) to (time of entry into mol2)
    for transition_traj in transition_traj_list:
        truncated = [i for i in transition_traj if in_state(i[1])]
        tt = truncated # only keep the entry times and mol states for non-trivial mols
        hit_yet=False

        if tt[0][1][0]:     # mol says we're in the first defined macrostate, i.e. initial state
            first_name = mol_name(tt[0][1])
        if tt[0][1][-1]:    # mol says we're in the last defined macrostate, i.e. final state (unexpected, but possible)
            last_name = mol_name(tt[0][1])
            hitting_times.append(tt[0][0])
            hit_yet=True

        for i in range(len(tt)-1):
            nm = trans_name(tt[i][1],tt[i+1][1])
            if nm in transition_dict:
                transition_dict[nm].append( tt[i+1][0] - tt[i][0] )
            else:
                transition_dict[nm] = [tt[i+1][0] - tt[i][0]]
            if tt[i+1][1][0]:   # mol says we're in the first defined macrostate, i.e. initial state
                first_name = mol_name(tt[i+1][1])
            if tt[i+1][1][-1]:  # mol says we're in the last defined macrostate, i.e. final state
                last_name = mol_name(tt[i+1][1])
                if not hit_yet:
                    hitting_times.append(tt[i+1][0])
                    hit_yet=True

    if first_name == None or last_name == None:
        print("Error: first state %s and/or last state %s were not visited" % (first_name, last_name))
        import pdb; pdb.set_trace()

    return (transition_dict, first_name, last_name, hitting_times)


def parse_transition_list( transition_traj_list ):
    return parse_transition_lists( [transition_traj_list] )


# For a second-order chain, each state is now the pair (previous mol, current mol), named e.g. "A,B;X,Y,Z".  
# This is useful if transitions between macrostates are not Markov; tracking the history of most recent macrostate may make the chain "more Markov".
# For this purpose, we want to skip X->X transitions, since the history is a surrogate for "from what direction we entered the macrostate".
# Because they are (must be!) exact microstates, the history doesn't matter for the initial state and for the final state, and these are left as-is
# (which is essential for them to be identified uniquely for the forward and reverse path computations).
def parse2_transition_lists( transition_traj_list ):
    transition_dict = {}
    first_name = None   # find the full name, indicating all defined macrostates that we are in, if we visit the initial exact state
    last_name = None    # find the full name, indicating all defined macrostates that we are in, if we visit the final exact state
    hitting_times = []

    # the mol1 --> mol2 transition times represent (time of entry into mol1) to (time of entry into mol2)
    for transition_traj in transition_traj_list:

        truncated = [i for i in transition_traj if in_state(i[1])] # only keep the entry times and mol states for non-trivial mols
        tt = truncated 
        hit_yet=False

        prev    = tt[0]                # this will track the previous mol state 
        prev_name = mol_name(prev[1])  # name soon becomes second-order name-pair
        if not prev[1][0]:             # prev[1] is the mol state, i.e. a boolean vector of which macrostates we are in.  
            print("Expected to start in state containing 'A', but actually started in %s." % (prev_name))
            import pdb; pdb.set_trace()

        for i in range(len(tt)-1):
            if not tt[i][1] == tt[i+1][1]:   # skip X->X transitions (i.e. excursions around the edges of X
                curr = tt[i+1]
                curr_name = mol_name(curr[1])
                if curr[1][0]:    # mol says we're in the first defined macrostate, i.e. initial state
                    first_name = curr_name
                elif curr[1][-1]: # mol says we're in the last defined macrostate, i.e. final state
                    last_name = curr_name
                    if not hit_yet:
                        hitting_times.append(curr[0])
                        hit_yet=True
                else:             # not in either initial or final state
                    curr_name = mol_name(prev[1]) + ";" + mol_name(curr[1])
                nm = prev_name + " -> " + curr_name
                if nm in transition_dict:
                    transition_dict[nm].append( curr[0] - prev[0] )
                else:
                    transition_dict[nm] = [curr[0] - prev[0]]
                prev_name = curr_name
                prev = curr

    if first_name == None or last_name == None:
        print("Error: first state %s and/or last state %s were not visited" % (first_name, last_name))
        import pdb; pdb.set_trace()
        
    return (transition_dict, first_name, last_name, hitting_times)


def parse2_transition_list( transition_traj_list ):
    return parse2_transition_lists( [transition_traj_list] )


# print the true names of the defined macrostates
def print_true_names( options, first, last ):
    print("True names of all %d defined macrostates and their abbreviations:" % (len(options.stop_conditions)))
    for i,idx in zip(options.stop_conditions,range(len(options.stop_conditions))):
        print(("{0}: {1}".format( i.tag, shorthand(idx))))
    print("Initial state (intersection of defined macrostates) is: %s" % first)
    print("Final state (intersection of defined macrostates) is: %s" % last)


def print_transition_dict( transition_dict ):
    k = list(transition_dict.keys())
    k.sort() 

    visited_states = set()
    for i in k:
        visited_states.add( i.split(" -> ")[0] )
        visited_states.add( i.split(" -> ")[1] )

    print("  Transitions averaged over all %d simulations (%d coarse-grained states visited):" % (o.num_simulations,len(visited_states)))
    for i in k:
        transition_times = np.array( transition_dict[i] )
        print(("{0}: {2:.2e} ({1})".format(i,len(transition_dict[i]),np.mean(transition_times))))
    

# The first defined macrostate must be the initial macrostate, and must be exact;
# the last defined macrostate must be the final macrostate, and must be exact.
# The returned list will begin with "first" (presumably the mol containing the initial macrostate)
# and end with "last" (presumably the mol containing the final macrostate).
# The order of the other visited macrostates (mols) is not defined.
# Recall that a mol is a boolean indicator of membership in the defined macrostates, specifying an intersection macrostate.
def visited_macrostates( transition_dict, first, last):

    # each key is a transition between multisets of macrostates, e.g. "A,B -> C,e,f"
    visited_macs = []
    for trans in list(transition_dict.keys()):
        visited_macs += trans.split(" -> ")
    visited_macs = list(set(visited_macs))

    # the initial state must be the first one, i.e. index 0.
    # the final stop state must be the last one, i.e. index n-1;  

    if (not first in visited_macs) or (not last in visited_macs):
        print("Error: initial or final macrostate not found in visited macrostate list.")
        import pdb; pdb.set_trace()

    visited_macs = [ first ] + list(set(visited_macs).difference(set([first, last]))) + [ last ]
    return visited_macs


def build_transition_matrices( transition_dict, visited_macs ):
    
    n   = len(visited_macs)  # set of states for chain --  note this isn't quite treated as a markov chain, since P_ij =!= (1/dt_ij) / sum_j (1/dt_ij)
    dt  = np.zeros((n,n))    # average time for transition from state i to j
    P   = np.zeros((n,n))    # probability of transitioning from state i to j, given in state i
    N   = np.zeros((n,n))    # number of transitions from state i to j  --  P_ij will be computed as N_ij / sum_j (N_ij)
    out = np.zeros(n)        # number of transitions out of state i -- computed as sum_j (N_ij)

    for trans in list(transition_dict.keys()):
        (mac_from, mac_to) = tuple( trans.split(" -> ") )
        if mac_from in visited_macs and mac_to in visited_macs:   # always true, unless subset of all macs are used
            i = visited_macs.index(mac_from)
            j = visited_macs.index(mac_to)
            transition_times = np.array( transition_dict[trans] )
            dt[i,j] = np.mean(transition_times)
            N[i,j]  = len(transition_dict[trans])
            out[i] += len(transition_dict[trans])
    
    # calculate P from N
    for i in range(n):
        for j in range(n):
            if out[i] == 0:
                P[i,j] = 0.0
            else:
                P[i,j] = float(N[i,j]) / float(out[i])

    return dt, P, N

def compute_hitting_times( transition_dict, visited_macs ):


    # visited macs must be sorted such that the first state is first, and the final state is last
    (dt, P, N) = build_transition_matrices( transition_dict, visited_macs )

    expected_dt = np.sum( dt * P, axis = 1)  #  == sum_j ( dt_ij * P_ij )
    # print "Expected time for transitions out of state %d = %g sec." % (len(expected_dt),expected_dt[-1])   # should be zero

    expected_dt = expected_dt[0:-1]  # omit the expected time out of the final state
    P_hat = P[0:-1,0:-1]             # omit transitions into and out of the final state

    n = len(expected_dt)             # Here, n refers to len(visited_macs)-1, because "last" is excluded (unlike in build_transition_matrices, above).
    id = np.identity(n)

    hitting = np.linalg.solve(id-P_hat,expected_dt)   #  (id-P_hat) x hitting = expected_dt

    if np.allclose( np.dot(id-P_hat,hitting), expected_dt ) and np.all( hitting >= 0 ):
        print("Solved linear equations successfully!")
    else:
        print("Linear equations not solved well.")
        if np.all( hitting >= 0 ):
            print("Solution is not close to posed equality.")
        else:
            print("Some hitting times are negative.")

    # convert hitting back into a regular list, and add 0.0 for the expected time from last to last
    hitting = list(hitting)+[0.0]

    # we sort both lists according to hitting time!  
    for (s,h) in sorted(zip(visited_macs,hitting+[0.0]), key=lambda pair: -pair[1]):
        print("Mean time from %15s to %15s is %g sec." % (s,visited_macs[-1],h))
    return hitting


def tally_residence_times( transition_traj_list, visited_macs):

    n = len(visited_macs)      # just the mol for 1st order chains;  previous+current mol for 2nd order chains (separated by ";")
    first = visited_macs[0]
    last = visited_macs[-1]
    residence_times = [0]*n
    trivial_time = 0

    second_order = ( n > 2 and ";" in visited_macs[1])   # in this case, we need to tally residence in second-order states.

    # the mol1 --> mol2 transition times represent (time of entry into mol1) to (time of entry into mol2)
    for transition_traj in transition_traj_list:
        tt = transition_traj  # entry times to any mol, either trivial or non-trivial

        for i in range(len(tt)-1):
            if in_state(tt[i][1]):  # tally times in (intersections of) defined macrostates.
                if not second_order:
                    w = visited_macs.index(mol_name(tt[i][1]))
                    residence_times[w] += (tt[i+1][0]-tt[i][0]) # time entering something else - time entering this mol
                else:
                    curr = mol_name(tt[i][1])
                    if curr == first or curr == last:
                        curr2 = curr
                    elif not prev == curr:
                        curr2 = prev + ";" + curr
                    else:
                        pass     # if the previous mol has not changed, then variable curr2 remains correct for the second-order chain
                    w = visited_macs.index(curr2)
                    residence_times[w] += (tt[i+1][0]-tt[i][0]) 
                    prev = curr
            else:
                trivial_time += (tt[i+1][0]-tt[i][0]) # time entering something else - time entering trivial mol (i.e. unlabeled state)

    return (residence_times, trivial_time)


def compute_most_likely_path( transition_dict, visited_macs, options = None ):
    from scipy.sparse.csgraph import shortest_path

    n = len(visited_macs)
    (dt, P, N) = build_transition_matrices( transition_dict, visited_macs )
    np.putmask(P,np.identity(n)>0,np.zeros([n,n])) # ignore self-transitions
    newP = []
    for i in range(n):                             # re-normalize out of each state
        rowsum = sum(P[i,:])
        if rowsum == 0.0:
            newP.append( list(P[i,:]) )
        else:
            newP.append( [v/rowsum for v in P[i,:]] )
    P = np.array(newP)

    # we want to find the argmax_path of prod_{(i,j) in path} prob(i to j)
    #   = argmin_path prod_{(i,j) in path} 1/prob(i to j)
    #   = argmin_path sum_{(i,j) in path} log(1/prob(i to j))
    old_err_state = np.seterr(divide='ignore')
    distance_graph = np.log(1/P)  # all non-negative numbers, with 'inf' for zero-probability transitions, and 0 for probability-one transitions.
    m_distance_graph = np.ma.masked_array( distance_graph, mask = (distance_graph==np.inf) )
    np.seterr(**old_err_state)    # explicitly masked arrays are necessary because probability=1 transitions have log=0, a value which shortest_path treats as a non-edge by default.

    path_matrix = shortest_path(m_distance_graph)   # value[i,j] == min distance from i to j

    i = 0
    path = [i]
    pathstring = visited_macs[i]
    probstring = ""
    while not i == n-1:
        neighbors = [j for j in range(n) if not distance_graph[i,j] == np.inf]
        j = neighbors[ np.argmin(path_matrix[neighbors,n-1]) ]
        # print "at %d, neighbors %s, distances %s, going to %d" % (i, str(neighbors), str(path_matrix[neighbors,n-1]), j)
        pathstring += " -> " + visited_macs[j]
        probstring += ("[%d/%d]" % (N[i,j],int(round(N[i,j]/P[i,j]))))
        if j in path:
            print("Most likely path has a loop?!?")
            import pdb; pdb.set_trace()
        i=j
        path += [i]

    print("The most likely path through macrostates:")
    print(pathstring)

    # also print the true names of the macrostates, if an Options object is provided
    if options:
        namedict = dict()  # reverse-lookup-table
        for m,idx in zip(options.stop_conditions,range(len(options.stop_conditions))):
            namedict[shorthand(idx)] = m.tag
        def lookup(s):  # is a shorthand representation of a subsets of defined macrosets (i.e. first-order shorthand name)
            parts = [namedict[n] for n in s.split(',')]
            return ','.join(parts)
        verbosepathstring = ''
        path_macs = pathstring.split(' -> ')
        for mac in path_macs:
            if not verbosepathstring == '':
                verbosepathstring += ' -> '
            verbosepathstring += '{' + ';'.join([ lookup(s) for s in mac.split(';') ]) + '}' # for second-order states with history
        print("a.k.a.")
        print(verbosepathstring)

    print("with probability %s" % probstring)

    back_to_start = path_matrix[:,0]<np.inf  # which states can reach back to start
    return (path,back_to_start)  # return the path, along with information about what the reverse paths might be found


if __name__ == '__main__':

    if not len(sys.argv) == 4:
        print("Error: First argument must be 'small' or 'medium' or 'large', specifying which molecule to simulate.") 
        print("Error: Second argument must be 'loose' or 'count' or 'distance' or 'exact' (or 'loose2' or 'count2' or 'distance2' or 'exact2') for second-order coarse-graining).")
        print("Error: Third argument must be a positive integer, the number of trajectories to be simulated.")
        sys.exit()
    if not sys.argv[1] in ['small', 'medium', 'large']:
        print("Error: First argument must be 'small' or 'medium' or 'large', specifying which molecule to simulate.") 
        sys.exit()
    if not sys.argv[2] in ['loose', 'count', 'distance', 'loose2', 'count2', 'distance2', 'exact', 'exact2']:
        print("Error: Second argument must be 'loose' or 'count' or 'distance' or 'exact' (or 'loose2' or 'count2' or 'distance2' or 'exact2') for second-order coarse-graining).")
        sys.exit()
    try:
        num_traj = int(sys.argv[3])
    except TypeError:
        print("Error: Third argument must be a positive integer, the number of trajectories to be simulated.")
        sys.exit()
    if num_traj<1:
        print("Error: Third argument must be a positive integer, the number of trajectories to be simulated.")
        sys.exit()
    molecule = str(sys.argv[1])
    macrotype = str(sys.argv[2])

    if molecule == "small":
        (toe1,stem,loop,toe2) = ("TC","TA","CTTCTAC","TAC")
    if molecule == "medium":
        (toe1,stem,loop,toe2) = ("ACTCCA","CATTGTC","TTTATATTGT","AAGCGT")
    if molecule == "large":
        (toe1,stem,loop,toe2) = ("CATTGTAACT","GGCGATGCTACCTGTA","TTTT","CCATTAACTC")

    # loose is faster than count;  takes maybe an hour for 500 simulations of 'medium'
    secondorder = (macrotype[-1]=='2')
    if macrotype in ['loose', 'loose2']:
        o = setup_options_hybridization_Loose(num_traj,toe5_seq=toe1, stem_seq=stem, loop_seq=loop, toe3_seq=toe2)
    elif macrotype in ['count', 'count2']:
        o = setup_options_hybridization_Count(num_traj,toe5_seq=toe1, stem_seq=stem, loop_seq=loop, toe3_seq=toe2)
    elif macrotype in ['distance', 'distance2']:
        o = setup_options_hybridization_Distance(num_traj,toe5_seq=toe1, stem_seq=stem, loop_seq=loop, toe3_seq=toe2)
    elif macrotype in ['exact', 'exact2']:
        o = setup_options_hybridization_Exact(num_traj,toe5_seq=toe1, stem_seq=stem, loop_seq=loop, toe3_seq=toe2)
    else:
        print("Error: Macrostate type not defined.")
        sys.exit()

    print("--- Running %d simulations for '%s' and its complement ---" % (num_traj,o.start_state[0].sequence))

    s=SimSystem(o)
    s.start()

    print()
    print("--- Analysis of simulations with transition macrostates ---")

    # print "  Coarse-grained trajectory of simulation #1:"
    # print_transitions(o.interface.transition_lists[0])
    # print "  Transitions from simulation #1:"
    # parsedlist = parse_transition_list(o.interface.transition_lists[0])
    # print_transition_dict(parsedlist)
    if secondorder:
        (parsedlist,first_state,last_state,hitting_times) = parse2_transition_lists(o.interface.transition_lists)
    else:
        (parsedlist,first_state,last_state,hitting_times) = parse_transition_lists(o.interface.transition_lists)
    print_true_names( o, first_state, last_state )
    print("Average hitting time is %g for %d out of %d trials that reached the final state." % (np.mean(hitting_times),len(hitting_times), len(o.interface.transition_lists)))

    macs = visited_macrostates( parsedlist, first_state, last_state )
    (residence_times, trivial_time) = tally_residence_times( o.interface.transition_lists, macs )
    total_time = np.sum(residence_times)+trivial_time
    print("Spent %g sec out of %g sec (i.e. %4.1f%%) in unlabeled states." % ( trivial_time, total_time, 100.0 * trivial_time / total_time ))
    print_transition_dict(parsedlist) 

    times = compute_hitting_times( parsedlist, macs )
    (path,back_to_start) = compute_most_likely_path( parsedlist, macs, o )

    # we're interested in going from (almost) the final state, to the initial state
    print()
    print("Now examining the reverse process, from duplex to single strands...")

    if not np.argmin(times) == len(times)-1:
        print("  The hitting time from the last state is not the fastest... huh?")

    if back_to_start[-1]:  # this fails if the final state cannot go back and reach the start state.
        rev_macs = [ last_state ] + list(set(macs).difference(set([first_state, last_state]))) + [ first_state ]
        rev_times = compute_hitting_times( parsedlist, rev_macs )    
        (rev_path,to_the_end) = compute_most_likely_path( parsedlist, rev_macs, o )

        # compare to expectation from formula k_r = k_f * exp(dG/RT)
        print("Compare to expectation:")
        print("  Forward hitting time %g corresponds to kf = %g" % (times[0], times[1]))

    else:
        print("  The final macrostate on the forward path has no observed transition pathways back to the initial state.")
        print("  Thus our current reverse-path algorithm is not safe to run.")
