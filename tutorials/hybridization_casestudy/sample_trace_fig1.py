# FD, May 17th, 2017. 
# This demonstrates the pair-type functionality
# For a given complex, Pairtype returns a unique representation. 
# This is important for hassing functions 

import sys
import numpy as np

from multistrand.objects import Complex, Domain, Strand, StopCondition
from multistrand.options import Options
from multistrand.system import SimSystem
from multistrand.utils import pairType
from multistrand.experiment import standardOptions, makeComplex

    

def printTrajectory(o):
    
    seqstring = ""
    
    for i in range(len(o.full_trajectory)):
    
        time = 1e9 * o.full_trajectory_times[i]
        states = o.full_trajectory[i]
        
        ids = []
        newseqs = []
        structs = []
        dG = 0.0;
        dGC = 0.0
        
        pairTypes = []
        
        for state in states: 
            
            ids += [ str(state[2]) ]
            newseqs += [ state[3] ]  # extract the strand sequences in each complex (joined by "+" for multistranded complexes)
            structs += [ state[4] ]  # similarly extract the secondary structures for each complex
            dG += dG + state[5]
            
            dGC += (state[5] - (  o._temperature_kelvin * 0.0019872036 * np.log(  1.0 / o.join_concentration) * state[4].count("+") )) 
            
#             print "count is  " +  str(state[4].count("+"))
#             print "join conc is " + str(o.join_concentration)
#             print "dG-Complex is " + "%.2f" % dGC + " kcal/mol  for " + str(state[3]) 
            
            
        newseqstring = ' '.join(newseqs)  # make a space-separated string of complexes, to represent the whole tube system sequence
        tubestruct = ' '.join(structs)  # give the dot-paren secondary structure for the whole test tube
                 
        
        if not newseqstring == seqstring : 
            print newseqstring
            seqstring = newseqstring  # because strand order can change upon association of dissociation, print it when it changes        

        print tubestruct + ('   t=%.3f ns,  dG=%3.2f kcal/mol, dGC=%3.2f kcal/mol   ' % (time, dG, dGC)) 

        

def doSims(strandSeq, numTraj=2):    

    o1 = standardOptions()
    
    o1.simulation_mode = Options.trajectory
    o1.num_simulations = numTraj
    o1.output_interval = 1
    o1.join_concentration = 1.0e-9
#     o1.join_concentration = 1.0
    o1.simulation_time = 0.000005
    o1.DNA23Metropolis()

    seq1 = "ACTGACTGACTG"
    seq2 = "ACTG"
    seq3 = "CATTCAGTACAGT"
    
    struct = "((((((.((.(.+).((+)).)).)).))))"
    
#     # Using domain representation makes it easier to write secondary structures.
#     domain1 = Domain(name="d1", sequence=seq1)
#     domain2 = Domain(name="d2", sequence=seq2)
#     domain3 = Domain(name="d3", sequence=seq3)
#     
#     strand1 = Strand(name="s1", domains=[domain1])
#     strand2 = Strand(name="s2", domains=[domain2])
#     strand3 = Strand(name="s3", domains=[domain3])

    myComplex = makeComplex([seq1, seq2, seq3], struct)

    print myComplex

    o1.start_state = [myComplex]

    # Note: no stopping conditions
    
    o1.initial_seed = 1777+6

    s = SimSystem(o1)
    s.start()
    printTrajectory(o1)        
        

    
# # The actual main method
if __name__ == '__main__':
    
    print sys.argv
    doSims( "GCGTTTCAC",1)
         
        
        

















