from __future__ import print_function
# FD, May 17th, 2017. 
# This demonstrates the pair-type functionality
# For a given complex, Pairtype returns a unique representation. 
# This is important for hassing functions 

import sys, time

from multistrand.objects import Complex, Domain, Strand, StopCondition
from multistrand.options import Options
from multistrand.system import SimSystem
from multistrand.utils import pairType
from multistrand.experiment import standardOptions, makeComplex

    
ATIME_OUT = 0.000001

def printTrajectory(o):
    
    seqstring = ""
    print("seed =  %i" % o.initial_seed)
    
    for i in range(len(o.full_trajectory)):
    
        time = 1e3 * o.full_trajectory_times[i]
        states = o.full_trajectory[i]
        
        ids = []
        newseqs = []
        structs = []
        dG = 0.0;
        
        pairTypes = []
        
        for state in states: 
            
            ids += [ str(state[2]) ]
            newseqs += [ state[3] ]  # extract the strand sequences in each complex (joined by "+" for multistranded complexes)
            structs += [ state[4] ]  # similarly extract the secondary structures for each complex
            dG += dG + state[5]
            
        newseqstring = ' '.join(newseqs)  # make a space-separated string of complexes, to represent the whole tube system sequence
        tubestruct = ' '.join(structs)  # give the dot-paren secondary structure for the whole test tube
                 
        
        #not printing anything
        
        if not newseqstring == seqstring : 
            print(newseqstring)
            seqstring = newseqstring  # because strand order can change upon association of dissociation, print it when it changes        

        print(tubestruct + ('   t=%.6f ms,  dG=%3.2f kcal/mol  ' % (time, dG))) 



def doSims(numTraj=2):    

    o1 = standardOptions()
    
    o1.simulation_mode = Options.firstStep
    o1.num_simulations = numTraj
    o1.output_interval = 1 
    o1.simulation_time = ATIME_OUT
    
    myComplex1 = makeComplex(["GTCACTGCTTTT"], "............")
    myComplex2 = makeComplex(["GTCACTGC", "GCAGTGAC"], ".(((((((+))))))).")
        
    o1.start_state = [myComplex1, myComplex2]

#     myComplex = makeComplex(["GTCACTGCTTTT","GTCACTGC","GCAGTGAC"], "..(.........+).((((((+))))))..")
#     o1.start_state = [myComplex]

    # no stop conditions, just timeout.

    o1.initial_seed = time.time() * 1000000
#     o1. initial_seed = 1501710503137097

    s = SimSystem(o1)
    s.start()
    printTrajectory(o1)        
        

    

    
# # The actual main method
if __name__ == '__main__':
    
    print(sys.argv)
    
#     doSims()
    
      
    for i in range(10000):
        time.sleep(0.000001)
        doSims(1)
         
        
        

















