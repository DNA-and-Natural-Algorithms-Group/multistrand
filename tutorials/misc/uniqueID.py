# FD, May 17th, 2017. 
# This demonstrates the pair-type functionality
# For a given complex, Pairtype returns a unique representation. 
# This is important for hashing functions 

import sys

from multistrand.objects import Complex, Domain, Strand, StopCondition
from multistrand.options import Options
from multistrand.system import SimSystem
from multistrand.utils import pairType
from multistrand.experiment import standardOptions, hybridization

    
ATIME_OUT = 0.0010

def printTrajectory(o):
    
    seqstring = ""
    
    for i in range(len(o.full_trajectory)):
    
        time = 1e6 * o.full_trajectory_times[i]
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
            
            uniqueID = pairType(state[2], state[4])
            
            pairTypes += [  uniqueID[0] + str(uniqueID[1])   ]
        
                
        identities2 = ' '.join(pairTypes)
        newseqstring = ' '.join(newseqs)  # make a space-separated string of complexes, to represent the whole tube system sequence
        tubestruct = ' '.join(structs)  # give the dot-paren secondary structure for the whole test tube
                 
        
        if not newseqstring == seqstring : 
            print newseqstring
            seqstring = newseqstring  # because strand order can change upon association of dissociation, print it when it changes        

        print tubestruct + (' t=%.4f ms, dG=%3.2f kcal/mol, uID  %s   ' % (time, dG, identities2)) 


def doSims(strandSeq, numTraj=2):    

    o1 = standardOptions()
    
    o1.num_simulations = numTraj
    o1.output_interval = 1 
    
    hybridization(o1, strandSeq)
    o1.initial_seed = 1778

    s = SimSystem(o1)
    s.start()
    printTrajectory(o1)        
        

    

    
# # The actual main method
if __name__ == '__main__':
    
    print sys.argv
    doSims( "GCGTTTCGC",2)
         
        
        

