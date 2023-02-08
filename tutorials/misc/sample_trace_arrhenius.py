"""
    Frits Dannenberg May 29th, 2018 
    This file demonstrates Arrhenius kinetics.
    The interface is not finalized.
"""
from __future__ import print_function
import sys

from multistrand.objects import Complex, Domain, Strand, StopCondition
from multistrand.experiment import standardOptions, hybridization
from multistrand.builder import codeToDesc
from multistrand.options import Options, Literals
from multistrand.system import SimSystem
from multistrand.utils import pairType
    
ATIME_OUT = 0.001


def printTrajectory(o):
    
    seqstring = ""
    
    for i in range(len(o.full_trajectory)):
    
        time = 1e3 * o.full_trajectory_times[i]
        states = o.full_trajectory[i]
        
        ids = []
        newseqs = []
        structs = []
        dG = 0.0;
        
        pairTypes = []

        arrType = o.full_trajectory_arrType[i]   

        for state in states: 
            
            ids += [ str(state[2]) ]
            newseqs += [ state[3] ]  # extract the strand sequences in each complex (joined by "+" for multistranded complexes)
            structs += [ state[4] ]  # similarly extract the secondary structures for each complex
            dG += dG + state[5]

        newseqstring = ' '.join(newseqs)  # make a space-separated string of complexes, to represent the whole tube system sequence
        tubestruct = ' '.join(structs)  # give the dot-paren secondary structure for the whole test tube
        
        if not newseqstring == seqstring : 
            print(newseqstring)
            seqstring = newseqstring  # because strand order can change upon association of dissociation, print it when it changes        

        print(tubestruct + ('   t=%.6f ms,  dG=%3.2f kcal/mol  Type:  %s %s' % (time, dG, codeToDesc(arrType)[0], codeToDesc(arrType)[1])))  
        

def doSims(strandSeq, numTraj=2):    

    o1 = standardOptions()
    
    o1.num_simulations = numTraj
    o1.output_interval = 1 
    o1.simulation_time = ATIME_OUT
    
    hybridization(o1, strandSeq)

    o1.DNA23Arrhenius()

    o1.initial_seed = 1777 + 6

    s = SimSystem(o1)
    s.start()
    printTrajectory(o1)        

    
# # The actual main method
if __name__ == '__main__':
    
    print(sys.argv)
    doSims("GCGTTTCAC", 1)

