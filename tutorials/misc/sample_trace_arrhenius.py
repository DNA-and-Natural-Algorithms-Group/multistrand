# FD, May 29th, 2018 
# This file demonstrates arrhenius kinetics

import sys

from multistrand.objects import Complex, Domain, Strand, StopCondition
from multistrand.options import Options
from multistrand.system import SimSystem
from multistrand.utils import pairType
from multistrand.experiment import standardOptions, hybridization

    
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
            print newseqstring
            seqstring = newseqstring  # because strand order can change upon association of dissociation, print it when it changes        

        print tubestruct + ('   t=%.6f ms,  dG=%3.2f kcal/mol, arrType %i' % (time, dG, arrType)) 

        

def doSims(strandSeq, numTraj=2):    

    o1 = standardOptions()
    
    o1.num_simulations = numTraj
    o1.output_interval = 1 
    o1.simulation_time = ATIME_OUT
    
    hybridization(o1, strandSeq )
   
    o1.useArr = True

    o1.lnAStack =     1.41839430e+01
    o1.EStack =     5.28692038e+00

    o1.lnAEnd =     1.64236969e+01
    o1.EEnd =     4.46143369e+00

    o1.lnALoop =    1.29648159e+01
    o1.ELoop =    3.49798154e+00

    o1.lnAStackLoop =  5.81061725e+00
    o1.EStackLoop =      -1.12763854e+00

    o1.lnAStackEnd =    1.75235569e+01
    o1.EStackEnd =    2.65589869e+00

    o1.lnALoopEnd =    2.42237267e+00
    o1.ELoopEnd =    8.49339120e-02

    o1.lnAStackStack =    8.04573830e+00
    o1.EStackStack =   -6.27121400e-01

    o1.bimolecular_scaling =    1.60062641e-02

    o1.initial_seed = 1777+6

    s = SimSystem(o1)
    s.start()
    printTrajectory(o1)        
        

    
# # The actual main method
if __name__ == '__main__':
    
    print sys.argv
    doSims( "GCGTTTCAC",1)
         
        
        

















