# FD, May 17th, 2017. 
# This demonstrates the pair-type functionality
# For a given complex, Pairtype returns a unique representation. 
# This is important for hassing functions 

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

        print tubestruct + ('   t=%.6f ms,  dG=%3.2f kcal/mol  ' % (time, dG)) 


def doSims(strandSeq, numTraj=2):    

    o1 = standardOptions()
    
    o1.simulation_mode = Options.trajectory
    o1.num_simulations = numTraj
    o1.output_interval = 1 
    o1.simulation_time = ATIME_OUT
    
    
    
    SHORT_SEQ1 = "ACCTCT"
    SHORT_SEQ2 = "TCTTTA"
    SHORT_SEQ7 = "ACATCC"
    SHORT_SEQ5 = "TACTAC"
    SHORT_SEQ6 = "ACCATT"
    SHORT_SEQT = "CTCT"
    
    LONG_SEQ2 = "CCAAACAAAACCTAT"
    LONG_SEQ5 = "AACCACCAAACTTAT"
    LONG_SEQ6 = "CCTAACACAATCACT"
    # some ive made up, but these shouldn't make much difference
    LONG_SEQ7 = "CCACAAAACAAAACT"
    LONG_SEQ1 = "CATCCATTCAACTAT"
    LONG_SEQT = SHORT_SEQT
    
    CL_LONG_S18 = "TCTTCTAACAT"
    CL_LONG_S5 = "CCACCAAACTT"
    CL_LONG_S6 = "TAACACAATCA"
    CL_LONG_S29 = "CCAATACTCCT"
    CL_LONG_S53 = "TATCTAATCTC"
    CL_LONG_S44 = "AAACTCTCTCT"
    CL_LONG_SEQT = "TCT"
    CLAMP_SEQ = "CA"
    
    SHORT_GATE_A_SEQ = [SHORT_SEQ1, SHORT_SEQ2, SHORT_SEQ5, SHORT_SEQ7, SHORT_SEQT]
    SHORT_GATE_B_SEQ = [SHORT_SEQ2, SHORT_SEQ5, SHORT_SEQ6, SHORT_SEQ7, SHORT_SEQT]
    LONG_GATE_A_SEQ = [LONG_SEQ1, LONG_SEQ2, LONG_SEQ5, LONG_SEQ7, LONG_SEQT]
    LONG_GATE_B_SEQ = [LONG_SEQ2, LONG_SEQ5, LONG_SEQ6, LONG_SEQ7, LONG_SEQT]
    CL_LONG_GATE_A_SEQ = [CL_LONG_S44, CL_LONG_S18,
                          CL_LONG_S5, CL_LONG_S29, CL_LONG_SEQT, CLAMP_SEQ]
    CL_LONG_GATE_B_SEQ = [CL_LONG_S53, CL_LONG_S5,
                          CL_LONG_S6, CL_LONG_S29, CL_LONG_SEQT, CLAMP_SEQ]
    
    
    clamped_gate_output_leak(o1, CL_LONG_GATE_A_SEQ, 10)
    
    
    
    
    o1.initial_seed = 1777

    s = SimSystem(o1)
    s.start()
    printTrajectory(o1)        
        

    

    
# # The actual main method
if __name__ == '__main__':
    
    print sys.argv
    doSims( "GCGTTTCAC",1)
         
        
        

















