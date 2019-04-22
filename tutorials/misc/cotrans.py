# FD, Oct 20th, 2017. 
# This demonstrates the co-transcriptional folding.

import sys, time

from multistrand.objects import Complex, Domain, Strand, StopCondition
from multistrand.options import Options, Literals
from multistrand.system import SimSystem
from multistrand.utils import pairType
from multistrand.experiment import standardOptions, hybridization

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

    curr = time.time()

    o1 = standardOptions(tempIn = 36.95)
    
    o1.num_simulations = numTraj
    o1.output_time = 0.0004       #       output every .4 ms
    o1.simulation_time = 0.5      #       unit: second
    o1.gt_enable = 1;
    o1.substrate_type = Literals.substrateRNA
    o1.simulation_mode = Literals.trajectory
    o1.cotranscriptional = True; # enables the strand growing on the 3' end.

    onedomain = Domain(name="mydomain", sequence=strandSeq) 
    
    top = Strand(name="top", domains=[onedomain])
    startTop = Complex(strands=[top], structure=".")
    o1.start_state = [startTop]
    o1.DNA23Arrhenius()

#     o1.initial_seed = 1777+6

    s = SimSystem(o1)
    s.start()
    printTrajectory(o1)     
    
    print "Exe time is " + str(time.time() - curr)   
        

    
# # The actual main method
if __name__ == '__main__':
    
    print sys.argv
    doSims( "ATTCCGGTTGATCCTGCCGGAGGTCATTGCTATTGGGGTCCGATTTAGCCATGCTAGTTGCACGAGTTCATACTCGTGGCGAAAAGCTCAGTAACACGTGGCCAAACTACCCTACAGAGAA",1)
         
        
        






#    onedomain = Domain(name="itall", sequence="ATTCCGGTTGATCCTGCCGGAGGTCATTGCTATTGGGGTCCGATTTAGCCATGCTAGTTGCACGAGTTCATACTCGTGGCGAAAAGCTCAGTAACACGTGGCCAAACTACCCTACAGAGAACGATAACCTCGGGAAACTGAGGCTAATAGTTCATACGGGAGTCATGCTGGAATGCCGACTCCCCGAAACGCTCAGGCGCTGTAGGATGTGGCTGCGGCCGATTAGGTAGACGGTGGGGTAACGGCCCACCGTGCCGATAATCGGTACGGGTTGTGAGAGCAAGAGCCCGGAGACGGAATCTGAGACAAGATTCCGGGCCCTACGGGGCGCAGCAGGCGCGAAACCTTTACACTGCACGCAAGTGCGATAAGGGGACCCCAAGTGCGAGGGCATATAGTCCTCGCTTTTCTCGACCGTAAGGCGGTCGAGGAATAAGAGCTGGGCAAGACCGGTGCCAGCCGCCGCGGTAATACCGGCAGCTCAAGTGATGA") #CCGATATTATTGGGCCTAAAGCGTCCGTAGCCGGCCACGAAGGTTCATCGGGAAATCCGCCAGCTCAACTGGCGGGCGTCCGGTGAAAACCACGTGGCTTGGGACCGGAAGGCTCGAGGGGTACGTCCGGGGTAGGAGTGAAATCCCGTAATCCTGGACGGACCACCGATGGCGAAAGCACCTCGAGAAGACGGATCCGACGGTGAGGGACGAAAGCTAGGGTCTCGAACCGGATTAGATACCCGGGTAGTCCTAGCTGTAAACGATGCTCGCTAGGTGTGACACAGGCTACGAGCCTGTGTTGTGCCGTAGGGAAGCCGAGAAGCGAGCCGCCTGGGAAGTACGTCCGCAAGGATGAAACTTAAAGGAATTGGCGGGGGAGCACTACAACCGGAGGAGCCTGCGGTTTAATTGGACTCAACGCCGGACATCTCACCAGCTCCGACTACAGTGATGACGATCAGGTTGATGACCTTATCACGACGCTGTAGAGAGGAGGTGCATGGCCGCCGTCAGCTCGTACCGTGAGGCGTCCTGTTAAGTCAGGCAACGAGCGAGACCCGCACTTCTAATTGCCAGCAGCAGTTTCGACTGGCTGGGTACATTAGAAGGACTGCCGCTGCTAAAGCGGAGGAAGGAACGGGCAACGGTAGGTCAGTATGCCCCGAATGAGCTGGGCTACACGCGGGCTACAATGGTCGAGACAATGGGTTGCTATCTCGAAAGAGAACGCTAATCTCCTAAACTCGATCGTAGTTCGGATTGAGGGCTGAAACTCGCCCTCATGAAGCTGGATTCGGTAGTAATCGCATTTCAATAGAGTGCGGTGAATACGTCCCTGCTCCTTGCACACACCGCCCGTCAAAGCACCCGAGTGAGGTCCGGATGAGGCCACCACACGGTGGTCGAATCTGGGCTTCGCAAGGGGGCTTAAGTCGTAACAAGGTAGCCGTAGGGGAATCTGCGGCTGGATCACCTCCTG")    
 










