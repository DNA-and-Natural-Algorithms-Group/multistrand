# FD, May 17th, 2017. 
# This demonstrates the pair-type functionality
# For a given complex, Pairtype returns a unique representation. 
# This is important for hassing functions 

import sys
from multistrand.objects import Complex, Domain, Strand, StopCondition
from multistrand.options import Options
from multistrand.system import SimSystem
from multistrand.utils import pairType

    
ATIME_OUT = 0.0010


def printTrajectory(o):
    
    print "in print_trajectory"
    print o.start_state[0].structure
    print o.start_state[1].structure
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
        
        

        print tubestruct + (' t=%11.9f microseconds, ids = %s  , dG=%6.2f kcal/mol' % (time, identities2, dG)) 



def showTrajectory(strandSeq, numTraj=2):


    onedomain = Domain(name="itall", sequence=strandSeq)
    top = Strand(name="top", domains=[onedomain])
    bot = top.C

    start_complex_top = Complex(strands=[top], structure=".")
    start_complex_bot = Complex(strands=[bot], structure=".")
     
    # Turns Boltzmann sampling on for these complexes.
    start_complex_top.boltzmann_count = numTraj
    start_complex_bot.boltzmann_count = numTraj

    start_complex_top.boltzmann_sample = True

    # Stop when the exact full duplex is achieved.
    success_complex = Complex(strands=[top, bot], structure="(+)")
    success_stop_condition = StopCondition("SUCCESS", [(success_complex, Options.exactMacrostate, 0)])

    o = Options(temperature=310.15 - 5.0,
                simulation_time=ATIME_OUT,
                num_simulations=numTraj,  # don't play it again, Sam
                output_interval=1,  # record every single step
                rate_method= Options.metropolis, 
                simulation_mode = Options.trajectory
                )  

    o.start_state = [start_complex_top, start_complex_bot]
    o.stop_conditions = [success_stop_condition]

    s = SimSystem(o)
    s.start()
    printTrajectory(o)        
        

        

def run_sims():

    showTrajectory("GCGTTTCGC", 2)

    
# # The actual main method
if __name__ == '__main__':
    
    print sys.argv
    run_sims()
         
        
        

















