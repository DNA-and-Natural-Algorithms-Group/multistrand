# FD, May 17th, 2017. 
# This demonstrates the pair-type functionality
# For a given complex, Pairtype returns a unique representation. 
# This is important for hassing functions 


from multistrand.system import SimSystem
from multistrand.utils import printTrajectory
from multistrand.experiment import standardOptions, hybridization

    
def doSims(strandSeq, numTraj=2):
    o1 = standardOptions()
    o1.num_simulations = numTraj
    o1.output_interval = 1 
    o1.simulation_time = 0.001
    hybridization(o1, strandSeq )
    o1.initial_seed = 1777+6

    s = SimSystem(o1)
    s.start()
    printTrajectory(o1)        

    
if __name__ == '__main__':
    doSims("GCGT", 1)
