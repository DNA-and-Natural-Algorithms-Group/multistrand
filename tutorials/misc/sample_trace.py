# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
This demonstrates the pair-type functionality
For a given complex, Pairtype returns a unique representation.
This is important for hashing functions
"""

from multistrand.system import SimSystem
from multistrand.utils import printTrajectory
from multistrand.experiment import standardOptions, hybridization

    
def doSims(strandSeq, numTraj=1):
    o = standardOptions()
    o.num_simulations = numTraj
    o.output_interval = 1
    o.simulation_time = 0.001
    hybridization(o, strandSeq )
    o.initial_seed = 1777+6

    s = SimSystem(o)
    s.start()
    printTrajectory(o)
    return o, s


def main():
    return doSims("GCGT")


if __name__ == '__main__':
    main()
