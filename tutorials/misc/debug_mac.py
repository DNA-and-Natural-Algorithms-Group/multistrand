# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
This demonstrates the pair-type functionality
For a given complex, Pairtype returns a unique representation.
This is important for hashing functions.
"""

import time

from multistrand.options import Literals
from multistrand.system import SimSystem
from multistrand.utils import printTrajectory
from multistrand.experiment import standardOptions, makeComplex

    
def doSims(numTraj=1):
    o = standardOptions()

    o.simulation_mode = Literals.first_step
    o.num_simulations = numTraj
    o.output_interval = 1
    o.simulation_time = 0.000001
    
    myComplex1 = makeComplex(["GTCACTGCTTTT"], "............")
    myComplex2 = makeComplex(["GTCACTGC", "GCAGTGAC"], ".(((((((+))))))).")
    o.start_state = [myComplex1, myComplex2]
    # myComplex = makeComplex(["GTCACTGCTTTT","GTCACTGC","GCAGTGAC"], "..(.........+).((((((+))))))..")
    # o.start_state = [myComplex]

    # no stop conditions, just timeout.
    o.initial_seed = time.time() * 1000000
    o.verbosity = 3
    # o.initial_seed = 1501710503137097

    s = SimSystem(o)
    s.start()
    printTrajectory(o)
    print()
    print()
    return o, s


def main():
    o, s = None, None
    for _ in range(1000):
        time.sleep(0.000001)
        o, s = doSims()
    return o, s


if __name__ == '__main__':
    main()
