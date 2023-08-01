# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
This demonstrates the pair-type functionality
For a given complex, Pairtype returns a unique representation.
This is important for hashing functions
"""

from multistrand.options import Literals
from multistrand.system import SimSystem
from multistrand.utils import printTrajectory, dGC_feature
from multistrand.experiment import standardOptions, makeComplex


def doSims(numTraj=2):

    o = standardOptions()
    
    o.simulation_mode = Literals.trajectory
    o.num_simulations = numTraj
    o.output_interval = 1
    o.join_concentration = 1.0e-9
#     o1.join_concentration = 1.0
    o.simulation_time = 0.000005
    o.DNA23Metropolis()

    seq1 = "ACTGACTGACTG"
    seq2 = "ACTG"
    seq3 = "CATTCAGTACAGT"
    
    struct = "((((((.((.(.+).((+)).)).)).))))"
    
#     # Using domain representation makes it easier to write secondary structures.
#     domain1 = Domain(name="d1", sequence=seq1)
#     domain2 = Domain(name="d2", sequence=seq2)
#     domain3 = Domain(name="d3", sequence=seq3)
#     
#     strand1 = Strand(name="s1", domains=[domain1])
#     strand2 = Strand(name="s2", domains=[domain2])
#     strand3 = Strand(name="s3", domains=[domain3])

    myComplex = makeComplex([seq1, seq2, seq3], struct)
    print(myComplex)

    o.start_state = [myComplex]
    # Note: no stopping conditions
    o.initial_seed = 1777+6
    s = SimSystem(o)
    s.start()
    printTrajectory(o, timescale=(1e9, "ns"), feature=dGC_feature)

    
if __name__ == '__main__':
    doSims(1)
