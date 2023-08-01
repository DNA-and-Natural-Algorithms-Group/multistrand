# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
This file demonstrates Arrhenius kinetics.
The interface is not finalized.
"""

from multistrand._options.options import Options
from multistrand.experiment import standardOptions, hybridization
from multistrand.builder import codeToDesc
from multistrand.system import SimSystem
from multistrand.utils import printTrajectory


def arr_move_type(o: Options, i: int):
    desc = codeToDesc(o.full_trajectory_arrType[i])
    return f"Type: {desc[0]} {desc[1]}"


def doSims(strandSeq, numTraj=1):
    o = standardOptions()
    o.num_simulations = numTraj
    o.output_interval = 1
    o.simulation_time = 0.001
    hybridization(o, strandSeq)
    o.DNA23Arrhenius()
    o.initial_seed = 1777 + 6

    s = SimSystem(o)
    s.start()
    printTrajectory(o, feature=arr_move_type)

    return o, s

    
def main():
    return doSims("GCGTTTCAC")


if __name__ == '__main__':
    main()
