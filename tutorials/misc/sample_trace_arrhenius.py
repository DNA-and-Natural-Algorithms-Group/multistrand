"""
    Frits Dannenberg May 29th, 2018 
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


def doSims(strandSeq, numTraj=2):
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

    
if __name__ == '__main__':
    doSims("GCGTTTCAC", 1)
