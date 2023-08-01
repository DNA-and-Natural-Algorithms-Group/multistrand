# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
This demonstrates the co-transcriptional folding.
"""

import time

from multistrand.objects import Complex, Domain, Strand
from multistrand.options import Literals
from multistrand.system import SimSystem
from multistrand.utils import printTrajectory
from multistrand.experiment import standardOptions


def doSims(strandSeq, numTraj=1):
    curr = time.time()

    o = standardOptions(tempIn = 36.95)
    o.num_simulations = numTraj
    o.output_time = 0.0004       #       output every .4 ms
    o.simulation_time = 0.35     #       unit: second
    o.gt_enable = 1
    o.substrate_type = Literals.substrateRNA
    o.simulation_mode = Literals.trajectory
    o.cotranscriptional = True; # enables the strand growing on the 3' end.

    onedomain = Domain(name="mydomain", sequence=strandSeq)
    top = Strand(name="top", domains=[onedomain])
    startTop = Complex(strands=[top], structure=".")
    o.start_state = [startTop]
    o.DNA23Arrhenius()
#     o1.initial_seed = 1777+6

    s = SimSystem(o)
    s.start()
    printTrajectory(o)
    print(f"Exe time is {time.time() - curr}")

    
if __name__ == '__main__':
    doSims("ATTCCGGTTGATCCTGCCGGAGGTCATTGCTATTGGGGTCCGATTTAGCCATGCTAGTTGCACGAGTTCATACTCGTGGCGAAAAGCTCAGTAACACGTGGCCAAACTACCCTACAGAGAA")
    print()
    doSims("ATTCCGGTTGATCCTGCCGGAGGTCATTGCTATTGGGGTCCGATTTAGCCATGCTAGTTGCACGAGTTCATACTCGTGGCGAAAAGCTCAGTAACACGTGGCCAAACTACCCTACAGAGAACGATAACCTCGGGAAACTGAGGCTAATAGTTCATACGGGAGTCATGCTGGAATGCCGACTCCCCGAAACGCTCAGGCGCTGTAGGATGTGGCTGCGGCCGATTAGGTAGACGGTGGGGTAACGGCCCACCGTGCCGATAATCGGTACGGGTTGTGAGAGCAAGAGCCCGGAGACGGAATCT")
