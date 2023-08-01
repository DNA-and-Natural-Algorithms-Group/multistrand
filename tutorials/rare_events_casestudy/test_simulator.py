# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from hybridization23 import \
    Settings, ResultsHybridization, suyamaT, suyamaC, \
    enum_hybridization, title_hybridization, testSeq, NUM_OF_REPEATS, \
    doReactionAssociation
from multistrand.concurrent import FirstPassageRate
from multistrand.system import SimSystem

import sys, os, time
import numpy as np

RESULT_DIR = "test_simulator"
A_TIME_OUT = 2e-1
CONVERGENCE_CRIT = 0.01


def timings(seq, nTrials):

    output = ResultsHybridization()

    def association_comparison(seq):
        return Settings(doReactionAssociation, [nTrials, suyamaT, suyamaC, seq],
                        title_hybridization, enum_hybridization)

    sett = association_comparison(seq)
    for i in range(NUM_OF_REPEATS):
        startTime = time.time()
        o1 = sett.function(sett.arguments)

        ssystem = SimSystem(o1)
        ssystem.start()
        
        myRates = FirstPassageRate(o1.interface.results)
        k1 = myRates.kEff(o1.join_concentration)

        output.buildTime.append(time.time() - startTime)
        output.rates.append(np.log10(k1))
        output.nStates.append(100)
        output.matrixTime.append(10.0)
    return output


def iterateResults(seqs, nTrials):
    
    with open(RESULT_DIR + "/comparison_" + toggle + "-" + str(nTrials) +".txt", "w") as mf:
        for seq in seqs:
            mf.write(seq + "\n");

            result = timings(seq, nTrials)
    #         print result
            mf.write(str(result))
            mf.write("\n\n")
            mf.flush()


# # The actual main method
if __name__ == '__main__':
    
    if not os.path.exists(RESULT_DIR):
        os.makedirs(RESULT_DIR)
    
    if len(sys.argv) > 2:
        toggle = sys.argv[1]
        nTrials = sys.argv[2]
    else:
        raise ValueError("Expected two input arguments, example:   morrison 500  ")

    if toggle == "morrison":
        iterateResults(morrison, nTrials)
    if toggle == "test":
        iterateResults(testSeq, nTrials)
