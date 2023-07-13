
from hybridization23 import Settings, ResultsHybridization, suyamaT, suyamaC, enum_hybridization, title_hybridization, testSeq, NUM_OF_REPEATS, doReactionAssociation
from multistrand.concurrent import FirstPassageRate
from multistrand.system import SimSystem
from multistrand.builder import hybridizationString, Builder, BuilderRate
from multistrand.options import Options, Literals
from multistrand.experiment import standardOptions
from multistrand.objects import StopCondition

import sys, os, time
import numpy as np

RESULT_DIR = "test_simulator"
A_TIME_OUT = 2e-1
CONVERGENCE_CRIT = 0.01


def timings(seq, nTrials):

    output = ResultsHybridization()

    def association_comparison(seq):
        return Settings(doReactionAssociation, [nTrials, suyamaT, suyamaC, seq], enum_hybridization, title_hybridization)

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
    
    mf = file(RESULT_DIR + "/comparison_" + toggle + "-" + str(nTrials) +".txt", "w")
    
    for seq in seqs:
        mf.write(seq + "\n");
        
        result = timings(seq, nTrials)
#         print result
        mf.write(str(result))
        mf.write("\n\n")
        mf.flush()
        
        
    mf.close()


# # The actual main method
if __name__ == '__main__':
    
    if not os.path.exists(RESULT_DIR):
        os.makedirs(RESULT_DIR)
    
    print(sys.argv)

    if len(sys.argv) > 2:
        
        toggle = sys.argv[1]
        nTrials = sys.argv[2]

    else :
        raise ValueError("Expected two input arguments, example:   morrison 500  ")

    if toggle == "morrison":
    
        iterateResults(morrison, nTrials)
    
    if toggle == "test":

        iterateResults(testSeq, nTrials)

#             
        
