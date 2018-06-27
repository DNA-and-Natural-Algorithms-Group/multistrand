
from hybridization23 import Settings, ResultsHybridization, suyamaT, suyamaC, enum_hybridization, title_hybridization, testSeq, NUM_OF_REPEATS, CONVERGENCE_CRIT
from multistrand.system import SimSystem
from multistrand.builder import hybridizationString, Builder, BuilderRate
from multistrand.options import Options, Literals
from multistrand.experiment import standardOptions
from multistrand.objects import StopCondition

import sys, os, time
import numpy as np

RESULT_DIR = "test_delta_pruning"
A_TIME_OUT = 4e-3

morrison = ["TTGGTGATCC", "AGATTAGCAGGTTTCCCACC", "GCCCACACTCTTACTTATCGACT", "AGAGGCTTATAACTGTGTCGGGT", "TGTTCTAAGATTATCCTCCCGCC", "GGCGGCTATAACAATTTCATCCA"]

morrison0 = ["TTGGTGATCC"]
morrison1 = ["AGATTAGCAGGTTTCCCACC"]
morrison15 = ["AGATTAGCAGGTTTC"]

longdomain = ["AGAGGCTTATAACTGTGTCGGGT"]

""" We will omit the starting state, somewhat unusual"""


def associationNoInit(arguments): 

    stdOptions = standardOptions()
    stdOptions.simulation_mode = Literals.trajectory
    stdOptions.verbosity = True

    stdOptions.num_simulations = arguments[0]
    stdOptions.temperature = arguments[1]
    stdOptions.join_concentration = arguments[2]
    
    endComplex = arguments[3]
    
    stdOptions.simulation_time = A_TIME_OUT
    
    stopSuccess = StopCondition(Literals.success, [(endComplex, Literals.exact_macrostate, 0)])
    stdOptions.stop_conditions = [stopSuccess]
   
    stdOptions.output_interval = 1
    
    return stdOptions


def timings(seq, nTrials, deltaPruning=None):

    output = ResultsHybridization()

    def association_comparison(seq, endComplex):
        return Settings(associationNoInit, [nTrials, suyamaT, suyamaC, endComplex], enum_hybridization, title_hybridization)

    startStates = hybridizationString(seq)
    
    endState = startStates[-1]
    # startState = startStates[0]

    sett = association_comparison(seq, endState[0])

    for i in range(NUM_OF_REPEATS):
    
        startTime = time.time()
        
        myBuilder = Builder(sett.function, sett.arguments)
        myBuilder.genUntilConvergenceWithInitialState(CONVERGENCE_CRIT , startStates[:(len(startStates) - 1)], printMeanTime=True)

        if not deltaPruning == None:
            print "Going to delta prune with %.2E" %deltaPruning
            myBuilder.deltaPruning(deltaPruning, printCount=True)

        builderRate = BuilderRate(myBuilder)
        output.buildTime.append(time.time() - startTime)
     
        startTime = time.time()
        output.rates.append(np.log10(1.0 / builderRate.averageTimeFromInitial(bimolecular=True)))
        pruned_mfpt = builderRate.averageTimeFromInitial(bimolecular = False )
        print "Rate = %.2E, MFPT = %.2E, compute_time =  %.2f \n\n " % (output.rates[-1], pruned_mfpt,  output.buildTime[-1])
        output.matrixTime.append(time.time() - startTime)
        output.nStates.append(len(builderRate.statespace))
     
    return output


def iterateResults(seqs, nTrials, deltaPruning=None):
    
    outPath = RESULT_DIR + "/comparison_" + toggle + "-" + str(nTrials) + ".txt"
    
    if not deltaPruning == None:
        outPath + "dp=%.2E" % deltaPruning
    
    mf = file(outPath, "w")
    mf.write("The timeout = %f" % A_TIME_OUT)
    
    if not deltaPruning == None:
        mf.write("Delta-pruning enabled at %.2E" % deltaPruning)
    else:
        mf.write("Delta-pruning disabled.")
    
    for seq in seqs:
        mf.write(seq + "\n")
        
        result = timings(seq, nTrials, deltaPruning)
#         print result
        mf.write(str(result))
        mf.write("\n\n")
        mf.flush()
        
    mf.close()


# # The actual main method
if __name__ == '__main__':
    
    if not os.path.exists(RESULT_DIR):
        os.makedirs(RESULT_DIR)
    
    print sys.argv
    deltaPruning = None

    if len(sys.argv) > 2:
        
        toggle = sys.argv[1]
        nTrials = sys.argv[2]

    else :
        raise ValueError("Expected two input arguments, example:   morrison 500 0.01 ")

    if len(sys.argv) > 3:
        
        deltaPruning = float(sys.argv[3])

    if toggle == "morrison":
    
        iterateResults(morrison, nTrials, deltaPruning=deltaPruning)
        
    if toggle == "morrison0":
    
        iterateResults(morrison0, nTrials, deltaPruning=deltaPruning)
    
    if toggle == "morrison1":
        
        iterateResults(morrison1, nTrials, deltaPruning=deltaPruning)

    if toggle == "morrison15":
        
        iterateResults(morrison15, nTrials, deltaPruning=deltaPruning)

    if toggle == "test":

        iterateResults(testSeq, nTrials, deltaPruning=deltaPruning)

    if toggle == "longdomain":

        iterateResults(longdomain, nTrials, deltaPruning=deltaPruning)

#             
        
