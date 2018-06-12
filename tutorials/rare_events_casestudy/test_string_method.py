
from hybridization23 import Settings, ResultsHybridization, suyamaT, suyamaC, enum_hybridization, title_hybridization, testSeq, NUM_OF_REPEATS, CONVERGENCE_CRIT
from multistrand.system import SimSystem
from multistrand.builder import hybridizationString, Builder, BuilderRate
from multistrand.options import Options, Literals
from multistrand.experiment import standardOptions
from multistrand.objects import StopCondition

import sys, os, time
import numpy as np

RESULT_DIR = "test_string_method"
A_TIME_OUT = 4e-3

morrison = ["TTGGTGATCC", "AGATTAGCAGGTTTCCCACC", "GCCCACACTCTTACTTATCGACT", "AGAGGCTTATAACTGTGTCGGGT", "TGTTCTAAGATTATCCTCCCGCC", "GGCGGCTATAACAATTTCATCCA"]

morrison0 = ["TTGGTGATCC"]
morrison1 = ["AGATTAGCAGGTTTCCCACC"]

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


def timings(seq, nTrials):

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
        myBuilder.genUntilConvergenceWithInitialState(CONVERGENCE_CRIT, startStates[:(len(startStates) - 1)], printMeanTime=True)

        builderRate = BuilderRate(myBuilder)
        output.buildTime.append(time.time() - startTime)
     
        startTime = time.time()
        output.rates.append(np.log10(1.0 / builderRate.averageTimeFromInitial(bimolecular=True)))
        print "Rate = %.2E, time =  %.2f" % (output.rates[-1], output.buildTime[-1])
        output.matrixTime.append(time.time() - startTime)
        output.nStates.append(len(builderRate.statespace))
     
    return output


def iterateResults(seqs, nTrials):
    
    mf = file(RESULT_DIR + "/comparison_" + toggle + "-" + str(nTrials) + ".txt", "w")
    mf.write("The timeout = %f" % A_TIME_OUT)
    
    for seq in seqs:
        mf.write(seq + "\n")
        
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
    
    print sys.argv

    if len(sys.argv) > 2:
        
        toggle = sys.argv[1]
        nTrials = sys.argv[2]

    else :
        raise ValueError("Expected two input arguments, example:   morrison 500  ")

    if toggle == "morrison":
    
        iterateResults(morrison, nTrials)
        
    if toggle == "morrison0":
    
        iterateResults(morrison0, nTrials)
    
    if toggle == "morrison1":
        
        iterateResults(morrison1, nTrials)

    if toggle == "test":

        iterateResults(testSeq, nTrials)

    if toggle == "longdomain":

        iterateResults(longdomain, nTrials)


#             
        
