# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from hybridization23 import \
    Settings, ResultsHybridization, suyamaT, suyamaC, \
    enum_hybridization, title_hybridization, NUM_OF_REPEATS
from multistrand.builder import \
    hybridizationString, Builder, BuilderRate, \
    threewaybmString, dissociationString
from multistrand.options import Literals
from multistrand.experiment import standardOptions
from multistrand.objects import StopCondition

import sys, os, time
import numpy as np
import scipy

RESULT_DIR = "test_delta_pruning"
A_TIME_OUT = 4e-3

EXP_TOGGLE = 1

morrison = ["TTGGTGATCC", "AGATTAGCAGGTTTCCCACC", "GCCCACACTCTTACTTATCGACT", "AGAGGCTTATAACTGTGTCGGGT", "TGTTCTAAGATTATCCTCCCGCC", "GGCGGCTATAACAATTTCATCCA"]

morrison0 = ["TTGGTGATCC"]
morrison1 = ["AGATTAGCAGGTTTCCCACC"]
morrison15 = ["AGATTAGCAGGTTTC"]
morrison13 = ["AGATTAGCAGGTT"]

longdomain = ["AGAGGCTTATAACTGTGTCGGGT"]

TEST_RATE_LIMIT = True

""" We will omit the starting state, somewhat unusual"""

ADD_TRANSITIONS = True;

print("Scipy version = " + str(scipy.__version__))


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
   
#     stdOptions.DNA23Arrhenius()
   
    return stdOptions


def timings(seq, nTrials, deltaPruning=None):

    output = ResultsHybridization()

    def association_comparison(seq, endComplex):
        return Settings(associationNoInit, [nTrials, suyamaT, suyamaC, endComplex], enum_hybridization, title_hybridization)

    if EXP_TOGGLE == 1:
        startStates = hybridizationString(seq)
    if EXP_TOGGLE == 2:
        startStates = threewaybmString("A", seq, "ACTAGG")
    if EXP_TOGGLE == 3:            
        startStates = dissociationString(seq)
    
    endState = startStates[-1]
    sett = association_comparison(seq, endState[0])

    for i in range(NUM_OF_REPEATS):
        startTime = time.time()
        
        BuilderRate.solveToggle = 0
        myBuilder = Builder(sett.function, sett.arguments)
        Builder.verbosity = True

        myBuilder.genAndSavePathsFromString(startStates[:(len(startStates) - 1)])
        print(myBuilder)
        
        startTime2 = time.time()
        myBuilder.fattenStateSpace()
        print("Searching transitions took " + str(  1000 * (time.time() - startTime2) / len(myBuilder.protoSpace)) + "ms per state")
        print(myBuilder)
        
#         myBuilder.genUntilConvergenceWithInitialState(10000, startStates[:(len(startStates) - 1)], printMeanTime=True)

        if not deltaPruning == None:
            print("Going to delta prune with %.2E" % deltaPruning)
            myBuilder.deltaPruning(deltaPruning, printCount=True)

        maxRange = 1
        
        if TEST_RATE_LIMIT:
            maxRange = 5

        for i in range(maxRange):
            builderRate = BuilderRate(myBuilder)

            builderRate.rateLimit = builderRate.rateLimit * (10 ** i)
            if i == (maxRange - 1):
                builderRate.rateLimit = 0.0
                
            print("rateLimit = " + str(builderRate.rateLimit))
                
            builderRate.setMatrix()
            output.buildTime.append(time.time() - startTime)
         
            startTime = time.time()
            
            biCheck = (EXP_TOGGLE == 1 or EXP_TOGGLE == 2)
            
            output.rates.append(np.log10(1.0 / builderRate.averageTimeFromInitial(bimolecular=biCheck)))
            pruned_mfpt = builderRate.averageTimeFromInitial(bimolecular=False)
            print(f"Rate = {output.rates[-1]:.2E}, MFPT = {pruned_mfpt:.2E}, "
                  f"compute_time =  {output.buildTime[-1]:.2f} \n\n")
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
    
    print(sys.argv)
    deltaPruning = None

    if len(sys.argv) > 2:
        toggle = sys.argv[1]
        nTrials = sys.argv[2]
    else:
        raise ValueError("Expected two input arguments, example:   test 500  ")

    if len(sys.argv) > 3:
        deltaPruning = float(sys.argv[3])

    if toggle == "morrison":
        iterateResults(morrison, nTrials, deltaPruning=deltaPruning)
    if toggle == "morrison0"  or toggle == "test":
        iterateResults(morrison0, nTrials, deltaPruning=deltaPruning)
    if toggle == "morrison1":
        iterateResults(morrison1, nTrials, deltaPruning=deltaPruning)
    if toggle == "morrison15":
        iterateResults(morrison15, nTrials, deltaPruning=deltaPruning)
    if toggle == "morrison13":
        iterateResults(morrison13, nTrials, deltaPruning=deltaPruning)
    if toggle == "longdomain":
        iterateResults(longdomain, nTrials, deltaPruning=deltaPruning)
