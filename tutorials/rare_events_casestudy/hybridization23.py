# -*- coding: utf-8 -*-

""" Frits Dannenberg, Jun 3rd 2018 
    Feasability study for builder code for experiments in 
    
    Influence of thermodynamically unfavorable
    secondary structures on DNA hybridization kinetics
    Hiroaki Hata, Tetsuro Kitajima and Akira Suyama         """

import sys, time, os

from multistrand.builder import Builder, BuilderRate

from multistrand.options import Options, Literals
from multistrand.experiment import standardOptions, hybridization

import numpy as np

A_TIME_OUT = 1000.0
RESULT_DIR = "suyama"
NUM_OF_REPEATS = 6
DO_CONVERGENCE = True
CONVERGENCE_CRIT = 0.05


enum_hybridization = "hybridization"
title_hybridization = "Hybridization of "

suyamaT = 25.0
suyamaC = 50e-7
suyamaNa = 0.195


# suyama16 = ["GCCCACACTCTTACTTATCGACT", "AGAGGCTTATAACTGTGTCGGGT", "TGTTCTAAGATTATCCTCCCGCC", "GGCGGCTATAACAATTTCATCCA", "TAGCCCAGTGATTTATGACATGC", "GCATCTACACTCAATACCCAGCC",
#             "GCCCGTACTGTTGAGATTATGGT", "GCACCTCCAAATAAAAACTCCGC", "AGATCAGAGATAGTTACGCCGCA", "TATGTTCCTTACCCCGTTTACCA", "TAGCCAACTCTAAATAACGGACG", "GAAGGAATGTTAAAATCGTCGCG",
#             "TTTGTTTTCCTTATGAGCCAGCC", "GCCCCGATATCTATTTTAGGACG", "CGCAGGAGAGTTAAACGAAAGCA", "GGCTCTATACGATTAAACTCCCC"]

suyama16 = ["GCCCACACTCTTACTTATCGACT", "AGAGGCTTATAACTGTGTCGGGT", "TGTTCTAAGATTATCCTCCCGCC", "GGCGGCTATAACAATTTCATCCA"] #"TAGCCCAGTGATTTATGACATGC", "GCATCTACACTCAATACCCAGCC"]

suyama0 = "GCCCACACTCTTACTTATCGACT"
suyama1 = "AGAGGCTTATAACTGTGTCGGGT"

morrison = ["TTGGTGATCC", "AGATTAGCAGGTTTCCCACC"]


test2 = ["GCCCACACGC", "AGAGGCTGC"]

""" 
    results holds three arrays:
    .rates            - the computed rate
    .buildTime        - Time to simulate, create Builder and BuilderRate 
    .matrixTime       - Time to solve the matrix.                
"""
class resultsHybridization(object):
    
    def __init__(self):
    
        self.rates = []
        self.buildTime = []
        self.matrixTime = []
        self.nStates = []
        
    def __str__(self):
        
        ratesSD = np.std(self.rates)
        buildSD = (np.std(self.buildTime), 100.0 *np.std(self.buildTime)*NUM_OF_REPEATS/sum(self.buildTime))
        matrixSD = (np.std(self.matrixTime), 100.0 * np.std(self.matrixTime)*NUM_OF_REPEATS/sum(self.matrixTime))    
        nStatesSD = (np.std(self.nStates), 100.0 * np.std(self.nStates)*NUM_OF_REPEATS/sum(self.nStates))
        
        averages = (sum(self.rates) / NUM_OF_REPEATS , sum(self.buildTime) / NUM_OF_REPEATS , sum(self.matrixTime) / NUM_OF_REPEATS, sum(self.nStates) / NUM_OF_REPEATS  )
        
        output = "Rates, build times, matrix solve times, number of states \n"
        output += "%.2f   -   %.2f   -   %.2f   -  %d  \n" % averages
        output += ",     ".join(["%.2f" % x for x in self.rates]) + "                s.d. = %0.2f,       " % ratesSD   + "\n"
        output += ",     ".join(["%.2f" % x for x in self.buildTime]) + "            s.d. = %0.2E,       rel s.d. = %.1f %%" % buildSD   + "\n"
        output += ",     ".join(["%.2f" % x for x in self.matrixTime]) + "           s.d. = %0.2E,       rel s.d. = %.1f %%" % matrixSD   + "\n"
        output += ",     ".join(["%d" % x for x in self.nStates]) + "                s.d. = %0.2E,       rel s.d. = %.1f %%" % nStatesSD   + "\n"
        
        
        output += ""
        
        return output

class Settings(object):
    
    def __init__(self, function, arguments, title="", enum="", tempOverMelt=0.0, multipleC=1.0):
    
        self.type = enum
        self.title = title

        self.function = function
        self.arguments = arguments  # often, just sequences
        
        self.bimolecular = False
        
        self.multipleC = multipleC
        self.tempOverMelt = tempOverMelt

    def filename(self):
        
        output = self.type + "_" + str(self.arguments) + "_T=" 
        
        if DO_CONVERGENCE:
            output += "convergence" + str(CONVERGENCE_CRIT)
        
        if not self.tempOverMelt == 0.0:
            output += "_" + "_dT=" + "%.1f" % self.tempOverMelt
        
        if not self.multipelC == 1.0:
            output += " mulC=" + str(self.multipleC)
        
        return output

    def __str__(self):
        
        return self.filename()


def doReactionAssociation(arguments, output=True): 

#     print str(arguments)

    stdOptions = standardOptions()

    stdOptions.simulation_mode = Literals.trajectory
    
    stdOptions.num_simulations = arguments[0]
    stdOptions.temperature = arguments[1]
    stdOptions.join_concentration = arguments[2]
    
    seq = arguments[3]
    hybridization(stdOptions, seq, doFirstPassage =True)

    if len(arguments) > 4:
            stdOptions.sodium = arguments[4]
            stdOptions.magnesium = arguments[5]

    if output:
        stdOptions.output_interval = 1
        
    stdOptions.simulation_time = A_TIME_OUT
    stdOptions.simulation_mode = Literals.trajectory

    i = 65
    
    for state in stdOptions.start_state:
        for strand in state.strand_list:
            strand.id = i
            i += 1
    
    return stdOptions


def timings(settings):
    
    output = resultsHybridization()
    
    for i in range(NUM_OF_REPEATS):
        
#             print "Starting to make the statespace \n"
            startTime = time.time()

            myBuilder = Builder(settings.function, settings.arguments)
            
            if DO_CONVERGENCE:
                myBuilder.genUntilConvergence(CONVERGENCE_CRIT)
            else:            
                myBuilder.genAndSavePathsFile()
                
            builderRate = BuilderRate(myBuilder) 
            
            output.buildTime.append(time.time() - startTime)

#             print "Solving the matrix equation \n"       
            startTime = time.time()
            
            output.rates.append(np.log10(1.0 / builderRate.averageTimeFromInitial(bimolecular=True)))
            output.matrixTime.append(time.time() - startTime)
            
            output.nStates.append(len(builderRate.statespace))
            

    
    return output
    
    """Produce the required plots here"""


def genFile(mySeqs, nTrials, toggle):
    
    def association_comparison(seq):
        return Settings(doReactionAssociation, [nTrials, suyamaT, suyamaC, seq], enum_hybridization, title_hybridization)
    
    mf = file(RESULT_DIR + "/comparison_" + toggle + "-" + str(nTrials) +".txt", "w")
    
    for seq in mySeqs:
        
        mf.write(seq + "\n");
        results = timings(association_comparison(seq))
        
#         print results.rates
        
        mf.write(str(results))
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

    else:
        raise ValueError( "Expected two input arguments, example:   test 500  ")
    
    if toggle == "test":
        genFile(test2, nTrials, toggle)
        
    elif toggle == "suyama":
        genFile(suyama16, nTrials, toggle)

    elif toggle == "suyama0":
        genFile([suyama0], nTrials, toggle)
        
    elif toggle == "suyama1":
        genFile([suyama1], nTrials, toggle)
        
    elif toggle == "morrison":
        genFile(morrison, nTrials, toggle)

    
