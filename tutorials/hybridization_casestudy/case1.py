# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
This is used for the slowdown study
and the comparison between first step and trajectory mode

Run this using arguments
plots 160
or
slowDownStudy 82000
"""

from multistrand.objects import Strand
from multistrand.experiment import standardOptions, hybridization
from multistrand.utils.utility import concentration_string, standardFileName
import multistrand.utils.thermo as thermo
from multistrand.concurrent import Bootstrap, MergeSim
from multistrand.options import Literals


from constantsgao import goa2006_P0, goa2006_P3, goa2006_P4, setSaltGao2006, colors


import matplotlib.pylab as plt
import numpy as np
import time, sys

SCRIPT_DIR = "case1_first_step"
TEMPERATURE = 20.0
ATIME_OUT = 100.0
BOOTSTRAP_RESAMPLE = 1000
POINT_SIZE = 38  # size of markers in the plot


markers = ["8", ">", "D", "s", "*", "<", "^"] 



def first_step_simulation(strand_seq, trials, T=20.0):

    myMS = MergeSim()    
    myMS.setNumOfThreads(8) 
    print(f"Running first step mode simulations for {strand_seq} (with Boltzmann sampling)...")
    
    def getOptions(trials):
        o = standardOptions(Literals.first_step, TEMPERATURE, trials, ATIME_OUT)
        hybridization(o, strand_seq, trials)
        setSaltGao2006(o)
        o.DNA23Metropolis()
        return o
    
    myMS.setOptionsFactory1(getOptions, trials)
    myMS.setFirstStepMode()  # ensure the right results object is set.
#     myMultistrand.setLeakMode()
    myMS.setTerminationCriteria(terminationCount=trials)
    myMS.run()
    return myMS


def first_passage_association(strand_seq, trials, concentration, T=20.0):

    thisMS = MergeSim()
    thisMS.setNumOfThreads(8)
    print(f"Running first passage time simulations for association of "
          f"{strand_seq} at {concentration_string(concentration)}...")
    
    def getOptions(trials):
        o = standardOptions(Literals.first_passage_time, TEMPERATURE, trials, ATIME_OUT)
        hybridization(o, strand_seq, trials, True)
        setSaltGao2006(o)
        o.join_concentration = concentration
        o.DNA23Metropolis()
        return o
    
    thisMS.setOptionsFactory1(getOptions, trials)
    thisMS.setPassageMode()
    thisMS.run()
    return thisMS


def doFirstStepMode(seq, concentrations, T=20.0, numOfRuns=500):
    
    # track time for each kind of simulation, using time.time(), which has units of second
    # do one "first step mode" run, get k1, k2, etc, from which z_crit and k_eff(z) can be computed

    myMultistrand = first_step_simulation(seq, numOfRuns, T=T) 
    myRates = myMultistrand.results
    
#     print myRates.dataset
    
    time2 = time.time()
    print(str(myRates))

    FSResult = list()
    for z in concentrations:
        kEff = myRates.kEff(z)
            
        myBootstrap = Bootstrap(myRates, N=BOOTSTRAP_RESAMPLE, concentration=z, computek1=False)
        
        low, high = myBootstrap.ninetyFivePercentiles()
        logStd = myBootstrap.logStd()

        print(f"keff = {kEff:g} /M/s at {concentration_string(z)}")
        
        myResult = (np.log10(kEff), np.log10(low), np.log10(high), logStd)
        FSResult.append(myResult)
        
    print()
    
    # call NUPACK for pfunc dG of the reaction, calculate krev based on keff
    print("Calculating dissociate rate constant based on "
          "NUPACK partition function energies and first step mode k_eff...")

    dG_top = thermo.complex_free_energy([seq], celsius=T)
    dG_bot = thermo.complex_free_energy([ Strand(sequence=seq).C.sequence ], celsius=T)
    dG_duplex = thermo.complex_free_energy([ seq, Strand(sequence=seq).C.sequence ], celsius=T)
    RT = 1.987e-3 * (T + thermo.C2K)
    time3 = time.time()
    time_nupack = time3 - time2
    krev_nupack = kEff * np.exp((dG_duplex - dG_top - dG_bot) / RT)
    print(f"krev = {krev_nupack:g} /s ({time_nupack:g} seconds)")

    times = list()
    for i in concentrations:
        myTime = (np.log10(myMultistrand.runTime), 0.0, 0.0)
        times.append(myTime)
    
    return FSResult, times
   

def doFirstPassageTimeAssocation(seq, concentrations, T=20, numOfRuns=500):

    # for each concentration z, do one "first passage time" run for association, and get k_eff(z)
    Result = []
    times = list()

    for concentration in concentrations:
 
        if len(seq) > 10 and concentration < 1e-5:
            keff = 5.5
            low = 1e+5
            high = 1e+6
            logStd = -1.0
             
        else:
            myMultistrand = first_passage_association(seq, numOfRuns, concentration=concentration, T=T)
            myRates = myMultistrand.results
            
            keff = myRates.log10KEff(concentration)
            
            myBootstrap = Bootstrap(myRates, N=BOOTSTRAP_RESAMPLE, concentration=concentration)
            low, high = myBootstrap.ninetyFivePercentiles()
            logStd = myBootstrap.logStd()
                
        Result.append((keff, np.log10(low), np.log10(high), logStd))
        times.append((np.log10(myMultistrand.runTime), 0.0, 0.0))

    print(f"keff = {keff:g} /M/s at {concentration_string(concentration)}")

    return Result, times


def fluffyPlot(ax, seqs, concentrations):
    
    myXTicks = list()
    for conc in concentrations:
        myXTicks.append(concentration_string(conc))
    
    plt.xticks(rotation=-40)
    plt.xticks(np.log10(concentrations), myXTicks)
    
    plt.gca().invert_xaxis()
    ax.set_xlabel('Concentration')
    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2, box.width * 0.8, 0.8 * box.height])
    

    ax.legend(seqs, loc='center left', bbox_to_anchor=(1, 0.5))
    

def addPoints(results, i, alp, seqs, concentrations, lineStyle, extraOptions=None):
    
#     print results
    xVal = np.log10(concentrations)
    yVal = [r[0] for r in results[i]]
    
    plt.scatter(x=xVal, y=yVal, color=colors[i], marker=markers[i], alpha=alp, s=POINT_SIZE, label=str(seqs[i]))
    plt.plot(xVal, yVal, marker=markers[i], linestyle=lineStyle, color=colors[i], alpha=alp , label=str(seqs[i]))
    
    if not extraOptions == "noErrorBars" :

        yLow = np.array([ (r[1])  for r in results[i]])
        yHigh = np.array([ (r[2]) for r in results[i]])
        plt.errorbar(xVal, yVal, yerr=[yHigh - yVal, yVal - yLow], color=colors[i], linestyle=lineStyle, alpha=alp, capsize = 6,  label=str(seqs[i]))


def doPlots(seqs, concentrations, results1, results2, trials):
       
    ax = plt.subplot(111)

    length = len(seqs)
    for i in range(length):
        addPoints(results2, i, 0.85, seqs, concentrations, '-')
    fluffyPlot(ax, seqs, concentrations)
    
    for i in range(length):
        addPoints(results1, i, 0.45, ["", "", "", "", "", "", "", ""], concentrations, '--')
    
    ax.set_title("Estimated hybridization rate (" + str(trials) + " trajectories)")
    ax.set_ylabel("k-effective (per second, log 10)")    
    plt.gca().invert_xaxis()
    plt.savefig(standardFileName(SCRIPT_DIR) + 'hybridizationRate.pdf')
    plt.close()
    
    
def doTimePlots(seqs, concentrations, results1, results2, trials):
                   
    ax = plt.subplot(111)
                       
    ax.set_title("Computation time for " + str(trials) + " trajectories")
    ax.set_ylabel("Time (seconds, log 10)")    
        
    length = len(seqs)
    for i in range(length):
        addPoints(results2, i, 0.85, seqs, concentrations, '-', extraOptions="noErrorBars")
    fluffyPlot(ax, seqs, concentrations)
    for i in range(length):
        addPoints(results1, i, 0.45, seqs, concentrations, '--', extraOptions="noErrorBars")

    plt.gca().invert_xaxis()
    plt.savefig(standardFileName(SCRIPT_DIR) + 'runTime.pdf')


def basicResults(results1, results2, runTime1, runTime2, seq, concentrations, trials):

    FSResult, times = doFirstStepMode(seq, concentrations, numOfRuns=trials)
    results1.append(FSResult)
    runTime1.append(times)

    FPResult, times = doFirstPassageTimeAssocation(seq, concentrations, numOfRuns=trials)
    results2.append(FPResult)
    runTime2.append(times)    


def doInference(concentrations, trials):

    results1 = list()
    results2 = list()
    
    runTime1 = list()
    runTime2 = list()

    seqs = list()
#     seqs.append('TCGATG')
#     seqs.append('TCGATGC')
#     seqs.append('AGTCCTTTTTGG')
    
    seqs.append('TAGTCCCTTTTTGGG')
#     seqs.append('TCGATGC')
    seqs.append('TCGATGCT')

    for seq in seqs:
        basicResults(results1, results2, runTime1, runTime2, seq, concentrations, trials)
         
    # results1, resutls2 are identical but first passage time and first step    

    doPlots(seqs, concentrations, results1, results2, trials)   
    doTimePlots(seqs, concentrations, runTime1, runTime2, trials)
   

def doSlowdownStudy(trials):
    
    def computeMeanStd(seq):
        result, times = doFirstStepMode(seq, [1.0e-6], T=20, numOfRuns=trials)
        return result[0]

    result0 = computeMeanStd(goa2006_P0)
    result3 = computeMeanStd(goa2006_P3)
    result4 = computeMeanStd(goa2006_P4)
    
    print(result0)
    print(result3)
    print(result4)
   
    output = ""
    names = ["P0", "P3", "P4"]
    myResults = [result0, result3, result4]
    output += "Rate   -   lowerbound   -   upperbound   -  logSD \n \n"
   
    for name, result in zip(names, myResults):
        output += name + "mean,  is " + str(result[0]) + "  " + str(result[1]) + "  " + str(result[2]) + "  " + str(result[3]) + "\n"

    factor3 = np.power(10, result0[0] - result3[0])
    factor4 = np.power(10, result0[0] - result4[0])

    dev3 = result0[3] + result3[3]
    dev4 = result0[3] + result4[3]

    dev3 = np.power(10, result0[0] - result3[0] + dev3) - factor3
    dev4 = np.power(10, result0[0] - result4[0] + dev4) - factor4                           
    
    output += "\n"
    output += "Slowdown P3" + " " + str(factor3) + "  " + str(dev3) + " \n"
    output += "Slowdown P4" + " " + str(factor4) + "  " + str(dev4) + " \n"

    f = open(standardFileName(SCRIPT_DIR) + "relativeRates.txt", 'w')
    f.write(output)
    f.close()    
    

# # The actual main method
if __name__ == '__main__':

    print(sys.argv)

    if len(sys.argv) > 1:

        toggle = str(sys.argv[1])
        trials = int(sys.argv[2])

        if toggle == "plots":     
            #doInference([1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6], trials)
             doInference([1e0, 1e-1, 1e-2], trials)

        if toggle == "slowDownStudy":
            doSlowdownStudy(trials)
