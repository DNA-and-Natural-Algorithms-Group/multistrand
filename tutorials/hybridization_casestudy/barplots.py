# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
This is the new large figure for the MS 2.0 paper.
"""

import sys

from multistrand.experiment import standardOptions, hairpinclosing, hairpinopening
from multistrand.objects import StopCondition, Complex, Domain, Strand
from multistrand.options import Literals
from multistrand.utils import standardFileName, printTrajectory, dGC_feature
from multistrand.concurrent import MergeSim
from multistrand.system import SimSystem

import matplotlib
matplotlib.use('agg')

import matplotlib.pylab as plt
import numpy as np

A_TIME_OUT = 2000.0  # 10 s timeout
NUM_PROCESS = 8
nTrialsMod = 4  # number of trials per process

HAIRPIN_STEM = "CCCAA"
HAIRPIN_LOOP = "T"*21

FLAMM_SEQ = "GGGATTTCTCGCTATTCCAGTGGGA"
YURK_T6E2003 = "ACTAATCCTCAGATCCAGCTAGTGTCCGTACT"

YURKE2_CONCENTRATION = 0.0001  # 100 microMolar    

FIGURE_SIZE = (6 * 0.93, 4 * 0.93)

enum_bonnet = "bonnet"
enum_flamm = "flamm"
enum_yurke = "yurke"
enum_yurke2 = "yurke2"  # this is to compute the hybridization rate of the toehold
enum_rickettsia = "rickettsia"

title_bonnet = "Hairpin closing and opening - Bonnet et al."
title_flamm = "RNA kinetic trap - Flamm et al."  # figure 8
title_yurke = "Threeway strand displacement - Yurke and Mills"  # Yurke and Mills -- T6 in table 1
title_yurke2 = "Toehold binding rate - Yurke and Mills"  # Yurke and Mills -- T6 in table 1
title_rickettsia = "An autonomous polymerization motor powered \n by DNA hybridization - Venkataraman et al. "  # An autonomous polymerization motor powered by DNA hybridization SUVIR VENKATARAMAN, ROBERT M. DIRKS, PAUL W. K. ROTHEMUND, ERIK WINFREE AND NILES A. PIERCE


class settings(object):
    
    def __init__(self, enum, title, reverseIn=False, nTrials=10):
    
        self.type = enum
        self.title = title
        self.nTrials = nTrials

        self.reverse = reverseIn

    def __str__(self):
        return self.title


def simulationHairpin(trialsIn, reverse):
    
    stdOptions = standardOptions(simMode=Literals.trajectory, trials=trialsIn)
#     stdOptions.JSDefault()
    stdOptions.DNA23Metropolis()
    stdOptions.simulation_time = A_TIME_OUT
    stdOptions.temperature = 50.0
    
    if reverse:
        hairpinopening(stdOptions, HAIRPIN_STEM, HAIRPIN_LOOP)
    else:
        hairpinclosing(stdOptions, HAIRPIN_STEM, HAIRPIN_LOOP)
    
    return stdOptions


# Figure 8 of Flamm 2000 -- bistable. Compute transition time S0 -> S1
def simulationFlamm2000(trialsIn):
    
    seq = "GGCCCCTTTGGGGGCCAGACCCCTAAAGGGGTC"
    
    structStart = "................................."
    struct0 = "((((((((((((((.....))))))))))))))" 
    struct1 = "((((((....)))))).((((((....))))))"
        
    stdOptions = standardOptions(simMode=Literals.trajectory, trials=trialsIn, tempIn=37.0)
    stdOptions.substrate_type = Literals.substrateRNA
    stdOptions.gt_enable = 1
    stdOptions.simulation_time = A_TIME_OUT
    stdOptions.DNA23Metropolis()
 
    stemdomain1 = Domain(name="stemdomain1", sequence=seq)
    strand = Strand(name="top", domains=[stemdomain1])
    
    startComplex = Complex(strands=[strand], structure=structStart)
    successComplex0 = Complex(strands=[strand], structure=struct0)
    successComplex1 = Complex(strands=[strand], structure=struct1)

    # Stop when the exact full duplex is achieved.
    stopSuccess0 = StopCondition(Literals.success, [(successComplex0, Literals.exact_macrostate, 0)])
    stopSuccess1 = StopCondition(Literals.alt_success, [(successComplex1, Literals.exact_macrostate, 0)])
    
    stdOptions.start_state = [startComplex]
    stdOptions.stop_conditions = [stopSuccess0, stopSuccess1]
    
    return stdOptions

    
# # FD: not using multistrand.experiment.threewayDisplacement 
# # because the toehold is on the 3' end
def simulationYurke(trialsIn):
    
    stdOptions = standardOptions(simMode=Literals.first_passage_time, trials=trialsIn)
    stdOptions.simulation_time = A_TIME_OUT
    stdOptions.DNA23Metropolis()
   
    stdOptions.temperature = 25.0

    domS = Domain(sequence="ACTAATCCTCAGATCCAGCTAGTGTC", name="d_S")
    domD = Domain(sequence="A", name="d_A")
    domT = Domain(sequence="CGTACT", name="d_T")
    
    strandQ = Strand(domains=[domS, domD])
    strandT = Strand(domains=[domT, domS])
    strandS = strandT.C

    complexStart = Complex(strands=[strandQ, strandS, strandT], structure="(.+)(+).")
    complexEndS = Complex(strands=[strandQ], structure="..")
    complexEndF = Complex(strands=[strandT], structure="..")  # # ALT_SUCCESS is dissociation
    
    stopSuccess = StopCondition(Literals.success, [(complexEndS, Literals.dissoc_macrostate, 3)])
    stopFailed = StopCondition(Literals.alt_success, [(complexEndF, Literals.dissoc_macrostate, 3)])
    
    stdOptions.start_state = [complexStart]
    stdOptions.stop_conditions = [stopSuccess, stopFailed]    
    
    return stdOptions


# # FD: not using multistrand.experiment.threewayDisplacement 
# # because the toehold is on the 3' end
def simulationYurke2(trialsIn):
    
    stdOptions = standardOptions(simMode=Literals.first_passage_time, trials=trialsIn)
    stdOptions.simulation_time = A_TIME_OUT
    stdOptions.DNA23Metropolis()

    domS = Domain(sequence="ACTAATCCTCAGATCCAGCTAGTGTC", name="d_S")
    domD = Domain(sequence="A", name="d_A")
    domT = Domain(sequence="CGTACT", name="d_T")
    
    strandQ = Strand(domains=[domS, domD])
    strandT = Strand(domains=[domT, domS])
    strandS = strandT.C

#     complexEndS = Complex(strands=[strandQ], structure="..")
    complexEndF = Complex(strands=[strandT], structure="..")
    complexEndFC = Complex(strands=[strandQ, strandS], structure="(.+).")
    
    complexAttached = Complex(strands=[strandQ, strandS, strandT], structure="**+*(+)*")
    
    stopSuccess = StopCondition(Literals.success, [(complexAttached, Literals.loose_macrostate, 1)])
    
    stdOptions.start_state = [complexEndF, complexEndFC]
    stdOptions.stop_conditions = [stopSuccess]
    
    stdOptions.join_concentration = 0.0001  # 100 microMolar    
    
    return stdOptions


'''
An autonomous polymerization motor
powered by DNA hybridization
SUVIR VENKATARAMAN,
ROBERT M. DIRKS,
PAUL W. K. ROTHEMUND,
ERIK WINFREE AND NILES A. PIERCE
'''

# domains a, b, c, x, y,         (lengths 6, 18, 6, 3, and 3)

# H1: a b c b* x*
# H1: ATTCAAGCGACACCGTGGACGTGCACCCACGCACGTCCACGGTGTCGCACC

# dom_a = ATTCAA
# dom b = GCGACACCGTGGACGTGC
# dom c = ACCCAC
# dom x = *(ACC) = GGT

# H2: y* b* a* b c* 
# H2: GTTGCACGTCCACGGTGTCGCTTGAATGCGACACCGTGGACGTGCGTGGGT

# dom y = c(GTT) = AAC

#  A: b* a*
#  A: GCACGTCCACGGTGTCGCTTGAAT

#  R: x b y
#  R: GGTGCGACACCGTGGACGTGCAAC 


def simulationRickettsia(trialsIn):
    
    stdOptions = standardOptions(simMode=Literals.first_passage_time, trials=trialsIn)
    stdOptions.simulation_time = A_TIME_OUT
    stdOptions.DNA23Metropolis()
    stdOptions.temperature = 25.0
    stdOptions.magnesium = 0.0125
    stdOptions.sodium = 0.1

    dom_a = Domain(sequence="ATTCAA", name="a")  # length 6
    dom_b = Domain(sequence="GCGACACCGTGGACGTGC", name="b")  # length 18
    dom_c = Domain(sequence="ACCCAC", name="c")  # length 6
    
    dom_x = Domain(sequence="GGT", name="x")  # length 3
    dom_y = Domain(sequence="AAC", name="y")  # length 3
        
    strand_H1 = Strand(domains=[dom_a, dom_b, dom_c, dom_b.C, dom_x.C])
    strand_H2 = Strand(domains=[dom_y.C, dom_b.C, dom_a.C, dom_b, dom_c.C])
    
    strand_A = Strand(domains=[dom_b.C, dom_a.C])
    strand_B = Strand(domains=[dom_x, dom_b, dom_y])
    strand_R = Strand(domains=[dom_x, dom_b, dom_y])

    H1 = Complex(strands=[strand_H1], structure=".(.).", name="H1")
    H2 = Complex(strands=[strand_H2], structure=".(.).", name="H2")
    
#     state1 = Complex(strands=[strand_H1, strand_R, strand_A], structure="((.)*+*(.+))")  # domain x does not have to be bound
    state2 = Complex(strands=[strand_H1, strand_R, strand_A], structure="((.((+)).+))", name="state2")
    state3 = Complex(strands=[strand_H1, strand_R, strand_H2, strand_A], structure="(((((+))(+)(.))+))")
#     state4 = Complex(strands=[strand_H1, strand_R, strand_H2, strand_A], structure="(((((+)((+)).))+))", name = "state4")
    state5 = Complex(strands=[strand_H1, strand_R, strand_H2, strand_A], structure="((((.+.((+)).))+))", name="state5")
#     state6 = Complex(strands=[strand_H1, strand_H1, strand_R, strand_H2, strand_A], structure="((((.+((.)*+*((+)))))+))")  # domain x does not have to be bound
#     state6 = Complex(strands=[strand_H1, strand_H1, strand_R, strand_H2, strand_A], structure="((((.+((.).+.((+)))))+))", name = "state6")

#     state7 = Complex(strands=[strand_H1, strand_H1, strand_R, strand_H2, strand_A], structure="((((.+((.((+))*+*))))+))", name = "state7")

    stopFailure = StopCondition(Literals.failure, [(state2, Literals.dissoc_macrostate, 0)])
    stopSuccess = StopCondition(Literals.success, [(state5, Literals.loose_macrostate, 6)])
    
    stdOptions.start_state = [state3]
    stdOptions.stop_conditions = [stopSuccess, stopFailure]
    
    stdOptions.join_concentration = 0.001
    
    return stdOptions


def computeHittingTimes(settings, reverse=False):
    
    myMultistrand = MergeSim()
    myMultistrand.setNumOfThreads(NUM_PROCESS)
    
    if settings.type == enum_yurke2:
        myMultistrand.setOptionsFactory1(simulationYurke2, settings.nTrials)
    
    if settings.type == enum_bonnet:
        myMultistrand.setOptionsFactory2(simulationHairpin, settings.nTrials, reverse)
            
    if settings.type == enum_flamm:
        myMultistrand.setOptionsFactory1(simulationFlamm2000, settings.nTrials)
    
    if settings.type == enum_yurke:
        myMultistrand.setOptionsFactory1(simulationYurke, settings.nTrials)
    
    if settings.type == enum_rickettsia:
        myMultistrand.setOptionsFactory1(simulationRickettsia, 12 * settings.nTrials)
        myMultistrand.setTerminationCriteria(terminationCount=settings.nTrials)
    
    if settings.type == enum_bonnet or settings.type == enum_yurke2:
        myMultistrand.setPassageMode()  # using the pre-set success / fail

    if settings.type == enum_flamm or settings.type == enum_yurke or settings.type == enum_rickettsia:
        # non-first stepping mode, no need to store trajectory information
        myMultistrand.setPassageMode()
        
    myMultistrand.run()
    
    return myMultistrand.results

 
def computeCompletionLine(results, N=None):
    
    if N == None:
        N = len(results)
    
    if N == 0:
        exit("Number of results in zero")
    
    results.sort()
    Y = (100.0 / N) + (100.0 / N) * np.arange(len(results))
    
    return results, Y

    
def setLabelAndClose(settings, plt, ax):
    
    fname = standardFileName("barplots", settings.type, "", settings.nTrials)

    ax.set_xlabel('Trajectory time (ms)')
    
    plt.xticks(rotation=-30)
    plt.tight_layout()
    plt.savefig(fname + "-bar" + "-" + settings.title + '.pdf')
    plt.close()
         

def removeOutliers(times):
    
    times = np.array(sorted(times))
    return times[np.arange(int(len(times) * 0.98))]
    

def doBarplot(times, settings):

    observations = str(len(times))
    print("Number of observations is " + str(observations))
    
    times = [1000 * ele for ele in times]
      
    fig = plt.figure(figsize=FIGURE_SIZE)
    ax = fig.gca()

    myMin = min(times)
    myMax = max(removeOutliers(times))
    # compute the max for outliers removed.
    binwidth = (myMax - myMin) / 40 
    
    myBins = np.arange(myMin, myMax + binwidth, binwidth)
    
    weights1 = np.empty_like(times)
    weights1.fill(100.0 / len(times))
    
    ax.hist(times, alpha=0.20, log=1, bins=myBins, weights=weights1)
    ax.set_title(settings.title)      
      
    ax = plt.gca()
    ax.set_ylabel('Trajectory % (total = ' + observations + ')')  
    ax.set_ylim([0.1, 40.0])
    ax.set_xlim([0.0, myMax])
    
    survX, survY = computeCompletionLine(times)
             
    ax2 = ax.twinx()
    ax2.plot(survX, survY, lw=2)
    ax2.set_ylabel('Cummulative completion %')  
    ax.set_xlim([0.0, myMax])

    setLabelAndClose(settings, plt, ax)
         
        
def doDoubleBarplot(times, times2, setting):
     
    observations = str(len(times))
    observations2 = str(len(times2))
     
    times = [1000 * ele for ele in times]
    times2 = [1000 * ele for ele in times2]     

    fig = plt.figure(figsize=FIGURE_SIZE)
    ax = fig.gca()
    
    myMin = min(min(times), min(times2))
    myMax = max(max(removeOutliers(times)), max(removeOutliers(times2)))
    # compute the max for outliers removed.
    binwidth = (myMax - myMin) / 40 
    
    myBins = np.arange(myMin, myMax + binwidth, binwidth)
    
    weights1 = np.empty_like(times)
    weights1.fill(100.0 / len(times))
    
    weights2 = np.empty_like(times2)
    weights2.fill(100.0 / len(times2))
    
    ax.hist(times, alpha=0.20, log=1, histtype='bar', bins=myBins, weights=weights1)
    ax.hist(times2, alpha=0.20, log=1, histtype='bar', bins=myBins, weights=weights2, stacked=True)
            
    N = len(times) + len(times2)
            
    survX, survY = computeCompletionLine(times, N)
    survX2, survY2 = computeCompletionLine(times2, N)
    
    ax.set_title(setting.title)      
    ax = plt.gca()
    ax.set_ylim([0.1, 40.0])
    ax.set_xlim([0, myMax])
    
    if setting.type == enum_flamm or setting.type == enum_yurke or setting.type == enum_rickettsia:
            ax.set_ylabel('Trajectory counts (' + observations + ' and ' + observations2 + ")")  
    else:
        ax.set_ylabel('Trajectory counts (total = ' + observations + ')')  
    
    ax2 = ax.twinx()
    ax2.plot(survX, survY, lw=2)
    ax2.plot(survX2, survY2, lw=2)
    ax2.set_ylim([0.0, 100.0])
    ax2.set_ylabel('Cummulative completion %')  
    ax2.set_xlim([0, myMax])
    
    setLabelAndClose(setting, plt, ax)


def makePlots(settings):

    results = computeHittingTimes(settings)
    times = [i.time for i in results.dataset if i.tag == Literals.success]
    
    if settings.type == enum_yurke2 :
        
        doBarplot(times, settings)
        print("rate toehold binding = " + str(settings.nTrials / (sum(times))))
    
    if settings.type == enum_bonnet :
        
        results2 = computeHittingTimes(settings, True)
        times2 = [i.time for i in results2.dataset if i.tag == Literals.success]
        
        doDoubleBarplot(times, times2, settings)
        
    if settings.type == enum_flamm or settings.type == enum_yurke:
        
        times = [i.time for i in results.dataset if i.tag == Literals.success]       
        times2 = [i.time for i in results.dataset if i.tag == Literals.alt_success]

        if not len(times) == 0 and not sum(times) == 0:        
            print("rate reaction 1 = " + str(len(times) / sum(times)))
        if not len(times2) == 0 and not sum(times2) == 0:
            print("rate reaction 2 = " + str(len(times2) / sum(times2)))

        doDoubleBarplot(times, times2, settings)

    if settings.type == enum_rickettsia:
            
        times = [i.time for i in results.dataset if i.tag == Literals.success]       
        times2 = [i.time for i in results.dataset if i.tag == Literals.failure]

        if not len(times) == 0 and not sum(times) == 0:        
            print("rate reaction 1 = " + str(len(times) / sum(times)))
        if not len(times2) == 0 and not sum(times2) == 0:
            print("rate reaction 2 = " + str(len(times2) / sum(times2)))
            
        if len(times) > 0:
            settings.type = enum_rickettsia + "-1"
            doBarplot(times, settings)
        if len(times2) > 0:
            settings.type = enum_rickettsia + "-2"
            doBarplot(times2, settings)
            
        settings.type = enum_rickettsia 


def debugTester():
    """" Debug tester. """
    options = simulationRickettsia(trialsIn=1)
    options.simulation_mode = Literals.trajectory
    options.output_interval = 10000
    options.temperature = 25.0
    options.simulation_time = 5.0e-1
    options.join_concentration = 0.1
    s = SimSystem(options)
    s.start()
    printTrajectory(options, timescale=(1e9, "ns"), feature=dGC_feature)


if __name__ == '__main__':
    if len(sys.argv) > 2:
        
        NUM_PROCESS = int(sys.argv[1])
        nTrials = int(sys.argv[2])
        type = sys.argv[3]

        # by default, the first two examples get 10x more trajectories
        # For rickettsia, nTrials is the number of succesful trajectoires.
        settings_bonnet = settings(enum_bonnet, title_bonnet, True, 5 * nTrials)
        settings_flamm = settings(enum_flamm, title_flamm, nTrials=5 * nTrials)
        settings_yurke = settings(enum_yurke, title_yurke, nTrials=nTrials)
        settings_yurke2 = settings(enum_yurke2, title_yurke2, nTrials=nTrials)
        settings_rickettsia = settings(enum_rickettsia, title_rickettsia, nTrials=0.2 * nTrials)
        
        switcher = {
            
            enum_bonnet : settings_bonnet,
            enum_flamm : settings_flamm,
            enum_yurke : settings_yurke,
            enum_yurke2 : settings_yurke2,
            enum_rickettsia : settings_rickettsia,
            
            }
        
#         debugTester()
        makePlots(switcher[type])
        
    else:
        
        print("Please supply the number of processes and total number of trajectories to simulate per case study, and the type \n")
        print("Example: python barplot.py 2 20 bonnet")
        print("Allowed <type> :     bonnet, flamm, yurke, yurke2, rickettsia")
