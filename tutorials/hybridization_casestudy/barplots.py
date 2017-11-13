# Frits Dannenberg, Caltech, 2016.
# fdann@caltech.edu

# # this is the new large figure for the MS 2.0 paper.

from multistrand.experiment import standardOptions, hairpinclosing, hairpinopening
from multistrand.objects import StopCondition, Complex, Domain, Strand
from multistrand.options import Options
from multistrand.utils import standardFileName
from multistrand.concurrent import myMultistrand

import matplotlib.pylab as plt
import numpy as np

A_TIME_OUT = 600.0


HAIRPIN_STEM = "CCCAA"
HAIRPIN_LOOP = "T"*21

FLAMM_SEQ = "GGGATTTCTCGCTATTCCAGTGGGA"
YURK_T6E2003 = "ACTAATCCTCAGATCCAGCTAGTGTCCGTACT"

myMultistrand.setNumOfThreads(8) 


enum_bonnet = "bonnet"
enum_flamm = "flamm"
enum_yurke = "yurke"
enum_yurke2 = "yurke2" # this is to compute the hybridization rate of the toehold

title_bonnet = "Hairpin closing (blue) and opening (orange) - Bonnet et al."
title_flamm = "RNA kinetic trap - Flamm et al."  # figure 8
title_yurke = "Threeway strand displacement - Yurke and Mills"  # Yurke and Mills -- T6 in table 1
title_yurke2 = "Toehold binding rate - Yurke and Mills"  # Yurke and Mills -- T6 in table 1


nTrialsMod = 100  # number of trials per process


class settings(object):
    
    def __init__(self, enum, title, reverseIn=False, nTrials=10):
        
    
        self.type = enum
        self.title = title
        self.nTrials = nTrials

        self.reverse = reverseIn


    def __str__(self):
        return self.title

setting_bonnet = settings(enum_bonnet, title_bonnet, True, 10 * myMultistrand.numOfThreads * nTrialsMod)
setting_flamm = settings(enum_flamm, title_flamm, nTrials=10 * myMultistrand.numOfThreads * nTrialsMod)
settings_yurke = settings(enum_yurke, title_yurke, nTrials=myMultistrand.numOfThreads * nTrialsMod)
settings_yurke2 = settings(enum_yurke2, title_yurke2, nTrials=myMultistrand.numOfThreads * nTrialsMod)


def simulationHairpin(trialsIn, reverse):
    
    stdOptions = standardOptions(simMode=Options.trajectory, trials=trialsIn)
#     stdOptions.JSMetropolis25()
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

        
    stdOptions = standardOptions(simMode=Options.trajectory, trials=trialsIn, tempIn=37.0)
    stdOptions.substrate_type = Options.substrateRNA
    stdOptions.gt_enable = 1
    stdOptions.simulation_time = A_TIME_OUT
    
 
    stemdomain1 = Domain(name="stemdomain1", sequence=seq)
    strand = Strand(name="top", domains=[stemdomain1])
    
    startComplex = Complex(strands=[strand], structure=structStart)
    successComplex0 = Complex(strands=[strand], structure=struct0)
    successComplex1 = Complex(strands=[strand], structure=struct1)

    # Stop when the exact full duplex is achieved.
    stopSuccess0 = StopCondition(Options.STR_SUCCESS, [(successComplex0, Options.exactMacrostate, 0)])
    stopSuccess1 = StopCondition(Options.STR_ALT_SUCCESS, [(successComplex1, Options.exactMacrostate, 0)])
    
    stdOptions.start_state = [startComplex]
    stdOptions.stop_conditions = [stopSuccess0, stopSuccess1]
    
    
    return stdOptions

    
# # FD: not using multistrand.experiment.threewayDisplacement 
# # because the toehold is on the 3' end
def simulationYurke(trialsIn):
    
    stdOptions = standardOptions(simMode=Options.firstPassageTime, trials=trialsIn)
    stdOptions.simulation_time = A_TIME_OUT

    domS = Domain(sequence="ACTAATCCTCAGATCCAGCTAGTGTC", name="d_S")
    domD = Domain(sequence="A", name="d_A")
    domT = Domain(sequence="CGTACT", name="d_T")
    
    
    strandQ = Strand(domains=[domS, domD])
    strandT = Strand(domains=[domT, domS])
    strandS = strandT.C

    complexStart = Complex(strands=[strandQ, strandS, strandT], structure="(.+)(+).")
    complexEndS = Complex(strands=[strandQ], structure="..")
    complexEndF = Complex(strands=[strandT], structure="..")
    
    stopSuccess = StopCondition(Options.STR_SUCCESS, [(complexEndS, Options.dissocMacrostate, 3)])
    stopFailed = StopCondition(Options.STR_ALT_SUCCESS, [(complexEndF, Options.dissocMacrostate, 3)])
    
    stdOptions.start_state = [complexStart]
    stdOptions.stop_conditions = [stopSuccess, stopFailed]    
    
    return stdOptions


# # FD: not using multistrand.experiment.threewayDisplacement 
# # because the toehold is on the 3' end
def simulationYurke2(trialsIn):
    
    stdOptions = standardOptions(simMode=Options.firstPassageTime, trials=trialsIn)
    stdOptions.simulation_time = A_TIME_OUT

    domS = Domain(sequence="ACTAATCCTCAGATCCAGCTAGTGTC", name="d_S")
    domD = Domain(sequence="A", name="d_A")
    domT = Domain(sequence="CGTACT", name="d_T")
    
    
    strandQ = Strand(domains=[domS, domD])
    strandT = Strand(domains=[domT, domS])
    strandS = strandT.C

    complexEndS = Complex(strands=[strandQ], structure="..")
    complexEndF = Complex(strands=[strandT], structure="..")
    complexEndFC = Complex(strands=[strandQ, strandS], structure="(.+).")
    
    stopSuccess = StopCondition(Options.STR_SUCCESS, [(complexEndS, Options.dissocMacrostate, 3)])
    
    stdOptions.start_state = [complexEndF, complexEndFC]
    stdOptions.stop_conditions = [stopSuccess]
    
    stdOptions.join_concentration = 0.0001 # 100 microMolar    
    
    return stdOptions



def computeHittingTimes(settings, reverse=False):
    
    if settings.type == enum_yurke2:
        myMultistrand.setOptionsFactory1(simulationYurke2, settings.nTrials)

    
    if settings.type == enum_bonnet:
        myMultistrand.setOptionsFactory2(simulationHairpin, settings.nTrials, reverse)
            
    if settings.type == enum_flamm:
        myMultistrand.setOptionsFactory1(simulationFlamm2000, settings.nTrials)
    
    if settings.type == enum_yurke:
        myMultistrand.setOptionsFactory1(simulationYurke, settings.nTrials)

    
    if settings.type == enum_bonnet or settings.type == enum_yurke2:
        myMultistrand.setPassageMode()  # using the pre-set success / fail

    if settings.type == enum_flamm or settings.type == enum_yurke:
        myMultistrand.setTrajectoryMode()  # non-first stepping mode, no need to store trajectory information
    

        
    myMultistrand.run()
    
    return myMultistrand.results


def doBarplot(times, settings):
     
    fname = standardFileName("barplots", settings.type, "", settings.nTrials)
      
    fig = plt.figure()
    ax = fig.gca()
    
    ax.hist(times, 50, alpha=0.75, log=1)
    ax.set_title(settings.title)      
      
    ax = plt.gca()
    ax.set_ylabel('Trajectory counts (total = ' + str(len(times)) + ')')  
    ax.set_xlabel('Trajectory time')
     
    plt.xticks(rotation=-40)
             
    plt.tight_layout()
    plt.savefig(fname + "-bar" + "-" + settings.title + '.pdf')
    plt.close()
         
         
        
def doDoubleBarplot(times, times2, setting):
     
    fname = standardFileName("barplots", setting.type, "", setting.nTrials)
 
    def computeCompletionLine(results):
        
        N = len(results)
        
        if N == 0:
            exit("Number of results in zero")
        
        results.sort()
        Y = (100.0 / N) + (100.0 / N) * np.array(range(N))

#         print str(Y)
        
        return results, Y
        
    fig = plt.figure()
    ax = fig.gca()
    
    ax.hist(times, 25, alpha=0.20, log=1, histtype='bar')
    ax.hist(times2, 25, alpha=0.20, log=1, histtype='bar', stacked=True)
            
    survX, survY = computeCompletionLine(times)
    survX2, survY2 = computeCompletionLine(times2)
    
    ax.set_title(setting.title)      
    ax = plt.gca()
    
    if setting.type == enum_flamm:
            ax.set_ylabel('Trajectory counts (T1= ' + str(len(times)) + ', T2=)' + str(len(times)) + ")")  
    else:
        ax.set_ylabel('Trajectory counts (total = ' + str(len(times)) + ')')  
    
    ax.set_xlabel('Trajectory time')
    
    ax2 = ax.twinx()
    ax2.plot(survX, survY, lw=2)
    ax2.plot(survX2, survY2, lw=2)
    ax2.set_ylabel('Cummulative completion pct')  
    
      
    
    plt.xticks(rotation=-30)
             
    plt.tight_layout()
    plt.savefig(fname + "-bar" + "-" + setting.title + '.pdf')
    plt.close()
         


def makePlots(settings):

    results = computeHittingTimes(settings)
    
    if settings.type == enum_yurke2 :
        
        doBarplot(results.times, settings)
    
    
    if settings.type == enum_bonnet :
        
        results2 = computeHittingTimes(settings, True)
        doDoubleBarplot(results.times, results2.times, settings)
        
    if settings.type == enum_flamm or settings.type == enum_yurke:
        
        times = [i.time for i in results.dataset if i.tag == Options.STR_SUCCESS]       
        times2 = [i.time for i in results.dataset if i.tag == Options.STR_ALT_SUCCESS]
                
        doDoubleBarplot(times, times2, settings)
    


# The actual main method
if __name__ == '__main__':

    makePlots(setting_bonnet)
    makePlots(setting_flamm)
    makePlots(settings_yurke)
    makePlots(settings_yurke2)

