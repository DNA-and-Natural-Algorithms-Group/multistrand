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

# HAIRPIN_STEM = "ACGTACGT"
# HAIRPIN_LOOP = "TTTTT"

HAIRPIN_STEM = "CCCAA"
HAIRPIN_LOOP = "T"*21


FLAMM_SEQ = "GGGATTTCTCGCTATTCCAGTGGGA"
YURK_T6E2003 = "ACTAATCCTCAGATCCAGCTAGTGTCCGTACT"

myMultistrand.setNumOfThreads(8) 


enum_bonnet = "bonnet"
enum_flamm = "flamm"
enum_yurke = "yurke"


title_bonnet = "Hairpin closing - Bonnet et al."
title_flamm = "Bistable conformation - Flamm et al."  # figure 8
title_yurke =  "Threeway strand displacement - Yurke and Mills" # Yurke and Mills -- T6 in table 1

title_hairpin_alt = "Hairpin opening - Bonnet et al."
# enum_flamm_alt = "Bistable conformation - Flamm et al."  # figure 8
# enum_yurke_alt =  "Threeway strand displacement - Yurke and Mills" # Yurke and Mills -- T6 in table 1


nTrialsMod = 30  # number of trials per process


class settings(object):
    
    def __init__(self, enum, seq, title,  title_alt = None, reverseIn = False, nTrials = 10):
        
    
        self.type = enum
        self.title = title
        self.title_alt = title_alt
        self.seq = seq
        self.nTrials = nTrials

        self.reverse = reverseIn


    def __str__(self):
        return self.seq

setting_bonnet = settings(enum_bonnet, HAIRPIN_STEM, title_bonnet, title_hairpin_alt, True, myMultistrand.numOfThreads*nTrialsMod)
setting_flamm = settings(enum_flamm, FLAMM_SEQ, title_flamm, nTrials = myMultistrand.numOfThreads*nTrialsMod)
settings_yurke = settings(enum_yurke, YURK_T6E2003, title_yurke, nTrials = myMultistrand.numOfThreads*nTrialsMod)



def simulationHairpin(trialsIn, reverse):
    
    stdOptions = standardOptions(simMode=Options.trajectory, trials=trialsIn)
    stdOptions.simulation_time = A_TIME_OUT
    stdOptions.temperature = 50.0
    
    
    if reverse:
        hairpinopening(stdOptions, HAIRPIN_STEM, HAIRPIN_LOOP)
    else:
        hairpinclosing(stdOptions, HAIRPIN_STEM, HAIRPIN_LOOP)
    
    
    return stdOptions



# Figure 8 of Flamm 2000 -- bistable. Compute transition time S0 -> S1
def simulationFlamm2000(trialsIn):
    
    seq =     "GGCCCCTTTGGGGGCCAGACCCCTAAAGGGGTC"
    
    structStart =   "................................."
    struct0 =       "((((((((((((((.....))))))))))))))" 
    struct1 =       "((((((....)))))).((((((....))))))"

# #     seq =     "GGGATTTCTCGCTATTCCAGTGGGA"
#     struct0 = "......(((((((.....)))))))"
#     struct1 = "(((....))).....(((....)))"
        
    stdOptions = standardOptions(simMode=Options.trajectory, trials=trialsIn, tempIn = 37.0)
    stdOptions.substrate_type = Options.substrateRNA
    stdOptions.gt_enable = 1
    stdOptions.simulation_time = A_TIME_OUT
    
 
    stemdomain1 = Domain(name="stemdomain1", sequence=   seq)
    strand = Strand(name="top", domains=[stemdomain1])
    
    startComplex = Complex(strands=[strand], structure=structStart)
    successComplex0 = Complex(strands=[strand], structure=struct0)
    successComplex1 = Complex(strands=[strand], structure=  struct1)

    # Stop when the exact full duplex is achieved.
    stopSuccess0 = StopCondition(Options.STR_SUCCESS, [(successComplex0, Options.exactMacrostate, 0)])
    stopSuccess1 = StopCondition(Options.STR_ALT_SUCCESS, [(successComplex1, Options.exactMacrostate, 0)])
    
    stdOptions.start_state = [startComplex]
    stdOptions.stop_conditions = [stopSuccess0, stopSuccess1]
    
    
    return stdOptions

    
## FD: not using multistrand.experiment.threewayDisplacement 
## because the toehold is on the 3' end
def simulationYurke(trialsIn):
    
    stdOptions = standardOptions(simMode=Options.firstPassageTime, trials=trialsIn)
    stdOptions.simulation_time = A_TIME_OUT

    domS    = Domain(sequence = "ACTAATCCTCAGATCCAGCTAGTGTC", name="d_S")
    domD    = Domain(sequence = "A", name="d_A")
    domT    = Domain(sequence = "CGTACT", name="d_T")
    
    
    strandQ = Strand(domains= [domS, domD])
    strandT = Strand(domains = [domT, domS])
    strandS = strandT.C

    complexStart = Complex(strands = [strandQ, strandS, strandT], structure="(.+)(+)." )
    complexEndS  = Complex(strands = [strandQ], structure="..")
    complexEndF  = Complex(strands = [strandT], structure="..")
    
    stopSuccess = StopCondition(Options.STR_SUCCESS, [(complexEndS, Options.dissocMacrostate, 3)])
    stopFailed = StopCondition(Options.STR_FAILURE, [(complexEndF, Options.dissocMacrostate, 3)])

    
    stdOptions.start_state = [complexStart]
    stdOptions.stop_conditions = [stopSuccess,stopFailed]    
    
    return stdOptions


def computeHittingTimes(settings, reverse = False):
    
    if settings.type == enum_bonnet:
        myMultistrand.setOptionsFactory2(simulationHairpin, settings.nTrials, reverse)
            
    if settings.type == enum_flamm:
        myMultistrand.setOptionsFactory1(simulationFlamm2000, settings.nTrials)
    
    
    if settings.type == enum_yurke:
        myMultistrand.setOptionsFactory1(simulationYurke, settings.nTrials)
    
    if settings.type == enum_bonnet or settings.title == enum_flamm:
        myMultistrand.setPassageMode()   # non-first stepping mode, no need to store trajectory information
    
    if settings.type == enum_yurke:
        myMultistrand.setTrajectoryMode()   # retain all information you would normally find in a SimSystem.results object
        
    myMultistrand.run()
    
    
    return myMultistrand.results



def doBarplot(results, setting):
     
    fname = standardFileName("barplots", setting.title, "", setting.nTrials)
      
    def makeFig(title, times):
         
        fig = plt.figure()
        ax = fig.gca()
 
        ax.hist(times, 50, alpha=0.75, log = 1)
        ax.set_title(title  )      
          
        ax = plt.gca()
        ax.set_ylabel('Trajectory counts (total = ' +  str(len(times)) + ')')  
        ax.set_xlabel('Trajectory time')
         
        plt.xticks(rotation=-40)
                 
        plt.tight_layout()
        plt.savefig(fname + "-bar" + "-" + title + '.pdf')
        plt.close()
         
         
    if setting.type == enum_bonnet:
        
        makeFig(setting.title, results.times)

    if setting.type == enum_flamm:
        
        makeFig(setting.title, results.times)


    if setting.type == enum_yurke:

        times = [ i.time for i in results.dataset if i.tag == Options.STR_SUCCESS]
        makeFig(setting.title, times)
        
        fname = fname + "-failed"
        times = [ i.time for i in results.dataset if i.tag == Options.STR_FAILURE]
        makeFig(setting.title, times)
        
def doDoubleBarplot(results, results2, setting):
     
    fname = standardFileName("barplots", setting.title, "", setting.nTrials)
 
    def computeCompletionLine(results):
        
        N= len(results)
        
        if N == 0:
            exit("Number of results in zero")
        
        results.sort()
        Y = (100.0/N) + (100.0/N) * np.array(range(N))

#         print str(Y)
        
        return results,Y
        
        
      
    def makeFig(title, times, times2):
         
        fig = plt.figure()
        ax = fig.gca()
 
        ax.hist(times, 25, alpha=0.20, log = 1, histtype='bar')
        ax.hist(times2, 25, alpha=0.20, log = 1, histtype='bar', stacked=True)
                
        survX, survY = computeCompletionLine(times)
        survX2, survY2 = computeCompletionLine(times2)

        ax.set_title( title  )      
        ax = plt.gca()
        ax.set_ylabel('Trajectory counts (total = ' +  str(len(times)) + ')')  
        ax.set_xlabel('Trajectory time')
        
        ax2 = ax.twinx()
        ax2.plot(survX, survY, lw=2)
        ax2.plot(survX2, survY2, lw=2)
        ax2.set_ylabel('Cummulative completion pct' +  str(len(times)) + ')')  

          
        
        plt.xticks(rotation=-30)
                 
        plt.tight_layout()
        plt.savefig(fname + "-bar" + "-" + title + '.pdf')
        plt.close()
         
         
    if setting.type == enum_bonnet:
        
        makeFig(setting.title, results.times, results2.times)

    if setting.type == enum_flamm:
        
        makeFig(setting.title, results.times)


    if setting.type == enum_yurke:

        times = [ i.time for i in results.dataset if i.tag == Options.STR_SUCCESS]
        makeFig(setting.title, times)
        
        fname = fname + "-failed"
        times = [ i.time for i in results.dataset if i.tag == Options.STR_FAILURE]
        makeFig(setting.title, times)


def makePlots(setting):

    results = computeHittingTimes(setting)
    
    if setting.reverse:
        results2 = computeHittingTimes(setting, True)
        doDoubleBarplot(results, results2, setting)
    
    else:
        doBarplot(results, setting)



# The actual main method
if __name__ == '__main__':

    makePlots(setting_bonnet)
    makePlots(setting_flamm)
    makePlots(settings_yurke)
