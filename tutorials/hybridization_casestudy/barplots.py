# Frits Dannenberg, Caltech, 2016.
# fdann@caltech.edu

# # this is the new large figure for the MS 2.0 paper.
from multistrand.experiment import standardOptions, hairpinclosing
from multistrand.objects import StopCondition, Complex, Domain, Strand
from multistrand.options import Options
from multistrand.utils import standardFileName
from multistrand.concurrent import myMultistrand

import matplotlib.pylab as plt

A_TIME_OUT = 3.0

HAIRPIN_STEM = "ACGTACGT"
HAIRPIN_LOOP = "TTTTT"
FLAMM_SEQ = "GGGATTTCTCGCTATTCCAGTGGGA"
YURK_T6E2003 = "ACTAATCCTCAGATCCAGCTAGTGTCCGTACT"

myMultistrand.setNumOfThreads(8) 

enum_hairpin = "Hairpin closing"
enum_flamm = "Bistable conformation - Flamm et al."  # figure 8
enum_yurke =  "Threeway strand displacement - Yurke and Mills" # Yurke and Mills -- T6 in table 1
    
nTrialsMod = 30  # number of trials per process


class settings(object):
    
    def __init__(self, title, seq, nTrials = 10):
        
        self.title = title
        self.seq = seq
        self.nTrials = nTrials


    def __str__(self):
        return self.seq

setting_hairpin = settings(enum_hairpin, HAIRPIN_STEM, myMultistrand.numOfThreads*nTrialsMod)
setting_flamm = settings(enum_flamm, FLAMM_SEQ, myMultistrand.numOfThreads*nTrialsMod)
settings_yurke = settings(enum_yurke, YURK_T6E2003, myMultistrand.numOfThreads*nTrialsMod)



def simulationHairpin(trialsIn):
    
    stdOptions = standardOptions(simMode=Options.trajectory, trials=trialsIn)
    stdOptions.simulation_time = A_TIME_OUT
    hairpinclosing(stdOptions, HAIRPIN_STEM, HAIRPIN_LOOP)
    
    return stdOptions



# Figure 8 of Flamm 2000 -- bistable. Compute transition time S0 -> S1
def simulationFlamm2000(trialsIn):
    
#     seq =     "GGCCCCTTTGGGGGCCAGACCCCTAAAGGGGTC"
#     struct0 = "((((((((((((((.....))))))))))))))" 
#     struct1 = "((((((....)))))).((((((....))))))"

#     seq =     "GGGATTTCTCGCTATTCCAGTGGGA"
    struct0 = "......(((((((.....)))))))"
    struct1 = "(((....))).....(((....)))"
        
    stdOptions = standardOptions(simMode=Options.trajectory, trials=trialsIn, tempIn = 37.0)
    stdOptions.substrate_type = Options.substrateRNA
    stdOptions.gt_enable = 1
    stdOptions.simulation_time = A_TIME_OUT
    
 
    stemdomain1 = Domain(name="stemdomain1", sequence=   FLAMM_SEQ)
    strand = Strand(name="top", domains=[stemdomain1])
    
    startComplex = Complex(strands=[strand], structure=  struct1)
    successComplex = Complex(strands=[strand], structure=struct0)

    # Stop when the exact full duplex is achieved.
    stopSuccess = StopCondition(Options.STR_SUCCESS, [(successComplex, Options.exactMacrostate, 0)])
    
    stdOptions.start_state = [startComplex]
    stdOptions.stop_conditions = [stopSuccess]
    
    
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


def computeHittingTimes(settings):
    
    if settings.title == enum_hairpin:
        myMultistrand.setOptionsFactory1(simulationHairpin, settings.nTrials)
            
    if settings.title == enum_flamm:
        myMultistrand.setOptionsFactory1(simulationFlamm2000, settings.nTrials)
    
    
    if settings.title == enum_yurke:
        myMultistrand.setOptionsFactory1(simulationYurke, settings.nTrials)
    
    if settings.title == enum_hairpin or settings.title == enum_flamm:
        myMultistrand.setPassageMode()   # non-first stepping mode, no need to store trajectory information
    
    if settings.title == enum_yurke:
        myMultistrand.setTrajectoryMode()   # retain all information you would normally find in a SimSystem.results object
        
    myMultistrand.run()
    
    
    return myMultistrand.results



def doBarplot(results, setting):
     
    fname = standardFileName("barplots", setting.title, "", setting.nTrials)
      
    def makeFig(title, times):
         
        fig = plt.figure()
        ax = fig.gca()
 
        ax.hist(times, 50, alpha=0.75, log = 1)
        ax.set_title("First passage time for " + title  )      
          
        ax = plt.gca()
        ax.set_ylabel('Trajectory counts (total = ' +  str(len(times)) + ')')  
        ax.set_xlabel('Trajectory time')
         
        plt.xticks(rotation=-40)
                 
        plt.tight_layout()
        plt.savefig(fname + "-bar" + "-" + title + '.pdf')
        plt.close()
         
         
    if setting.title == enum_hairpin or setting.title == enum_flamm:
        
        makeFig(setting.seq, results.times)

    if setting.title == enum_yurke:

        times = [ i.time for i in results.dataset if i.tag == Options.STR_SUCCESS]
        makeFig(setting.seq, times)
        
        fname = fname + "-failed"
        times = [ i.time for i in results.dataset if i.tag == Options.STR_FAILURE]
        makeFig(setting.seq, times)
        
    


def makePlots(setting):

    results = computeHittingTimes(setting)
    doBarplot(results, setting)



# The actual main method
if __name__ == '__main__':

    makePlots(setting_hairpin)
    makePlots(setting_flamm)
    makePlots(settings_yurke)
