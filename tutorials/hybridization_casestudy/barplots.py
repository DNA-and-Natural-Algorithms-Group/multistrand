# Frits Dannenberg, Caltech, 2016.
# fdann@caltech.edu

# # this is the new large figure for the MS 2.0 paper.
from multistrand.experiment import standardOptions, hairpinclosing, threewayDisplacement
from multistrand.objects import StopCondition, Complex, Domain, Strand
from multistrand.options import Options
from multistrand.utils import standardFileName
from multistrand.concurrent import myMultistrand

import matplotlib.pylab as plt

A_TIME_OUT = 3.0

HAIRPIN_STEM = "ACGTACGT"
HAIRPIN_LOOP = "TTTTT"

FLAMM_SEQ = "GGGAUUUCUCGCUAUUCCAGUGGGA"

YURK_T6E2003 = "ACTAATCCTCAGATCCAGCTAGTGTCCGTACT"


myMultistrand.setNumOfThreads(8) 


enum_hairpin = "hairpin"
enum_flamm = "flamm2000"  # figure 8
enum_yurke =  "yurke2003" # Yurke and Mills -- T6 in table 1
    
nTrialsMod = 1


class settings(object):
    
    def __init__(self, title, seq, nTrials = 10):
        
        self.title = title
        self.seq = seq
        self.nTrials = nTrials


    def __str__(self):
        return self.seq

setting_hairpin = settings(enum_hairpin, HAIRPIN_STEM, 20* myMultistrand.numOfThreads*nTrialsMod)
setting_flamm = settings(enum_flamm, FLAMM_SEQ, 5 * myMultistrand.numOfThreads*nTrialsMod)
settings_yurke = settings(enum_yurke, YURK_T6E2003,  myMultistrand.numOfThreads*nTrialsMod)



def simulationHairpin(trialsIn):
    
    stdOptions = standardOptions(simMode=Options.trajectory, trials=trialsIn)
    stdOptions.simulation_time = A_TIME_OUT
    hairpinclosing(stdOptions, HAIRPIN_STEM, HAIRPIN_LOOP)
    
    return stdOptions



# Figure 8 of Flamm 2000 -- bistable. Compute transition time S0 -> S1
def simulationFlamm2000(trialsIn):
    
#     seq = "GGGAUUUCUCGCUAUUCCAGUGGGA" ## FD: morphing this to be T's
#     seq =     "GGCCCCTTTGGGGGCCAGACCCCTAAAGGGGTC"
#     struct0 = "((((((((((((((.....))))))))))))))" 
#     struct1 = "((((((....)))))).((((((....))))))"

    #          GGGAUUUCUCGCUAUUCCAGUGGGA
    seq =     "GGGATTTCTCGCTATTCCAGTGGGA"
    struct0 = "......(((((((.....)))))))"
    struct1 = "(((....))).....(((....)))"
        
    stdOptions = standardOptions(simMode=Options.trajectory, trials=trialsIn, tempIn = 37.0)
    stdOptions.substrate_type = Options.substrateRNA
    stdOptions.gt_enable = 1
    stdOptions.simulation_time = A_TIME_OUT
    
 
    stemdomain1 = Domain(name="stemdomain1", sequence=   seq)
    strand = Strand(name="top", domains=[stemdomain1])
    
    startComplex = Complex(strands=[strand], structure=  struct1)
    successComplex = Complex(strands=[strand], structure=struct0)


    # Stop when the exact full duplex is achieved.
    stopSuccess = StopCondition(Options.STR_SUCCESS, [(successComplex, Options.exactMacrostate, 0)])

    
    stdOptions.start_state = [startComplex]
    stdOptions.stop_conditions = [stopSuccess]
    
    
    return stdOptions
    

def simulationYurke(trialsIn):
    
    stdOptions = standardOptions(simMode=Options.trajectory, trials=trialsIn)
    stdOptions.simulation_time = A_TIME_OUT
    hairpinclosing(stdOptions, HAIRPIN_STEM, HAIRPIN_LOOP, myTrials=0)
    
    return stdOptions


def computeHittingTimes(settings):
    
    if settings.title == enum_hairpin:
        myMultistrand.setOptionsFactory1(simulationHairpin, settings.nTrials)
            
    if settings.title == enum_flamm:
        myMultistrand.setOptionsFactory1(simulationFlamm2000, settings.nTrials)

    myMultistrand.setPassageMode()
    myMultistrand.run()
    
    return myMultistrand.results.times



def doBarplot(times, setting):
     
    fname = standardFileName("barplots", setting.title, "", setting.nTrials)
      
    def makeFig(title, times):
         
        fig = plt.figure()
        ax = fig.gca()
 
        ax.hist(times, 50, alpha=0.75, log = 1)
        ax.set_title(("Trajectory times " + title + " count=" + str(len(times))))      
          
        ax = plt.gca()
        ax.set_ylabel('Trajectory counts')    
        ax.set_xlabel('Trajectory time')
         
        plt.xticks(rotation=-40)
                 
        plt.tight_layout()
        plt.savefig(fname + "-bar" + "-" + title + '.pdf')
        plt.close()
         
    makeFig(setting.seq, times)

    


def makePlots(setting):

    times = computeHittingTimes(setting)
    doBarplot(times, setting)



# The actual main method
if __name__ == '__main__':

    makePlots(setting_hairpin)
    makePlots(setting_flamm)
