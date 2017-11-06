# Frits Dannenberg, Caltech, 2016.
# fdann@caltech.edu

# # this is the new large figure for the MS 2.0 paper.
from multistrand.experiment import standardOptions, hairpinclosing, threewayDisplacement, makeComplex, setBoltzmann
from multistrand.objects import StopCondition
from multistrand.options import Options
from multistrand.utils import standardFileName
from multistrand.concurrent import myMultistrand

import matplotlib.pylab as plt
import enum

A_TIME_OUT = 4.0

HAIRPIN_STEM = "ACGTACGT"
HAIRPIN_LOOP = "TTTTT"

FLAMM_SEQ = "GGGAUUUCUCGCUAUUCCAGUGGGA"

YURK_T6E2003 = ""


myMultistrand.setNumOfThreads(8) 


enum_hairpin = "hairpin"
enum_flamm = "flamm2000"  # figure 8
enum_yurke =  "yurke2003" # Yurke and Mills -- T6 in table 1
    
nTrialsMod = 5


class settings(object):
    
    def __init__(self, title, seq, nTrials = 10):
        
        self.title = title
        self.seq = seq
        self.nTrials = nTrials


    def __str__(self):
        return self.seq

setting_hairpin = settings(enum_hairpin, HAIRPIN_STEM, 20* myMultistrand.numOfThreads*nTrialsMod)
setting_flamm = settings(enum_flamm, FLAMM_SEQ, 2 * myMultistrand.numOfThreads*nTrialsMod)



def simulationHairpin(trialsIn):
    
    stdOptions = standardOptions(simMode=Options.trajectory, trials=trialsIn)
    stdOptions.simulation_time = A_TIME_OUT
    hairpinclosing(stdOptions, HAIRPIN_STEM, HAIRPIN_LOOP, myTrials=0)
    
    return stdOptions
    

# Figure 8 of Flamm 2000 -- bistable. Compute transition time S0 -> S1
def simulationFlamm2000(trialsIn):
    
#     seq = "GGGAUUUCUCGCUAUUCCAGUGGGA" ## FD: morphing this to be T's
    seq =     "GGCCCCTTTGGGGGCCAGACCCCTAAAGGGGTC"
    struct0 = "((((((((((((((.....))))))))))))))" 
    struct1 = "((((((....)))))).((((((....))))))"
        
    stdOptions = standardOptions(simMode=Options.trajectory, trials=trialsIn, tempIn = 36.95)
    stdOptions.substrate_type = Options.substrateRNA
    stdOptions.simulation_time = A_TIME_OUT
    
    startComplex = makeComplex([seq], struct0)
    endComplex = makeComplex([seq], struct1)
    
    if(trialsIn > 0):
        setBoltzmann(startComplex, trialsIn)

    # stop when the invasion is complete, or when the invader dissociates
    stopSuccess = StopCondition(Options.STR_SUCCESS, [(endComplex, Options.countMacrostate , 3)])  # # three disagreements
    
    stdOptions.start_state = [startComplex]
    stdOptions.stop_conditions = [stopSuccess]
    
    
    return stdOptions
    




def computeHittingTimes(settings):
    
    if settings.title == enum_hairpin:
        myMultistrand.setOptionsFactory1(simulationHairpin, settings.nTrials)
        myMultistrand.setTrajectoryMode()
            
    if settings.title == enum_flamm:
        myMultistrand.setOptionsFactory1(simulationFlamm2000, settings.nTrials)
        myMultistrand.setTrajectoryMode()
        
            
    myMultistrand.run()
    
    return myMultistrand.results.dataset



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

    dataset = computeHittingTimes(setting)
    times = []
    
    for e in dataset:
        times.append(e.time)        

    print str(times)
    doBarplot(times, setting)




# The actual main method
if __name__ == '__main__':

    makePlots(setting_hairpin)
    makePlots(setting_flamm)
