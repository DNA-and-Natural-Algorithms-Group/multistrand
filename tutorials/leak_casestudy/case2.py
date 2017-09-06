# Mrinank Sharma, Summer 2017, ms2314@cam.ac.uk
# Clamped Seesaw Gate Case Study
import time
from enum import Enum

from multistrand.concurrent import myMultistrand, MergeSim, FirstStepRate, Bootstrap
from multistrand.objects import StopCondition
from multistrand.options import Options
from multistrand.experiment import ClampedSeesawGate, seesaw_gate_fuel_catalysis, seesaw_gate_gate_leak, seesaw_gate_output_production, seesaw_gate_fuel_leak, standardOptions


class Experiment(Enum):
    GATE_OUTPUT_PRODUCTION = 0
    GATE_FUEL_REGEN = 1
    GATE_FUEL_LEAK = 2
    GATE_GATE_LEAK = 3


# Sequence design from Thubagere, 2017
CL_LONG_S18 = "TCTTCTAACAT"
CL_LONG_S5 = "CCACCAAACTT"
CL_LONG_S6 = "TAACACAATCA"
CL_LONG_S29 = "CCAATACTCCT"
CL_LONG_S53 = "TATCTAATCTC"
CL_LONG_S44 = "AAACTCTCTCT"
CL_LONG_SEQT = "TCT"
CLAMP_SEQ = "CA"


CL_LONG_GATE_A_SEQ = [CL_LONG_S44, CL_LONG_S18,
                      CL_LONG_S5, CL_LONG_S29, CL_LONG_SEQT, CLAMP_SEQ]
CL_LONG_GATE_B_SEQ = [CL_LONG_S53, CL_LONG_S5,
                      CL_LONG_S6, CL_LONG_S29, CL_LONG_SEQT, CLAMP_SEQ]


# ATIME_OUT **must** be a float, or an error is thrown !
ATIME_OUT = 1.0
DNA = "DNA"

start_time = time.time()


# Setup simulations options here (all at once - perhaps not so great)
def setupSimulationOptions(numThreads=8, leakMode=True, minimumSuccess=0, maxTrials=0, nBoot = 0):
    myMultistrand.setNumOfThreads(numThreads)
    if leakMode:
        # Might be running a lot of trials here - save memory
        myMultistrand.setLeakMode()
    else:
        pass

    if minimumSuccess != 0:
        myMultistrand.setTerminationCriteria(minimumSuccess)

    if maxTrials != 0:
        myMultistrand.settings.max_trials = maxTrials
    
    if nBoot != 0:
        myMultistrand.setBootstrap(True, nBoot)

    


# Useful Printer Function for how long it's been since we started (if running a batch of trials)
def printTimeElapsed():
    curr_time = time.time()
    elap_time = curr_time - start_time
    if(elap_time < 1000):
        print "Time since start: {} s\n".format(elap_time)
    elif(elap_time < 6000):
        print "Time since start: {} min\n".format(elap_time / 60)
    else:
        print "Time since start: {} hr\n".format(elap_time / 3600)


def getExperiment(selIn):
    exp = selIn
    if exp == Experiment.GATE_OUTPUT_PRODUCTION:
        myMultistrand.setExperimentTag("Gate Output Production")
        experiment = seesaw_gate_output_production
    elif exp == Experiment.GATE_FUEL_REGEN:
        myMultistrand.setExperimentTag("Fuel Input Regeneration")
        experiment = seesaw_gate_fuel_catalysis
    elif exp == Experiment.GATE_FUEL_LEAK:
        myMultistrand.setExperimentTag("Gate Fuel Leak")
        experiment = seesaw_gate_fuel_leak
    elif exp == Experiment.GATE_GATE_LEAK:
        myMultistrand.setExperimentTag("Gate Gate Leak")
        experiment = seesaw_gate_gate_leak

    return experiment

# trials has to be the first argument.


def genOptions(trialsIn, gateA, sel, supersample=25, gateB=None):
    stdOptions = standardOptions(
        Options.firstStep, tempIn=25.0, trials=trialsIn, timeOut=ATIME_OUT)
    if gateB == None:
        getExperiment(sel)(stdOptions, gateA, trialsIn, supersample)
    else:
        getExperiment(sel)(stdOptions, gateA, gateB, trialsIn, supersample)
    stdOptions.DNA23Metropolis()
    return stdOptions


def runExperiment(trialsIn, gateA, sel, gateB=None, supersample=25):
    myMultistrand.setOptionsFactory5(
        genOptions, trialsIn, gateA, sel, supersample, gateB)
    myMultistrand.run()
    

# Trials in will determine the increment of the extra number of trials ran each time
#  (here, we keep running until we have a minimum number of succesful trials)

def runSimulations(trialsIn=1000):
    # uncomment for logging
    # myMultistrand.setOutputFile("case2")
    
    # Here we say we are going to use 2 threads, storing only succesful data. We will require at least 2 succesful trials with a maximum number of trials of 2500000
    setupSimulationOptions(2, True, 2, 2.5e6, 1000)

    # Here we create two clamped seesaw gates, according to the defined interface.
    gateA = ClampedSeesawGate(*CL_LONG_GATE_A_SEQ)
    gateB = ClampedSeesawGate(*CL_LONG_GATE_B_SEQ)

    # In order to run different reactions, use a different enum for the experiment type!
    runExperiment(trialsIn, gateA, Experiment.GATE_OUTPUT_PRODUCTION)


if __name__ == "__main__":
    runSimulations(1000)
