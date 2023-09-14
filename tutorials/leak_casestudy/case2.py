# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
Clamped Seesaw Gate Case Study
"""

import time
from enum import Enum

from multistrand.concurrent import MergeSim
from multistrand.options import Options, Literals
from multistrand.experiment import \
    ClampedSeesawGate, seesaw_gate_fuel_catalysis, seesaw_gate_gate_leak, \
    seesaw_gate_output_production, seesaw_gate_fuel_leak, standardOptions

myMultistrand = MergeSim()

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
        print(f"Time since start: {elap_time} s\n")
    elif(elap_time < 6000):
        print(f"Time since start: {elap_time / 60} min\n")
    else:
        print(f"Time since start: {elap_time / 3600} hr\n")


def getExperiment(selIn):
    exp = selIn
    fileName = ""
    
    if exp == Experiment.GATE_OUTPUT_PRODUCTION:
        print("Gate Output Production")
        experiment = seesaw_gate_output_production
    elif exp == Experiment.GATE_FUEL_REGEN:
        print("Fuel Input Regeneration")
        experiment = seesaw_gate_fuel_catalysis
    elif exp == Experiment.GATE_FUEL_LEAK:
        print("Gate Fuel Leak")
        experiment = seesaw_gate_fuel_leak
    elif exp == Experiment.GATE_GATE_LEAK:
        print("Gate Gate Leak")
        experiment = seesaw_gate_gate_leak

    return experiment

# trials has to be the first argument.


def genOptions(trialsIn, gateA, sel, supersample=25, gateB=None):
    
    stdOptions = standardOptions(
        Literals.first_step, tempIn=45.0, trials=trialsIn, timeOut=ATIME_OUT)
    if gateB == None:
        getExperiment(sel)(stdOptions, gateA, trialsIn, supersample)
    else:
        getExperiment(sel)(stdOptions, gateA, gateB, trialsIn, supersample)
    stdOptions.DNA23Metropolis()
    return stdOptions


def runExperiment(trialsIn, gateA, sel, gateB=None, supersample=25):
    myMultistrand.setOptionsFactory5(
        genOptions, trialsIn, gateA, sel, supersample, gateB)
    myMultistrand.setFirstStepMode()
    myMultistrand.run()
    

# Trials in will determine the increment of the extra number of trials ran each time
#  (here, we keep running until we have a minimum number of succesful trials)

def runSimulations(trialsIn=1000):
    myMultistrand.setOutputFile("case2")
    
    # Here we say we are going to use 2 threads, storing only succesful data. We will require at least 2 succesful trials with a maximum number of trials of 2500000
    setupSimulationOptions(2, True, 2, 2.5e6, 1000)

    # Here we create two clamped seesaw gates, according to the defined interface.
    gateA = ClampedSeesawGate(*CL_LONG_GATE_A_SEQ)
    gateB = ClampedSeesawGate(*CL_LONG_GATE_B_SEQ)

    # In order to run different reactions, use a different enum for the experiment type!
    runExperiment(trialsIn, gateA, Experiment.GATE_GATE_LEAK,  gateB = gateB)


if __name__ == "__main__":
    runSimulations(1000)
