# Frits Dannenberg, Caltech, 2016.
# Mrinank Sharma, Caltech SURF, 2017.
# fdann@caltech.edu , ms2314@cam.ac.uk
#
# Leak toolkit for seesaw circuits. Used to calculate relevant
# leak and non-leak rates for evaluation of the 'relative leak'
# rate.
#
# MS: Model adapted to use the DNA23 Arrhenius constants
import sys
from os.path import expanduser
from datetime import datetime

dirs = ["~/workspace/multistrand", "~/workspace/multistrandPy",
        "~/multistrand", "~/multistrandPy"]
for x in dirs:
    i = expanduser(x)
    sys.path.append(i)

from multistrand.concurrent import myMultistrand, MergeSim, FirstStepRate, Bootstrap
from multistrand.objects import StopCondition
from multistrand.options import Options
from msArrhenius import setArrheniusConstantsDNA23

import numpy as np


ATIME_OUT = 1.0
# lets see the error bars I get here....
MINIMUM_FORWARD = 2
A_CONCENTRATION = 50e-9
INCREMENT_TRIALS = 50000
DNA = "DNA"
DEFAULT_TEMPERATURE = 35


myMultistrand.setNumOfThreads(8)
myMultistrand.setTerminationCriteria(MINIMUM_FORWARD)
myMultistrand.setBootstrap(True, 10000)



def setMinimumSuccess(n):
    myMultistrand.setTerminationCriteria(n)

def setMaxTrials(n):
    myMultistrand.settings.max_trials = n

def setOutputFile(name):
    myMultistrand.setOutputFile(name)


def getOptions(trials, material, complexes,
               success_stop_conditions, failed_stop_conditions, T=25, simulationModeIn="First Step", supersample=25):

    o = Options(simulation_mode=simulationModeIn, substrate_type=material,
                rate_method="Metropolis", num_simulations=trials,
                simulation_time=ATIME_OUT, temperature=T)

    for x in complexes:
        x.boltzmann_supersample = supersample
        x.boltzmann_count = trials
        x.boltzmann_sample = True

    conds = []
    for x in [success_stop_conditions, failed_stop_conditions]:
        try:
            # x is a list
            conds += x
        except TypeError:
            # x is a single input
            conds.append(x)

    o.start_state = complexes
    o.stop_conditions = conds
    # Using new parameters.
    setArrheniusConstantsDNA23(o)
    return o


def calculateGateInputRate(gate_complex, input_complex, output_complex, trials=INCREMENT_TRIALS, material="DNA",  alt_output_complex=None):
    myMultistrand.setLeakMode()
    success_stop_condition = StopCondition(
        Options.STR_SUCCESS, [(output_complex, Options.dissocMacrostate, 0)])
    failed_stop_condition = StopCondition(
        Options.STR_FAILURE, [(input_complex, Options.dissocMacrostate, 0)])

    try:
        if alt_output_complex is None:
            raise TypeError
        alt_success_stop_condition = StopCondition(
            Options.STR_ALT_SUCCESS, [(output_complex, Options.dissocMacrostate, 0)])
        myMultistrand.setOptionsFactory6(getOptions, trials, material,
                                         [gate_complex, input_complex],
                                         [success_stop_condition,
                                             alt_success_stop_condition],
                                         [failed_stop_condition], DEFAULT_TEMPERATURE)
    except TypeError:
        myMultistrand.setOptionsFactory6(getOptions, trials, material, [gate_complex, input_complex], [
                                         success_stop_condition], [failed_stop_condition], DEFAULT_TEMPERATURE)

    return getRates()


def calculateGateDissocRate(start_complex, output_complex, trials=MINIMUM_FORWARD, material="DNA"):
    myMultistrand.setPassageMode()
    success_stop_condition = StopCondition(
        Options.STR_SUCCESS, [(output_complex, Options.dissocMacrostate, 0)])
    # no failure stop condition here...
    print trials
    myMultistrand.setOptionsFactory7(getOptions, trials, material, [start_complex], [
        success_stop_condition], [], DEFAULT_TEMPERATURE, "First Passage Time")

    return getRates()

# MS: Adapted to return two rates within the rates object returned:
#     Both the rate for one complex being released but also the other,
#     as an alternative success condition for first step mode
def calculateGateGateLeak(gateA, gateB, trials=INCREMENT_TRIALS, material="DNA"):
    myMultistrand.setLeakMode()
    # define stop conditions
    gateA_complex = gateA.gate_output_complex
    gateB_complex = gateB.gate_output_complex
    leak_complex = gateB.output_complex
    alt_leak_complex = gateA.output_complex

    success_stop_condition = StopCondition(
        Options.STR_SUCCESS, [(leak_complex,
                               Options.dissocMacrostate, 0)])

    alt_success_stop_condition = StopCondition(
        Options.STR_ALT_SUCCESS, [(alt_leak_complex, Options.dissocMacrostate, 0)])

    failed_stop_condition = StopCondition(
        Options.STR_FAILURE, [(gateA_complex,
                               Options.dissocMacrostate, 0)])

    myMultistrand.setOptionsFactory6(getOptions, trials, material,
                                     [gateA_complex, gateB_complex],
                                     [success_stop_condition,
                                         alt_success_stop_condition],
                                     [failed_stop_condition], DEFAULT_TEMPERATURE)

    return getRates()


def calculateBaseOutputRate(gate, trials=INCREMENT_TRIALS):
    myMultistrand.setExperimentTag("Forward Output Reaction")
    rates = calculateGateInputRate(
        gate.gate_output_complex, gate.input_complex, gate.output_complex, trials)
    return rates


def calculateReverseOutputRate(gate, trials=INCREMENT_TRIALS):
    myMultistrand.setExperimentTag("Reverse Output Reaction")
    rates = calculateGateInputRate(
        gate.gate_input_complex, gate.output_complex, gate.input_complex, trials)
    return rates


def calculateBaseFuelRate(gate, trials=INCREMENT_TRIALS):
    myMultistrand.setExperimentTag("Forward Fuel Reaction")
    rates = calculateGateInputRate(
        gate.gate_input_complex, gate.fuel_complex, gate.input_complex, trials)
    return rates


def calculateReverseFuelRate(gate, trials=INCREMENT_TRIALS,
                             ):
    myMultistrand.setExperimentTag("Reverse Fuel Reaction")
    rates = calculateGateInputRate(
        gate.gate_fuel_complex, gate.input_complex, gate.fuel_complex, trials)
    return rates


def calculateBaseThresholdRate(gate, trials=INCREMENT_TRIALS):
    myMultistrand.setExperimentTag("Thresholding")
    rates = calculateGateInputRate(
        gate.threshold_complex, gate.input_complex, gate.threshold_free_waste_complex, trials)
    return rates


def calculateGateFuelLeak(gate, trials=INCREMENT_TRIALS, material="DNA"):
    myMultistrand.setExperimentTag("Gate Fuel Leak")
    rates = calculateGateInputRate(
        gate.gate_output_complex, gate.fuel_complex, gate.output_complex, trials)
    return rates


def calculateOutputThresholdOcclusion(gate, trials=INCREMENT_TRIALS):
    myMultistrand.setExperimentTag("Output Threshold Occlusion")
    rates = calculateGateInputRate(
        gate.threshold_complex, gate.output_complex, gate.threshold_complex_output_occluded, trials)
    return rates


def calculateOutputThresholdOcclusionUnbind(gate, trials=MINIMUM_FORWARD):
    myMultistrand.setExperimentTag("Output Threshold Occlusion Unbind")
    rates = calculateGateDissocRate(
        gate.threshold_complex_output_occluded, gate.output_complex, trials)
    return rates


def calculateOutputGateOcclusion(gate, trials=INCREMENT_TRIALS):
    myMultistrand.setExperimentTag("Output Gate Occlusion")
    rates = calculateGateInputRate(
        gate.gate_output_complex, gate.output_complex, gate.gate_output_complex_output_occluded, trials)
    return rates


def calculateOutputGateOcclusionUnbind(gate, trials=MINIMUM_FORWARD):
    myMultistrand.setExperimentTag("Output Gate Occlusion Unbind")
    rates = calculateGateDissocRate(
        gate.gate_output_complex_output_occluded, gate.output_complex, trials)
    return rates


def calculateInputGateOcclusion(gate, trials=INCREMENT_TRIALS):
    myMultistrand.setExperimentTag("Input Gate Occlusion")
    rates = calculateGateInputRate(
        gate.gate_input_complex, gate.input_complex, gate.gate_input_complex_input_occluded, trials)
    return rates


def calculateInputGateOcclusionUnbind(gate, trials=MINIMUM_FORWARD):
    myMultistrand.setExperimentTag("Input Gate Occlusion Unbind")
    rates = calculateGateDissocRate(
        gate.gate_input_complex_input_occluded, gate.input_complex, trials)
    return rates


def getRates():
    myMultistrand.run()
    return myMultistrand.results
