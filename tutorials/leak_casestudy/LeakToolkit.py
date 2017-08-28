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

from SeesawGate import SeesawRates
import numpy as np


ATIME_OUT = 0.1
# lets see the error bars I get here....
MINIMUM_FORWARD = 2
A_CONCENTRATION = 50e-9
INCREMENT_TRIALS = 5000
DNA = "DNA"


myMultistrand.setNumOfThreads(8)
myMultistrand.setTerminationCriteria(MINIMUM_FORWARD)
myMultistrand.setLeakMode()

outputFile = "data.txt"


def setMinimumSuccess(n):
    myMultistrand.setTerminationCriteria(n)


def setMaxTrials(n):
    myMultistrand.settings.max_trials = n


def setOutputFile(title, info="No Additional Information Provied"):
    date_str = datetime.now().strftime(' %d-%m+%H:%M')
    outputFile = title + date_str + ".txt"
    f = open(outputFile, 'w')
    f.write("File Created using Multistrand 2.1.\n" + info)
    f.close()

# it is intended that rate is a SeesawRate object


def writeRate(title, rate):
    date_str = datetime.now().strftime(' %d-%m+%H:%M')
    outputFile = title + date_str
    f = open(outputFile, 'a')
    f.write("Rate " + title + "\n")
    f.write(rate.str() + "\n\n")
    f.close()


def getOptions(trials, material, complexes,
               success_stop_conditions, failed_stop_conditions, T=25, supersample=25):

    o = Options(simulation_mode="First Step", substrate_type=material,
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


def calculateGateInputRate(gate_complex, input_complex, output_complex, trials=INCREMENT_TRIALS, title="Gate w/ Input", material="DNA",  alt_output_complex=None):
    success_stop_condition = StopCondition(
        Options.STR_SUCCESS, [(output_complex, Options.dissocMacrostate, 0)])
    failed_stop_condition = StopCondition(
        Options.STR_FAILURE, [(input_complex, Options.dissocMacrostate, 0)])

    try:
        if alt_output_complex is None:
            raise TypeError
        alt_success_stop_condition = StopCondition(
            Options.STR_ALT_SUCCESS, [(output_complex, Options.dissocMacrostate, 0)])
        myMultistrand.setOptionsFactory5(getOptions, trials, material,
                                         [gate_complex, input_complex],
                                         [success_stop_condition,
                                             alt_success_stop_condition],
                                         [failed_stop_condition])
    except TypeError:
        myMultistrand.setOptionsFactory5(getOptions, trials, material, [gate_complex, input_complex], [
                                         success_stop_condition], [failed_stop_condition])

    return getRates(title)


def calculateGateDissocRate(start_complex, output_complex, trials=INCREMENT_TRIALS, title="Gate Dissociation", material="DNA"):
    success_stop_condition = StopCondition(
        Options.STR_SUCCESS, [(output_complex, Options.dissocMacrostate, 0)])
    failed_stop_condition = StopCondition(
        Options.STR_FAILURE, [(start_complex, Options.dissocMacrostate, 0)])

    myMultistrand.setOptionsFactory5(getOptions, trials, material, [start_complex], [
        success_stop_condition], [failed_stop_condition])

    return getRates(title)

# MS: Adapted to return two rates within the rates object returned:
#     Both the rate for one complex being released but also the other,
#     as an alternative success condition for first step mode


def calculateGateGateLeak(gateA, gateB, trials=INCREMENT_TRIALS, title="Gate Gate Fuel", material="DNA"):
    # define stop conditions
    gateA_complex = gateA.gate_output_complex
    gateB_complex = gateB.gate_output_complex
    leak_complex = gateB.output_complex
    alt_leak_complex = gateA.output_complex

    success_stop_condition = StopCo tndition(
        Options.STR_SUCCESS, [(leak_complex,
                               Options.dissocMacrostate, 0)])

    alt_success_stop_condition = StopCondition(
        Options.STR_ALT_SUCCESS, [(alt_leak_complex, Options.dissocMacrostate, 0)])

    failed_stop_condition = StopCondition(
        Options.STR_FAILURE, [(gateA_complex,
                               Options.dissocMacrostate, 0)])

    myMultistrand.setOptionsFactory5(getOptions, trials, material,
                                     [gateA_complex, gateB_complex],
                                     [success_stop_condition,
                                         alt_success_stop_condition],
                                     [failed_stop_condition])

    return getRates(title)


def calculateBaseOutputRate(gate, trials=INCREMENT_TRIALS):
    title = "Forward Output Reaction"
    rates = calculateGateInputRate(
        gate.gate_output_complex, gate.input_complex, gate.output_complex, trials, title)
    return rates


def calculateReverseOutputRate(gate, trials=INCREMENT_TRIALS):
    title = "Reverse Output Reaction"
    rates = calculateGateInputRate(
        gate.gate_input_complex, gate.output_complex, gate.input_complex, trials, title)
    return rates


def calculateBaseFuelRate(gate, trials=INCREMENT_TRIALS,
                          ):
    title = "Forward Fuel Reaction"
    rates = calculateGateInputRate(
        gate.gate_input_complex, gate.fuel_complex, gate.input_complex, trials, title)
    return rates


def calculateReverseFuelRate(gate, trials=INCREMENT_TRIALS,
                             ):
    title = "Reverse Fuel Reaction"
    rates = calculateGateInputRate(
        gate.gate_fuel_complex, gate.input_complex, gate.fuel_complex, trials, title)
    return rates


def calculateBaseThresholdRate(gate, trials=INCREMENT_TRIALS):
    title = "Thresholding"
    rates = calculateGateInputRate(
        gate.threshold_complex, gate.input_complex, gate.threshold_free_waste_complex, trials, title)
    return rates

# MS: Calculates the leak between a gate and its fuel


def calculateGateFuelLeak(gate, trials=INCREMENT_TRIALS, material="DNA"):
    title = "Gate Fuel Leak"
    rates = calculateGateInputRate(
        gate.gate_output_complex, gate.fuel_complex, gate.output_complex, trials, title)
    return rates


def calculateOutputThresholdOcclusion(gate, trials=INCREMENT_TRIALS, material="DNA"):
    title = "Output Threshold Occlusion"
    rates = calculateGateInputRate(
        gate.threshold_complex, gate.output_complex, gate.threshold_complex_output_occluded, trials, title)
    return rates


def calculateOutputThresholdOcclusionUnbind(gate, trials=INCREMENT_TRIALS, material="DNA"):
    title = "Output Threshold Occlusion Unbind"
    rates = calculateGateDissocRate(
        gate.threshold_complex_output_occluded, gate.output_complex, trials, title)
    return rates


def calculateOutputGateOcclusion(gate, trials=INCREMENT_TRIALS, material="DNA"):
    title = "Output Gate Occlusion"
    rates = calculateGateInputRate(
        gate.gate_output_complex, gate.output_complex, gate.gate_output_complex_output_occluded, trials, title)
    return rates


def calculateOutputGateOcclusionUnbind(gate, trials=INCREMENT_TRIALS, material="DNA"):
    title = "Output Gate Occlusion Unbind"
    rates = calculateGateDissocRate(
        gate.gate_output_complex_output_occluded, gate.output_complex, trials, title)
    return rates


def calculateInputGateOcclusion(gate, trials=INCREMENT_TRIALS, material="DNA"):
    title = "Input Gate Occlusion"
    rates = calculateGateInputRate(
        gate.gate_input_complex, gate.input_complex, gate.gate_input_complex_input_occluded, trials, title)
    return rates


def calculateInputGateOcclusionUnbind(gate, trials=INCREMENT_TRIALS, material="DNA"):
    title = "Input Gate Occlusion Unbind"
    rates = calculateGateDissocRate(
        gate.gate_input_complex, gate.input_complex, trials, title)
    return rates


def getRates(tag):
    myMultistrand.run()
    rates = SeesawRates(myMultistrand.results)
    writeRate(tag, rates)
    print rates
    return rates
