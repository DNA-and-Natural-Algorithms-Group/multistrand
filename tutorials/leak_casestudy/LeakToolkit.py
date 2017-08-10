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


ATIME_OUT = 10.0
# lets see the error bars I get here....
MINIMUM_FORWARD = 2
A_CONCENTRATION = 50e-9
INCREMENT_TRIALS = 2000
DNA = "DNA"


myMultistrand.setNumOfThreads(2)
#myMultistrand.setTerminationCriteria(MINIMUM_FORWARD)
myMultistrand.setLeakMode()


def getOptions(trials, material, complex1, complex2,
               success_stop_conditions, failed_stop_conditions, T=25):

    o = Options(simulation_mode="First Step", substrate_type=material,
                rate_method="Metropolis", num_simulations=trials,
                simulation_time=ATIME_OUT, temperature=T)

    o.start_state = [complex1, complex2]
    conds = []

    for x in [success_stop_conditions, failed_stop_conditions]:
        try:
            # x is a list
            conds += x
        except TypeError:
            # x is a single input
            conds.append(x)

    o.stop_conditions = conds
    # Using new parameters.
    setArrheniusConstantsDNA23(o)
    return o


def calculateGateInputRate(gate_complex, input_complex, output_complex, trials=INCREMENT_TRIALS,  alt_output_complex=None):
    success_stop_condition = StopCondition(
        Options.STR_SUCCESS, [(output_complex, Options.dissocMacrostate, 0)])
    failed_stop_condition = StopCondition(
        Options.STR_FAILURE, [(input_complex, Options.dissocMacrostate, 0)])

    for x in [gate_complex, input_complex]:
        x.boltzmann_supersample = 25
        x.boltzmann_count = trials
        x.boltzmann_sample = True
        

    try:
        if alt_output_complex is None:
            raise TypeError
        alt_success_stop_condition = StopCondition(
            Options.STR_ALT_SUCCESS, [(output_complex, Options.dissocMacrostate, 0)])
        myMultistrand.setOptionsFactory6(getOptions, trials, DNA,
                                         gate_complex, input_complex,
                                         [success_stop_condition,
                                             alt_success_stop_condition],
                                         [failed_stop_condition])
    except TypeError:
        myMultistrand.setOptionsFactory6(getOptions, trials, DNA, gate_complex, input_complex, [
                                         success_stop_condition], [failed_stop_condition])

    return getRates()


def calculateBaseOutputRate(gate, trials=INCREMENT_TRIALS):
    rates = calculateGateInputRate(
        gate.gate_output_complex, gate.input_complex, gate.output_complex, trials)
    return rates


def calculateBaseFuelRate(gate, trials=INCREMENT_TRIALS,
                          ):
    rates = calculateGateInputRate(
        gate.gate_input_complex, gate.fuel_complex, gate.input_complex, trials)
    return rates


def calculateBaseThresholdRate(gate, trials=INCREMENT_TRIALS):
    rates = calculateGateInputRate(
        gate.threshold_complex, gate.input_complex, gate.threshold_free_waste_complex, trials)
    return rates


# MS: Adapted to return two rates within the rates object returned:
#     Both the rate for one complex being released but also the other,
#     as an alternative success condition for first step mode
def calculateGateGateLeak(gateA, gateB, trials=INCREMENT_TRIALS, material="DNA"):
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

    for x in [gateA_complex, gateB_complex]:
        x.boltzmann_supersample = 25
        x.boltzmann_count = trials
        x.boltzmann_sample = True
        

    myMultistrand.setOptionsFactory6(getOptions, trials, material,
                                     gateA_complex, gateB_complex,
                                     [success_stop_condition,
                                         alt_success_stop_condition],
                                     [failed_stop_condition])

    return getRates()

def getRates():
    myMultistrand.run()
    rates = SeesawRates(myMultistrand.results)
    print rates
    return rates


# MS: Calculates the leak between a gate and its fuel
def calculateGateFuelLeak(gate, trials=INCREMENT_TRIALS, material="DNA"):
    rates = calculateGateInputRate(gate.gate_output_complex, gate.fuel_complex, gate.output_complex, trials)
    return rates
