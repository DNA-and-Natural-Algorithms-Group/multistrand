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

HOME = expanduser("~/workspace/multistrand")
MODEL = expanduser("~/workspace/multistrandPy")
sys.path.append(MODEL)
sys.path.append(HOME)

from multistrand.concurrent import myMultistrand, FirstStepRate
from multistrand.objects import StopCondition, Domain, Complex, Strand
from multistrand.options import Options
from multistrandPy.msArrhenius import setArrheniusConstantsDNA23
import numpy as np


ATIME_OUT = 10.0

myMultistrand.setNumOfThreads(4)


def getOptions(trials, material, complex1, complex2,
               success_stop_condition, failed_stop_condition, T=25):

    o = Options(simulation_mode="First Step", substrate_type=material,
                rate_method="Metropolis", num_simulations=trials,
                simulation_time=ATIME_OUT, temperature=T)

    o.start_state = [complex1, complex2]
    o.stop_conditions = [success_stop_condition, failed_stop_condition]
    # Using new parameters.
    setArrheniusConstantsDNA23(o)
    return o


def calculateBaseOutputRate(gate, trials=500,
                            material="DNA"):
    output_complex = gate.output_complex
    input_complex = gate.input_complex
    gate_complex = gate.gate_output_complex
    success_stop_condition = StopCondition(
        Options.STR_SUCCESS, [(output_complex, Options.dissocMacrostate, 0)])
    failed_stop_condition = StopCondition(
        Options.STR_FAILURE, [(input_complex, Options.dissocMacrostate, 0)])

    for x in [gate_complex, input_complex]:
        x.boltzmann_count = trials
        x.boltzmann_sample = True

    myMultistrand.setOptionsFactory6(getOptions, trials, material,
                                     gate_complex, input_complex,
                                     success_stop_condition,
                                     failed_stop_condition, )

    myMultistrand.run()
    dataset = myMultistrand.results
    myFSR = FirstStepRate(dataset, 5e-9)
    print("Was success :  %i  " % myFSR.nForward)
    print("Was failure :  %i  " % myFSR.nReverse)
    # print("Total runs  :  %i  " % myFSR.nTotal)
    print("k1          :  %.2f /M /s \n " % myFSR.k1())
    return myFSR.k1()


def calculateBaseFuelRate(gate, trials=500,
                          material="DNA"):
    # Not fully sure if this will work, but structure of the reaction
    # is similar between the two.
    k = calculateBaseOutputRate(gate, trials, material)
    return k


def calculateBaseThresholdRate(gate, trials=500, material="DNA"):
    # MS: Similarly to above, this needs to be tested - I am unsure as
    #     to whether it will work
    k = calculateBaseOutputRate(gate, trials, material)
    return k


# MS: Please note that this method assume that the output strand
#     of gateB, note gate A is released and returns that respective rate
def calculateGateGateLeak(gateA, gateB, trials=500, material="DNA"):
    # define stop conditions
    gateA_complex = gateA.gate_output_complex
    gateB_complex = gateB.gate_output_complex
    leak_complex = gateB.output_complex
    success_stop_condition = StopCondition(
        Options.STR_SUCCESS, [(leak_complex,
                               Options.dissocMacrostate, 0)])
    failed_stop_condition = StopCondition(
        Options.STR_FAILURE, [(gateA_complex,
                               Options.dissocMacrostate, 0)])

    for x in [gateA_complex, gateB_complex]:
        x.boltzmann_count = trials
        x.boltzmann_sample = True

    myMultistrand.setOptionsFactory6(getOptions, trials, material,
                                     gateA_complex, gateB_complex,
                                     success_stop_condition,
                                     failed_stop_condition)

    myMultistrand.run()
    dataset = myMultistrand.results
    myFSR = FirstStepRate(dataset, 5e-9)
    print("Was success :  %i  " % myFSR.nForward)
    print("Was failure :  %i  " % myFSR.nReverse)
    # print("Total runs  :  %i  " % myFSR.nTotal)
    print("k1          :  %.2f /M /s \n " % myFSR.k1())
    return myFSR.k1()


# MS: Calculates the leak between a gate and its fuel
def calculateGateFuelLeak(gate, trials=500, material="DNA"):
    gate_complex = gate.gate_output_complex
    fuel_complex = gate.fuel_complex
    leak_complex = gate.gate_output_complex
    success_stop_condition = StopCondition(
        Options.STR_SUCCESS, [(leak_complex,
                               Options.dissocMacrostate, 0)])
    failed_stop_condition = StopCondition(
        Options.STR_FAILURE, [(fuel_complex,
                               Options.dissocMacrostate, 0)])

    for x in [gate_complex, fuel_complex]:
        x.boltzmann_count = trials
        x.boltzmann_sample = True

    myMultistrand.setOptionsFactory6(getOptions, trials, material,
                                     gate_complex, fuel_complex,
                                     success_stop_condition,
                                     failed_stop_condition)

    myMultistrand.run()
    dataset = myMultistrand.results
    myFSR = FirstStepRate(dataset, 5e-9)
    print("Was success :  %i  " % myFSR.nForward)
    print("Was failure :  %i  " % myFSR.nReverse)
    # print("Total runs  :  %i  " % myFSR.nTotal)
    print("k1          :  %.2f /M /s \n " % myFSR.k1())
    return myFSR.k1()
