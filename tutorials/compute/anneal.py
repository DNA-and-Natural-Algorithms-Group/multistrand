# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from multistrand.concurrent import MergeSim, FirstStepRate, Bootstrap
from multistrand.experiment import standardOptions, hybridization
from multistrand.options import Literals


num_threads = 10
num_trials = 1000
num_success = 500
myMultistrand = MergeSim()
A_CONCENTRATION = 50e-9


def first_step_simulation(strand_seq: str, trials: int, timeout: float,
                          temperature: float, sodium: float, material="DNA") -> FirstStepRate:
 
    print(f"\nRunning first step mode simulations for {strand_seq} "
          "(with Boltzmann sampling)...\n")

    def getOptions(_trials):
        o = standardOptions(simMode=Literals.first_step,
                            tempIn=temperature, trials=_trials, timeOut=timeout)
        o.sodium = sodium
        hybridization(o, strand_seq, _trials)
        # o.DNA23Arrhenius()
        o.JSDefault()
        o.substrate_type = material
        return o
      
    myMultistrand.setNumOfThreads(num_threads)
    myMultistrand.setOptionsFactory1(getOptions, trials)
    myMultistrand.setTerminationCriteria(num_success)
    myMultistrand.setFirstStepMode()
    # myMultistrand.setLeakMode()
    myMultistrand.run()
    return myMultistrand.results


def compute(strand_seq, temperature=25.0, sodium=1.0):
    return first_step_simulation(strand_seq, trials=num_trials, timeout=1.0,
                                 temperature=temperature, sodium=sodium)


def computeAndWriteToCL(strand_seq, doBootstrap, temperature=25.0, sodium=1.0):
    result = first_step_simulation(strand_seq, trials=num_trials, timeout=1.0,
                                   temperature=temperature, sodium=sodium)
    print(f"The hybridization rate of {strand_seq} and the reverse complement is "
          f"{result.k1():.2e} /M /s")
    
    # check for two-stateness
#     result.testForTwoStateness(100e-9)
    if(doBootstrap):
        bootstrap = Bootstrap(result, concentration=A_CONCENTRATION, N=1200, computek1=True)
        bounds = bootstrap.ninetyFivePercentiles()
        
        print(f"Estimated 95% confidence interval: [{bounds[0]:.2e},{bounds[1]:.2e}]")
