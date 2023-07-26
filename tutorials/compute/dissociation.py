
from multistrand.options import Literals
from multistrand.concurrent import MergeSim, FirstPassageRate, Bootstrap
from multistrand.experiment import standardOptions, dissociation

"""
Frits Dannenberg, Aug 2017.

In order to compute dissociation rates for duplex, we can either compute the
forward rate k+ and simply compute k- from k+ / k- = exp ( - dG / RT )

In the following file, we simply simulate the actual dissocation time and
compute k- = 1/t.
"""

num_threads = 10
num_trials = 1000
num_success = 500
myMultistrand = MergeSim()
   

def first_step_simulation(strand_seq: str, trials: int, timeout: float,
                          temperature: float, sodium: float, material="DNA") -> FirstPassageRate:
 
    print(f"\nRunning first passage time simulations for {strand_seq} "
          "(with Boltzmann sampling)...\n")
        
    def getOptions(_trials):
        o = standardOptions(simMode=Literals.first_passage_time,
                            tempIn=temperature, trials=_trials, timeOut=timeout)
        o.sodium = sodium
        dissociation(o, strand_seq, _trials)
        o.DNA23Arrhenius()
        o.substrate_type = material
        return o
      
    myMultistrand.setNumOfThreads(num_threads)
    myMultistrand.setOptionsFactory1(getOptions, trials)
    myMultistrand.setTerminationCriteria(num_success)
    myMultistrand.setPassageMode()
    myMultistrand.run()
    return myMultistrand.results


def compute(strand_seq, temperature=25.0, sodium=1.0):
    return first_step_simulation(strand_seq, trials=num_trials, timeout=10.0,
                                 temperature=temperature, sodium=sodium)


def computeAndWriteToCL(strand_seq, doBootstrap, temperature=25.0, sodium=1.0):
    result = first_step_simulation(strand_seq, trials=num_trials, timeout=10.0,
                                   temperature=temperature, sodium=sodium)
    print(f"The dissociation rate of {strand_seq} and the reverse complement is "
          f"{result.k1():.2e} /s")
    
    if(doBootstrap):
        bootstrap = Bootstrap(result, 50e-9, computek1=True)
        bounds = bootstrap.ninetyFivePercentiles()
        
        print(f"Estimated 95% confidence interval: [{bounds[0]:.2e},{bounds[1]:.2e}]")
