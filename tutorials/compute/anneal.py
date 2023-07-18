
import math

from multistrand.concurrent import MergeSim, Bootstrap
from multistrand.experiment import standardOptions, hybridization
from multistrand.options import Literals
from multistrand.utils import seqComplement


A_CONCENTRATION = 50e-9;
GAS_CONSTANT_R = 0.0019872036

myMultistrand = MergeSim()
   
 
def first_step_simulation(strand_seq, trials, temperature=25.0, sodium = 1.0, material="DNA"):
 
    print(f"Running first step mode simulations for {strand_seq} (with Boltzmann sampling)...")
        
    def getOptions(trials, material, temperature=25.0, sodium = 1.0):
        o = standardOptions(Literals.first_step, tempIn=temperature, trials=200, timeOut = 1.0)
        o.sodium = sodium 
        hybridization(o, strand_seq, trials)
        o.DNA23Metropolis()
#         setArrParams(o, 92) # the best DNA23 parameter set 
        return o
      
    myMultistrand.setNumOfThreads(6)
    myMultistrand.setOptionsFactory4(getOptions, trials, material, temperature, sodium)
    myMultistrand.setTerminationCriteria(500)
    myMultistrand.setLeakMode()
    myMultistrand.run()
    
    return myMultistrand.results    # this is a first step rate object


def compute(strand_seq, temperature=25.0, sodium = 1.0):
    return first_step_simulation(strand_seq, 240, temperature,  sodium, material="DNA")


def computeAndWriteToCL(strand_seq, doBootstrap):
    
    result = first_step_simulation(strand_seq, 1200, material="DNA")
    print(f"The hybridization rate of {strand_seq} and the reverse complement is "
          f"{result.k1():.2e} /M /s")
    
    # check for two-stateness
#     result.testForTwoStateness(100e-9)
    if(doBootstrap):
        bootstrap = Bootstrap(result, concentration=A_CONCENTRATION, N=1200, computek1=True)
        bounds = bootstrap.ninetyFivePercentiles()
        
        print(f"Estimated 95% confidence interval: [{bounds[0]:.2e},{bounds[1]:.2e}]")
        
        
# compute dissociation from association rate
def computeDissociationAndWriteToCL(strand_seq, doBootstrap):

    result = first_step_simulation(strand_seq, 1200, material="DNA")
    
    seq = strand_seq
    seqC = seqComplement(seq)
    
    temp = 273.15+ 25.0 # this is just for NUPACK calls, setting temperature is not yet implemented properly.
    
    ## We only import nupack bindings here because it will print an welcom message 
    from nupack import pfunc   
    dG = pfunc([seq, seqC], [1,2], T=(temp-273.15), material="dna")
    print(f"Using dG = {dG:.2e} kcal/mol, and k+ = {result.k1():.2e} /M /s "
          f"to compute the dissociation rate.")
        
    kMinus = result.k1() * math.exp(dG / (GAS_CONSTANT_R * temp))
    print(f"The dissociation rate of {strand_seq} and the reverse complement is "
          f"{kMinus:.2e} /M /s \n")
    
    # check for two-stateness
#     result.testForTwoStateness(100e-9)
    
    if(doBootstrap):
        low, high = result.doBootstrap(NIn=1200)
        kMinusLow = low * math.exp( dG / ( GAS_CONSTANT_R * temp) ) 
        kMinusHigh = high * math.exp( dG / ( GAS_CONSTANT_R * temp) ) 

        print(f"Estimated 95% confidence interval: [{kMinusLow:.2e},{kMinusHigh:.2e}]")
