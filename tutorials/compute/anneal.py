from __future__ import print_function

from multistrand.concurrent import MergeSim, FirstStepRate, Bootstrap
from multistrand.experiment import standardOptions, hybridization
from multistrand.options import Options
from multistrand.utils import seqComplement

import sys, time, math

A_CONCENTRATION = 50e-9;
GAS_CONSTANT_R = 0.0019872036

myMultistrand = MergeSim()
   
 
def first_step_simulation(strand_seq, trials, temperature=25.0, sodium = 1.0, material="DNA"):
 
    print ("Running first step mode simulations for %s (with Boltzmann sampling)..." % (strand_seq))
        
    def getOptions(trials, material, temperature=25.0, sodium = 1.0):
         
        o = standardOptions(Options.firstStep, tempIn=temperature, trials=200, timeOut = 1.0)
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
    print("The hybridization rate of ", strand_seq, " and the reverse complement is ", "{:.2e}".format(result.k1()), " /M /s", sep="")
    
    # check for two-stateness
#     result.testForTwoStateness(100e-9)
    
    if(doBootstrap):
        
        bootstrap = Bootstrap(result, concentration=A_CONCENTRATION, N=1200, computek1=True)
        bounds = bootstrap.ninetyFivePercentiles()
        
        print("Estimated 95% confidence interval: [","{:.2e}".format(bounds[0]),",","{:.2e}".format(bounds[1]),"] ", sep="")
        
        
# compute dissociation from association rate
def computeDissociationAndWriteToCL(strand_seq, doBootstrap):

    
    
    
    result = first_step_simulation(strand_seq, 1200, material="DNA")
    
    seq = strand_seq
    seqC = seqComplement(seq)
    
    temp = 273.15+ 25.0 # this is just for NUPACK calls, setting temperature is not yet implemented properly.
    
    ## We only import nupack bindings here because it will print an welcom message 
    from nupack import pfunc   
    dG = pfunc([seq, seqC], [1,2], T=(temp-273.15), material="dna")
        
    print ("Using dG = " + "{:.2e}".format(dG) + " kcal/mol, and k+ = " + "{:.2e}".format(result.k1()) + " /M /s to compute the dissociation rate." )
        
    kMinus = result.k1() * math.exp( dG / ( GAS_CONSTANT_R * temp) ) 

       
    print("The dissociation rate of ", strand_seq, " and the reverse complement is ", "{:.2e}".format(kMinus), " /M /s \n", sep="")
    
    # check for two-stateness
#     result.testForTwoStateness(100e-9)
    
    if(doBootstrap):
        
        low, high = result.doBootstrap(NIn=1200)
        kMinusLow = low * math.exp( dG / ( GAS_CONSTANT_R * temp) ) 
        kMinusHigh = high * math.exp( dG / ( GAS_CONSTANT_R * temp) ) 

        print("Estimated 95% confidence interval: [","{:.2e}".format(kMinusLow),",","{:.2e}".format(kMinusHigh),"] ", sep="")
