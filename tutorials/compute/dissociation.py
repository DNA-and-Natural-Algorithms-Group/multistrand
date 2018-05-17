from __future__ import print_function

from multistrand.options import Options, Literals
from multistrand.concurrent import MergeSim, FirstStepRate, Bootstrap
from multistrand.experiment import standardOptions, dissociation

import sys, time

myMultistrand = MergeSim()
   
# Frits Dannenberg, Aug 2017.
# In order to compute dissociation rates for duplex, we can either compute the forward rate k+
# and simply compute k- from k+ / k-  = exp ( - dG / RT )
   
# In the following file, we simply simulate the actual dissocation time and compute k- = 1/t. 
 
def first_step_simulation(strand_seq, trials, temperature):
 
    print ("Running first passage time simulations for %s (with Boltzmann sampling)..." % (strand_seq))
        
    def getOptions(trials):
         
        o = standardOptions(Literals.first_passage_time, temperature, trials, timeOut = 1000.0) 
        dissociation(o, strand_seq, trials)
        o.DNA23Metropolis()
#         setArrheniusConstantsDNA23(o) # unreleased parameterization
          
        return o
      
    myMultistrand.setNumOfThreads(8)
    myMultistrand.setOptionsFactory1(getOptions, trials)
    myMultistrand.setTerminationCriteria(500)
    myMultistrand.setPassageMode()
    
    myMultistrand.run()
    
    return myMultistrand.results    # this is a first passage rate object


def compute(strand_seq, temperature=25):
    
    return first_step_simulation(strand_seq, 16,temperature)




def computeAndWriteToCL(strand_seq, doBootstrap):
    
    result = first_step_simulation(strand_seq, 64, 25.0)
    print("The dissociation rate of ", strand_seq, " and the reverse complement is ", "{:.2e}".format(result.k1()), " /s", sep="")
    
    if(doBootstrap):
        
        bootstrap = Bootstrap(result, 50e-9, computek1=True)
        bounds = bootstrap.ninetyFivePercentiles()
        
        print("Estimated 95% confidence interval: [","{:.2e}".format(bounds[0]),",","{:.2e}".format(bounds[1]),"] ", sep="")




