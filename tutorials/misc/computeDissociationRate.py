from __future__ import print_function

from multistrand.options import Options
from multistrand.concurrent import myMultistrand, FirstStepRate, Bootstrap
from multistrand.experiment import standardOptions, dissociation
from msArrhenius import setArrheniusConstantsDNA23

import sys, time

   
# Frits Dannenberg, Aug 2017.
# In order to compute dissociation rates for duplex, we can either compute the forward rate k+
# and simply compute k- from k+ / k-  = exp ( - dG / RT )
   
# In the following file, we simply simulate the actual dissocation time and compute k- = 1/t. 
 
def first_step_simulation(strand_seq, trials, temperature, material="DNA"):
 
    print ("Running first step mode simulations for %s (with Boltzmann sampling)..." % (strand_seq))
        
    def getOptions(trials, material):
         
        o = standardOptions(Options.firstPassageTime, temperature, trials, timeOut = 100.0) 
        dissociation(o, strand_seq, trials)
#         o.DNA23Metropolis()
        setArrheniusConstantsDNA23(o) # unreleased parameterization
          
        return o
      
    myMultistrand.setNumOfThreads(8)
    myMultistrand.setOptionsFactory2(getOptions, trials, material)
    myMultistrand.setTerminationCriteria(30)
    myMultistrand.setPassageMode()
    
    myMultistrand.run()
    
    return myMultistrand.results    # this is a first passage rate object


def compute(strand_seq):
    
    result = first_step_simulation(strand_seq, 16, 25.0, material="DNA")

    return result.k1()



def computeAndWriteToCL(strand_seq, doBootstrap):
    
    result = first_step_simulation(strand_seq, 16, 25.0, material="DNA")
    print("The dissociation rate of ", strand_seq, " and the reverse complement is ", "{:.2e}".format(result.k1()), " /M /s", sep="")
    
    if(doBootstrap):
        
        bootstrap = Bootstrap(result, 50e-9, computek1=True)
        bounds = bootstrap.ninetyFivePercentiles()
        
        print("Estimated 95% confidence interval: [","{:.2e}".format(bounds[0]),",","{:.2e}".format(bounds[1]),"] ", sep="")
        

if(len(sys.argv) < 2):
    print("Please provide a DNA sequence as commandline argument")
    print("Add -bootstrap to do a boostrap ")
    print("Example: computeDisscociationRate.py ATGCAGT -bootstrap")
    exit()

start_time = time.time()

mySequence = sys.argv[1]
doBootstrap = False
if(len(sys.argv) > 2):
    if(str(sys.argv[2])=="-bootstrap"):
        doBootstrap = True

result = computeAndWriteToCL(mySequence, doBootstrap )


print ("Computing took %.4f s" % (time.time() - start_time))



