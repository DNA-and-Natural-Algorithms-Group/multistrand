from __future__ import print_function

from multistrand.concurrent import myMultistrand, FirstStepRate, Bootstrap
from multistrand.experiment import standardOptions, hybridization
from multistrand.options import Options
from msArrhenius import setArrheniusConstantsDNA23

import sys, time

A_CONCENTRATION = 50e-9;

   
 
def first_step_simulation(strand_seq, trials, temperature=25.0, sodium = 1.0, material="DNA"):
 
    print ("Running first step mode simulations for %s (with Boltzmann sampling)..." % (strand_seq))
        
    def getOptions(trials, material, temperature=25.0, sodium = 1.0):
         
        o = standardOptions(Options.firstStep, tempIn=temperature, trials=200, timeOut = 1.0)
        o.sodium = sodium 
        hybridization(o, strand_seq, trials)
#         o.DNA23Metropolis()
        setArrheniusConstantsDNA23(o)
          
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
    
    result = first_step_simulation(strand_seq, 1200, T=25.0, material="DNA")
    print("The hybridization rate of ", strand_seq, " and the reverse complement is ", "{:.2e}".format(result.k1()), " /M /s", sep="")
    
    # check for two-stateness
    result.testForTwoStateness(100e-9)
    
    if(doBootstrap):
        
        bootstrap = Bootstrap(result, A_CONCENTRATION, computek1=True)
        bounds = bootstrap.ninetyFivePercentiles()
        
        print("Estimated 95% confidence interval: [","{:.2e}".format(bounds[0]),",","{:.2e}".format(bounds[1]),"] ", sep="")
        
        # print("The hybridization rate of ", mySequence, "and the reverse complement is ", myRates[0], "/M /s")


# if(len(sys.argv) < 2):
#     print("Please provide a DNA sequence as commandline argument")
#     print("Add -bootstrap to do a boostrap ")
#     print("Example: computeAnnealRate.py ATGCAGT -bootstrap")
#     exit()
# 
# start_time = time.time()
# 
# mySequence = sys.argv[1]
# doBootstrap = False
# if(len(sys.argv) > 2):
#     if(str(sys.argv[2])=="-bootstrap"):
#         doBootstrap = True
# 
# result = computeAndWriteToCL(mySequence, doBootstrap )
# 
# 
# print ("Computing took %.4f s" % (time.time() - start_time))



