from __future__ import print_function

from multistrand.concurrent import myMultistrand, FirstStepRate, Bootstrap
from multistrand.experiment import standardOptions, hybridization

import sys, time

A_CONCENTRATION = 50e-9;


class customResult(object):
    
    def __init__(self):
    
        self.thing = 'x'
    
 
def first_step_simulation(strand_seq, trials, T=20.0, material="DNA"):
 
    print ("Running first step mode simulations for %s (with Boltzmann sampling)..." % (strand_seq))
        
    def getOptions(trials, material):
         
         
        o = standardOptions("First Step", tempIn=25.0, trials=200, timeOut = 0.1) 
        hybridization(o, strand_seq, trials)
        o.DNA23Metropolis()
        
          
        return o
      
    myMultistrand.setNumOfThreads(2)
    myMultistrand.setOptionsFactory2(getOptions, trials, material)
    myMultistrand.setTerminationCriteria(1000)
    myMultistrand.setLeakMode()
    myMultistrand.run()
    
    return myMultistrand.results    # this is a first step rate object


def compute(strand_seq):
    
    result = first_step_simulation(strand_seq, 200, T=25.0, material="DNA")
    
    rRate = result.k1()

    return "{:.2e}".format(float(rRate)), '999.0'


def computeAndWriteToCL(strand_seq, doBootstrap):
    
    result = first_step_simulation(strand_seq, 12000, T=25.0, material="DNA")
    print("The hybridization rate of ", strand_seq, " and the reverse complement is ", "{:.2e}".format(result.k1()), " /M /s", sep="")
    
    # check for two-stateness
    result.testForTwoStateness(100e-9)
    
    if(doBootstrap):
        
        bootstrap = Bootstrap(result, A_CONCENTRATION, computek1=True)
        bounds = bootstrap.ninetyFivePercentiles()
        
        print("Estimated 95% confidence interval: [","{:.2e}".format(bounds[0]),",","{:.2e}".format(bounds[1]),"] ", sep="")
        
        # print("The hybridization rate of ", mySequence, "and the reverse complement is ", myRates[0], "/M /s")


if(len(sys.argv) < 2):
    print("Please provide a DNA sequence as commandline argument")
    print("Add -bootstrap to do a boostrap ")
    print("Example: computeAnnealRate.py ATGCAGT -boostrap")
    exit()

start_time = time.time()

mySequence = sys.argv[1]
doBootstrap = False
if(len(sys.argv) > 2):
    if(str(sys.argv[2])=="-bootstrap"):
        doBootstrap = True

result = computeAndWriteToCL(mySequence, doBootstrap )


print ("Computing took %.4f s" % (time.time() - start_time))



