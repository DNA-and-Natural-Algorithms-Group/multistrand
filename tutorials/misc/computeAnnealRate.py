from __future__ import print_function

from multistrand.concurrent import myMultistrand, FirstStepRate, Bootstrap
from multistrand.experiment import setSaltGao2006, standardOptions, hybridization
from multistrand.utils import DNA23Metropolis


import sys

A_CONCENTRATION = 50e-9;


class customResult(object):
    
    def __init__(self):
    
        self.thing = 'x'
    
 
def first_step_simulation(strand_seq, trials, T=20.0, material="DNA"):
 
    print ("Running first step mode simulations for %s (with Boltzmann sampling)..." % (strand_seq))
        
    def getOptions(trials, material):
         
         
        o = standardOptions("First Step", tempIn=25.0, trials=200, timeOut = 0.001) 
        hybridization(o, strand_seq, trials)
        setSaltGao2006(o)
        DNA23Metropolis(o)
        
          
        return o
      
    myMultistrand.setOptionsFactory2(getOptions, trials, material)
    myMultistrand.setTerminationCriteria(FirstStepRate())
    myMultistrand.run()
    dataset = myMultistrand.results

    
    return FirstStepRate(dataset, A_CONCENTRATION) #, myMultistrand.runTime


def compute(strand_seq):
    
    result = first_step_simulation(strand_seq, 200, T=25.0, material="DNA")
    
    rRate = result.k1()

    return "{:.2e}".format(float(rRate)), '999.0'


def computeAndWriteToCL(strand_seq, doErrorBars):
    

    # myRates = computeWError(mySequence)
    
    result = first_step_simulation(strand_seq, 800, T=25.0, material="DNA")
    print("The hybridization rate of ", strand_seq, " and the reverse complement is ", "{:.2e}".format(result.kEff()), " /M /s", sep="")
    

    # check for two-stateness
    result.testForTwoStateness()
    
    if(doErrorBars):
        
        bootstrap = Bootstrap(result, A_CONCENTRATION, computek1=True)
        bounds = bootstrap.ninetyFivePercentiles()
        
        print("Estimated 95% confidence interval: [","{:.2e}".format(bounds[0]),",","{:.2e}".format(bounds[1]),"] ", sep="")
        
        # print("The hybridization rate of ", mySequence, "and the reverse complement is ", myRates[0], "/M /s")


if(len(sys.argv) < 2):
    print("Please provide a DNA sequence as commandline argument")
    print("Adding an additional argument toggles errorbound computation")
    print("Example: computeAnnealRate.py ATGCAGT bounds")
    exit()

mySequence = sys.argv[1]
result = computeAndWriteToCL(mySequence, len(sys.argv) > 2)




