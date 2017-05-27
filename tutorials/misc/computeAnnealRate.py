from __future__ import print_function

from multistrand.concurrent import myMultistrand, FirstStepRate
from multistrand.experiment import setSaltGao2006, standardOptions, hybridization
from multistrand.utils import DNA23Metropolis


import sys

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
    myMultistrand.run()
    dataset = myMultistrand.results

    
    return FirstStepRate(dataset, 50e-9).k1() #, myMultistrand.runTime


def compute(strand_seq):
    
    result = first_step_simulation(strand_seq, 200, T=25.0, material="DNA")
    
    rRate = result    

    return "{:.2e}".format(float(rRate)), '999.0'


if(len(sys.argv) < 2):
    print("Please provide a DNA sequence as commandline argument")
    exit()

# print("Input: ", sys.argv[1])
mySequence = sys.argv[1]
myRate = compute(mySequence)[0]
print("The hybridization rate of ", mySequence, "and it's complement is ", myRate, "/M /s")

