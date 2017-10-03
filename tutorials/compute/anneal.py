from __future__ import print_function

from multistrand.concurrent import myMultistrand, FirstStepRate, Bootstrap
from multistrand.experiment import standardOptions, hybridization
from multistrand.options import Options
from multistrand.objects import Strand
# from msArrhenius import setArrheniusConstantsDNA23
from rawdata.readensemble import setArrParams

import sys, time

A_CONCENTRATION = 50e-9;

   
 
def first_step_simulation(strand_seq, trials, temperature=25.0, sodium = 1.0, material="DNA"):
 
    print ("Running first step mode simulations for %s (with Boltzmann sampling)..." % (strand_seq))
        
    def getOptions(trials, material, temperature=25.0, sodium = 1.0):
         
        o = standardOptions(Options.firstStep, tempIn=temperature, trials=200, timeOut = 1.0)
        o.sodium = sodium 
        hybridization(o, strand_seq, trials)
        o.DNA23Metropolis()
#         setArrheniusConstantsDNA23(o)
        setArrParams(o, 92)
          
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
    result.testForTwoStateness(100e-9)
    
    if(doBootstrap):
        
        bootstrap = Bootstrap(result, concentration=A_CONCENTRATION, N=1200, computek1=True)
        bounds = bootstrap.ninetyFivePercentiles()
        
        print("Estimated 95% confidence interval: [","{:.2e}".format(bounds[0]),",","{:.2e}".format(bounds[1]),"] ", sep="")
        
        
# compute dissociation from association rate
def computeDissociationAndWriteToCL(strand_seq, doBootstrap):
    
    result = first_step_simulation(strand_seq, 1200, material="DNA")
    
    seq = strand_seq
     
    topStrand = Strand(seq)
    seqC = (topStrand.C).seq
    
    
    dG = pfunc([seq, seqC], [1,2], T=(temp-273.15), material="dna")
    print (str(dG)) 
        
    kMinus = predicted.k1() * math.exp( dG / ( GAS_CONSTANT_R * temp) ) 
    kMinusLow = low * math.exp( dG / ( GAS_CONSTANT_R * temp) ) 
    kMinusHigh = high * math.exp( dG / ( GAS_CONSTANT_R * temp) ) 
    
       
    print("The dissociation rate of ", strand_seq, " and the reverse complement is ", "{:.2e}".format(result.k1()), " /M /s", sep="")
    
    # check for two-stateness
    result.testForTwoStateness(100e-9)
    
    if(doBootstrap):
        
        bootstrap = Bootstrap(result, concentration=A_CONCENTRATION, N=1200, computek1=True)
        bounds = bootstrap.ninetyFivePercentiles()
        
        print("Estimated 95% confidence interval: [","{:.2e}".format(bounds[0]),",","{:.2e}".format(bounds[1]),"] ", sep="")
        

