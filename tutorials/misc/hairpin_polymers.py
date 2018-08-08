"""

 Frits Dannenberg, July 2018


"""

import sys, os

import matplotlib
matplotlib.use('Agg')

import numpy as np

from matplotlib.ticker import ScalarFormatter

from multistrand.concurrent import MergeSim, FirstStepRate, Bootstrap
from multistrand.experiment import standardOptions, setBoltzmann
from multistrand.objects import StopCondition, Domain, Complex, Strand
from multistrand.options import Options, Literals

"""
Some global variables for this one-of script -- put in struct later
"""

NUMOFTHREADS = 2
NUMOFPATHS = 400
TERMINATIONCRIT = 2

RAVAN_H1 = "ACAAGTCAAACGCAGAGCGAGTGGAAGTAACGTGAGGCCACCTATTACCCTACCCTCACGTTACTTCCACTCGCTCTG"
RAVAN_H2 = "GAGTGGAAGTAACGTGAGGGTAGGGTAATAGGTGGCGCTCTGCGTTTGACTTGTGCCACCTATTACCCTACCCTCACG"
RAVAN_H3 = "GGTAGGGTAATAGGTGGCACAAGTCAAACGCAGAGCCTCACGTTACTTCCACTCGCTCTGCGTTTGACTTGTGCCACC"
RAVAN_I = "CTCACGTTACTTCCACTCGCTCTGCGTTTGACTT"

PENG_YIN = False


""" Note that the success complex has nonsense secondary structure (dotparen) 
    specification. This does not matter for dissociation stopping condition """


def ravan(options, trialsIn, selector):
    # we only allow first step mode at this point.
    
    if PENG_YIN == True:
        
        RAVAN_H1 = "GCTTGAGATGTTAGGGAGTAGTGCTCCAATCACAACGCACTACTCCCTAACATC"
        RAVAN_H2 = "AGGGAGTAGTGCGTTGTGATTGGAAACATCTCAAGCTCCAATCACAACGCACTA"
        RAVAN_H3 = "GTTGTGATTGGAGCTTGAGATGTTGCACTACTCCCTAACATCTCAAGCTCCAAT"
        RAVAN_I =  "GCACTACTCCCTAACATCTCAAGC"
    
    
    strand_H1 = Strand(name="H1", sequence=RAVAN_H1)
    strand_H2 = Strand(name="H2", sequence=RAVAN_H2)
    strand_H3 = Strand(name="H3", sequence=RAVAN_H3)
    strand_I = Strand(name="I", sequence=RAVAN_I)
    
    dotparen_H1 = "." * len(RAVAN_H1)
    dotparen_H2 = "." * len(RAVAN_H2)
    dotparen_H3 = "." * len(RAVAN_H3)
    dotparen_I = "." * len(RAVAN_I)
    
    """ No special secondary structure, just a hack to make a connected 
    complex composed of the correct set of strands """
    
    dotparen_I_H1 = '(' + ("." * (len(RAVAN_I) - 1)) + "+" + ')' + ("." * (len(RAVAN_H1) - 1)) 
    dotparen_I_H1_H2 = '(' + ("." * (len(RAVAN_I) - 1)) + "+" + ')' + ("." * (len(RAVAN_H1) - 2)) + "(" + "+" + ")" + ("." * (len(RAVAN_H2) - 1))
    dotparen_H1_H2_H3 = '(' + ("." * (len(RAVAN_H1) - 1)) + "+" + ')' + ("." * (len(RAVAN_H2) - 2)) + '(' + "+" + ')' + ("." * (len(RAVAN_H3) - 1))
    
    if selector == 0:
    
        target_complex = Complex(strands=[strand_H1], structure=dotparen_H1)
        complex_trigger = Complex(strands=[strand_I], structure=dotparen_I)
        
        success_complex = Complex(strands=[strand_I, strand_H1], structure=dotparen_I_H1)

    if selector == 1:

        target_complex = Complex(strands=[strand_I, strand_H1], structure=dotparen_I_H1)
        complex_trigger = Complex(strands=[strand_H2], structure=dotparen_H2)
        
        success_complex = Complex(strands=[strand_I, strand_H1, strand_H2], structure=dotparen_I_H1_H2)

    if selector == 2:

        target_complex = Complex(strands=[strand_I, strand_H1, strand_H2], structure=dotparen_I_H1_H2)
        complex_trigger = Complex(strands=[strand_H3], structure=dotparen_H3)
        
        success_complex = Complex(strands=[strand_H1, strand_H2, strand_H3], structure=dotparen_H1_H2_H3)

    setBoltzmann(complex_trigger, trialsIn)
    setBoltzmann(target_complex, trialsIn)
    
    stopSuccess = StopCondition(Literals.success, [(success_complex, Literals.dissoc_macrostate, 0)])
    stopFailed = StopCondition(Literals.failure, [(complex_trigger, Literals.dissoc_macrostate, 0)])

    # actually set the intial and stopping states    
    options.start_state = [complex_trigger, target_complex]
    options.stop_conditions = [stopSuccess, stopFailed]         
    

# trials has to be the first argument.
def genOptions(trialsIn, selector):

    stdOptions = standardOptions(Literals.first_step, tempIn=23.0, trials=trialsIn, timeOut=100.0)
    
    stdOptions.sodium = 0.05
    stdOptions.magnesium = 12.5e-3
    
    stdOptions.join_concentration = 5e-9
    
    stdOptions.simulation_time = 5.0
    
    stdOptions.DNA23Metropolis()
    ravan(stdOptions, trialsIn, selector)
    
    return stdOptions


def computeRate(selector):
    
    myMultistrand = MergeSim()
    
    myMultistrand.setOptionsFactory2(genOptions, NUMOFPATHS, selector) 
    myMultistrand.setNumOfThreads(NUMOFTHREADS)
    myMultistrand.setTerminationCriteria(TERMINATIONCRIT)
    myMultistrand.setLeakMode()  # the new leak object -- faster bootstrapping.
    
    myMultistrand.printTrajectory()
     
    return 0.0, 0.0, 0.0
    
    myMultistrand.run()
     
    myFSR = myMultistrand.results  
    low, high = myFSR.doBootstrap()
     
    return myFSR.k1(), low, high 


def simulateHairpinAssembly():

    for i in range(3):

        rate, low, high = computeRate(i)
        
        print " \n "
        print " reaction " + str(i)
        print " Rate = %.2e /M /s" % rate
        print " Low  = %.2e /M /s" % low
        print " High = %.2e /M /s" % high
        print " \n "    

    
# The actual main method
if __name__ == '__main__':

    if len(sys.argv) < 4:
        print("Please specify the number of successful trials, and number of threads and number of trials.")
        print("Example: python hadi_ravan.py  2 2 80")
        exit()

    TERMINATIONCRIT = int(sys.argv[1])
    NUMOFTHREADS = int(sys.argv[2])
    NUMOFPATHS = int(sys.argv[3])
    PENG_YIN = bool(sys.argv[4])
    
    simulateHairpinAssembly()

