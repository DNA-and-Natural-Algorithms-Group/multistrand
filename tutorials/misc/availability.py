# Frits Dannenberg, Aug 2017

# simulates ONLY the base reaction.

# Availability: A Metric for Nucleic Acid Strand Displacement Systems. Xiaoping Olson, Shohei Kotani, Jennifer E. Padilla, 
# Natalya Hallstrom, Sara Goltry, Jeunghoon Lee, Bernard Yurke, William L. Hughes, and Elton Graugnard


from multistrand.objects import StopCondition, Domain, Complex, Strand
from multistrand.options import Options
from multistrand.concurrent import myMultistrand
from multistrand.experiment import makeComplex, standardOptions, setBoltzmann

from msArrhenius import setArrheniusConstantsDNA23

myMultistrand.setNumOfThreads(5)

def doExperiment(trials):

    # complexes
    
    seqOutput = "CTACTTTCACCCTACGTCTCCAACTAACTTACGG"
    seqSignal = "CCACATACATCATATTCCCTCATTCAATACCCTACG"
    seqBackbone = "TGGAGACGTAGGGTATTGAATGAGGGCCGTAAGTT"
    
    # output + signal + backbone
    complex_dotparen = "."*10 + "("* (len(seqOutput)-10) 
    complex_dotparen += "+" + "."*16 + "("* (len(seqSignal)-16)  
    complex_dotparen += "+" + "." * 6  + ")"* (len(seqBackbone)-6)
    myGate = makeComplex([seqOutput, seqSignal, seqBackbone], complex_dotparen)
    
    
    seqFuel = "CCTACGTCTCCAACTAACTTACGGCCCTCATTCAATACCCTACG"
    myFuel = makeComplex([seqFuel],"."*len(seqFuel))
    
    myLeakedSignal = makeComplex([seqSignal], "."*len(seqSignal))
    
    for x in [myGate, myFuel]:
        setBoltzmann(x, trials, supersample=25)
    
    # stopping states
    
    # the failed state is when the fuel releases again
    # the success state is when the leaked signal is observed.
    # we use dissocMacrostate - this only checks the presence of strands in the complex, and 
    # does not depend on dotparens structure.
    
    success_stop_condition = StopCondition(Options.STR_SUCCESS, [(myLeakedSignal, Options.dissocMacrostate, 0)])
    failed_stop_condition = StopCondition(Options.STR_FAILURE, [(myFuel, Options.dissocMacrostate, 0)])
    
    
    # options
    
    stdOptions = standardOptions(trials=trials)
     
    stdOptions.start_state = [myGate, myFuel]
    stdOptions.stop_conditions = [success_stop_condition, failed_stop_condition]

    # buffer
    
    stdOptions.sodium = 0.05
    stdOptions.magnesium = 0.0125 ##  believed to be 11.5 mM effectively -- using 12.5 mM anyway
    
    # rate model
    setArrheniusConstantsDNA23(stdOptions)
    
    return stdOptions

    
    
# actually calling multistrand

myMultistrand.setOptionsFactory1(doExperiment, 15000)
myMultistrand.setTerminationCriteria(terminationCount=20)
myMultistrand.setLeakMode()

myMultistrand.run()
    

