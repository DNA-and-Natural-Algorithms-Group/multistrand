# Frits Dannenberg, Aug 2017

# simulates ONLY the base reaction.

# Availability: A Metric for Nucleic Acid Strand Displacement Systems. Xiaoping Olson, Shohei Kotani, Jennifer E. Padilla, 
# Natalya Hallstrom, Sara Goltry, Jeunghoon Lee, Bernard Yurke, William L. Hughes, and Elton Graugnard


from multistrand.objects import StopCondition, Domain, Complex, Strand
from multistrand.options import Options
from multistrand.concurrent import myMultistrand
from multistrand.experiment import makeComplex, standardOptions, setBoltzmann


myMultistrand.setNumOfThreads(8)

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
        setBoltzmann(x, trials, supersample=10)
    
    # stopping states
    
    # the failed state is when the fuel releases again
    # the success state is when the leaked signal is observed.
    # we use dissocMacrostate - this only checks the presence of strands in the complex, and 
    # does not depend on dotparens structure.
    
    successStopping = StopCondition(Options.STR_SUCCESS, [(myLeakedSignal, Options.dissocMacrostate, 0)])
    failedStopping = StopCondition(Options.STR_FAILURE, [(myFuel, Options.dissocMacrostate, 0)])
    
    
    # options
    
    stdOptions = standardOptions(trials=trials)
     
    stdOptions.start_state = [myGate, myFuel]
    stdOptions.stop_conditions = [successStopping, failedStopping]

    # buffer and temperature
    stdOptions.sodium = 0.05
    stdOptions.magnesium = 0.0125 ##  believed to be 11.5 mM effectively -- using 12.5 mM anyway
    
    stdOptions.temperature = 25.0  # can run at higher temperature to increase leak rate.
    
    # rate model
    stdOptions.DNA23Metropolis()
    #setArrheniusConstantsDNA23(stdOptions)
    stdOptions.simulation_time = 10.0
    
    
    return stdOptions

    
    
# actually calling multistrand

myMultistrand.setOptionsFactory1(doExperiment, 5*20*4)
myMultistrand.setTerminationCriteria(terminationCount=6)
myMultistrand.setLeakMode()

myMultistrand.run()
    

