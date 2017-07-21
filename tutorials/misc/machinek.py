# Frits Dannenberg, July 2017

# Simulation of Figure 2d in 
# Programmable energy landscapes for kinetic control of DNA strand displacement
# RRF Machinek, TE Ouldridge, NEC Haley, J Bath & AJ Turberfield. Nature communications

import sys

from multistrand.concurrent import myMultistrand, FirstStepRate, Bootstrap
from multistrand.experiment import standardOptions, setBoltzmann
from multistrand.objects import StopCondition, Domain, Complex, Strand
from multistrand.options import Options



# Figure 2d has 3x12 = 36 rates plotted. 
# Input: 0 <= selector < 36  

# range 0 -11:  6 nt toehold
# range 12-23:  7 nt toehold
# range 24-36: 10 nt toehold

# order of mismatch position:
# perfect - 2 - 3 - 4 - 5 - 6 - 7 - 8 - 9 -10 - 12 - 14
   
def machinek2014(options, selector, trialsIn):
    # we only allow first step mode at this point.
    
    # these are the sequences we need to build the dot-parens
    incumbent = ""
    target = ""
    invader = ""
    
    toeholdSelect =  selector / 12 
    mismatchSelect = selector % 12
    
    positionSelector = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14]    
    mismatchSelect = positionSelector[mismatchSelect]
    
    # decide on toehold sequence
    toeholdSeq = "ATGTGG"       # 6 nt toehold option
    
    if toeholdSelect == 1 :
        toeholdSeq = "ATGTGGA"  # 7 nt toehold option
    if toeholdSelect == 2 :
        toeholdSeq = "ATGTGGAGGG" # 10 nt toehold option


    # determine the incumbent, target and invader sequences
    # FD: copy-pasting supplementary Table 6 directly

    if mismatchSelect == 0 or mismatchSelect == 2 or mismatchSelect == 12 or mismatchSelect == 14  :
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTTTGAGGTTGA"
        target = "CCCTCCACATTCAACCTCAAACTCACC"
        
        if mismatchSelect == 0: #perfect
            invader = "GGTGAGTTTGAGGTTGA"

        if mismatchSelect == 2:
            invader = "GGTGAGTTTGAGGTTCA"
        
        if mismatchSelect == 12:
            invader = "GGTGACTTTGAGGTTGA"
            
        if mismatchSelect == 14:
            invader = "GGTCAGTTTGAGGTTGA"
    
    if mismatchSelect == 3:
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTTTGAGGTGAT"
        target = "CCCTCCACATATCACCTCAAACTCACC"
        invader = "GGTGAGTTTGAGGTCAT"
        
    if mismatchSelect == 4:
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTTTGAGTGAGT"
        target = "CCCTCCACATACTCACTCAAACTCACC"
        invader = "GGTGAGTTTGAGTCAGT"


    if mismatchSelect == 5:
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTTTGATGAGGT"
        target = "CCCTCCACATACCTCATCAAACTCACC"
        invader = "GGTGAGTTTGATCAGGT"

    if mismatchSelect == 6:
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTTTGTGAAGGT"
        target = "CCCTCCACATACCTTCACAAACTCACC"
        invader = "GGTGAGTTTGTCAAGGT"
        
    if mismatchSelect == 7:
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTTTTGAGAGGT"
        target = "CCCTCCACATACCTCTCAAAACTCACC"
        invader = "GGTGAGTTTTCAGAGGT"

    if mismatchSelect == 8:
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTTTGATGAGGT"
        target = "CCCTCCACATACCTCATCAAACTCACC"
        invader = "GGTGAGTTTCATGAGGT"
        
    if mismatchSelect == 9:
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTTGATTGAGGT"
        target = "CCCTCCACATACCTCAATCAACTCACC"
        invader = "GGTGAGTTCATTGAGGT"
        
    if mismatchSelect == 10:
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTGATTTGAGGT"
        target = "CCCTCCACATACCTCAAATCACTCACC"
        invader = "GGTGAGTCATTTGAGGT"

    invader = invader + toeholdSeq
    
#     print "incumb:" +  incumbent
#     print "target:" + target
#     print "invader:" + invader
    
    # set up the actual complexes
    strandIncumbent = Strand(name="incumbent", sequence=incumbent)
    strandTarget = Strand(name="target", sequence=target)
    strandInvader = Strand(name="invader", sequence=invader)
    
    intialDotParen = '.'*16 + '('*17 + "+" + '.'*10 + ')'* 17  
    intialInvaderDotParen = '.'*len(invader)
    successDotParen = '.'*33
    
    initialComplex = Complex(strands=[strandIncumbent, strandTarget], structure=intialDotParen)
    initialInvader = Complex(strands=[strandInvader], structure=intialInvaderDotParen)
    successComplex = Complex(strands=[strandIncumbent], structure=successDotParen)
    
        # Turns Boltzmann sampling on for this complex and also does sampling more efficiently by sampling 'trials' states.
    if(trialsIn > 0):
        setBoltzmann(initialComplex, trialsIn)
        setBoltzmann(initialInvader, trialsIn)

    stopSuccess = StopCondition(Options.STR_SUCCESS, [(successComplex, Options.dissocMacrostate, 0)])
    stopFailed = StopCondition(Options.STR_FAILURE, [(initialComplex, Options.dissocMacrostate, 0)])
    
    # actually set the intial and stopping states    
    options.start_state = [initialComplex, initialInvader]
    options.stop_conditions = [stopSuccess, stopFailed]         
    
        
    

# trials has to be the first argument.
def genOptions( trialsIn, select):

    stdOptions = standardOptions(Options.firstStep, tempIn=25.0, trials=trialsIn, timeOut=0.1)
    machinek2014(stdOptions, select, trialsIn)
    return stdOptions


def generateGraph(trials = 15):
    
    myMultistrand.setOptionsFactory2(genOptions, trials, 0) 
    myMultistrand.setNumOfThreads(4)
    
    myMultistrand.run()
    myFSR = FirstStepRate(myMultistrand.results, 5e-9) # 5 nM concentration
    
    print myFSR
    
    
    
# The actual main method
if __name__ == '__main__':

    if len(sys.argv) > 1:
        print("No commandline arguments required. Aborting program.")
        exit()
    
    generateGraph(150)



























