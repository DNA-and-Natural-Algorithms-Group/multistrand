
## Often recurring experimental setups
from multistrand.objects import Complex, Domain, Strand, StopCondition
from multistrand.options import Options


def setBoltzmann(complexIn, trials):
    
    complexIn.boltzmann_count = trials
    complexIn.boltzmann_sample = True


# easy handle for options creation
def standardOptions(simMode = Options.firstStep, tempIn = 25.0, trials=10, timeOut= 0.1):
    
    output = Options(simulation_mode=simMode,
                      rate_method=Options.metropolis,
                      num_simulations=trials,
                      simulation_time=timeOut,
                      temperature=tempIn
                      )

    output.DNA23Metropolis()
    
    return output



def hybridization(options, mySeq, myTrials=0, doFirstPassage=False):
                
    # Using domain representation makes it easier to write secondary structures.
    onedomain = Domain(name="itall", sequence=mySeq)
    top = Strand(name="top", domains=[onedomain])
    bot = top.C
    
    # Note that the structure is specified to be single stranded, but this will be over-ridden when Boltzmann sampling is turned on.
    startTop = Complex(strands=[top], structure=".")
    startBot = Complex(strands=[bot], structure=".")
    
    
    # Turns Boltzmann sampling on for this complex and also does sampling more efficiently by sampling 'trials' states.
    if(myTrials > 0):
        setBoltzmann(startTop, myTrials)
        setBoltzmann(startBot, myTrials)
    
    # Stop when the exact full duplex is achieved.
    success_complex = Complex(strands=[top, bot], structure="(+)")
    stopSuccess = StopCondition(Options.STR_SUCCESS, [(success_complex, Options.exactMacrostate, 0)])
    
    # Declare the simulation unproductive if the strands become single-stranded again.
    failed_complex = Complex(strands=[top], structure=".")
    stopFailed = StopCondition(Options.STR_FAILURE, [(failed_complex, Options.dissocMacrostate, 0)])
    
    options.start_state = [startTop, startBot]

    # Point the options to the right objects
    if not doFirstPassage:
    
        options.stop_conditions = [stopSuccess, stopFailed]           
    
    else :

        options.stop_conditions = [stopSuccess]         
        
        
        
        
# Figure 2d has 3x12 = 36 rates plotted. 
# Input: 0 <= selector < 36  

# range 0 -11:  6 nt toehold
# range 12-23:  7 nt toehold
# range 24-36: 10 nt toehold

# order of mismatch position:
# perfect - 2 - 3 - 4 - 5 - 6 - 7 - 8 - 9 -10 - 12 - 14
   
def machinek2014(selector):
    
    # these are the sequences we need to build the dot-parens
    incumbent = ""
    target = ""
    invader = ""
    
    toeholdSelect = floor(selector / 12)
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

    if mismatchSelect == 0 | mismatchSelect == 2 | mismatchSelect == 12 | mismatchSelect == 14  :
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
        incumbent = "TGG TGT TTG TGG GTG TGG TGA GTT TGT GAA GGT"
        target = "CCC TCC ACA TAC CTT CAC AAA CTC ACC"
        invader = "GGT GAG TTT GTC AAG GTA TGT GG"
        
    if mismatchSelect == 7:
        incumbent = "TGG TGT TTG TGG GTG TGG TGA GTT TTG AGA GGT"
        target = "CCC TCC ACA TAC CTC TCA AAA CTC ACC"
        invader = "GGT GAG TTT TCA GAG GTA TGT GG"

    if mismatchSelect == 8:
        incumbent = "TGG TGT TTG TGG GTG TGG TGA GTT TGA TGA GGT"
        target = "CCC TCC ACA TAC CTC ATC AAA CTC ACC"
        invader = "GGT GAG TTT CAT GAG GTA TGT GG"
        
    if mismatchSelect == 9:
        incumbent = "TGG TGT TTG TGG GTG TGG TGAG TT GAT TGA GGT"
        target = "CCC TCC ACA TAC CTC AAT CAA CTC ACC"
        invader = "GGT GAG TTC ATT GAG GTA TGT GG"
        
    if mismatchSelect == 10:
        incumbent = "TGG TGT TTG TGG GTG TGG TGA GTG ATT TGA GGT"
        target = "CCC TCC ACA TAC CTC AAA TCA CTC ACC"
        invader = "GGT GAG TCA TTT GAG GTA TGT GG"
       

    invader = invader + toeholdSeq
    
 