from __future__ import print_function

from multistrand.experiment import standardOptions, hybridization, hairpinclosing, dissociation, seesaw_gate_fuel_leak, threewayDisplacement
from multistrand.concurrent import FirstPassageRate
from multistrand.builder import Builder, BuilderRate
from multistrand.options import Options, Literals
from multistrand.system import SimSystem

str_association = "association"
str_hairpin_closing = "hairpinclosing"
str_threeway_strand_displacement = "threewaystranddisplacement"
str_dissociation = "dissociation"

test3mer = "TTT"
test6mer = "TTGGTG"
test8mer = "TTGGTGAT"
test10mer = "TTGGTGATCC"
test20mer = "AGATTAGCAGGTTTCCCACC"
  

def doReaction(arguments): 
    
    stdOptions = standardOptions()
    
    stdOptions.output_interval = 1
    stdOptions.simulation_mode = Literals.trajectory
    stdOptions.join_concentration = 0.001
    stdOptions.simulation_time = 400.0

    stdOptions.num_simulations = arguments[0]
    stdOptions.temperature = arguments[2]

    toggle = arguments[3]

    if toggle == str_association:
        hybridization(stdOptions, arguments[1], myTrials=arguments[0], doFirstPassage=True)    

    
    elif toggle == str_hairpin_closing:
        hairpinclosing(stdOptions, arguments[1], "TTTTT")
        
    elif toggle == str_threeway_strand_displacement:
        threewayDisplacement(stdOptions, "ACTACG", "AACATGAAGTA", myTrials=arguments[0])

    elif toggle == str_dissociation:
        dissociation(stdOptions, arguments[1], myTrials=arguments[0])
    
    """ Make sure that multistrand IDs equal, everytime.
        This is required for correctly combining multiple runs.  
        The Multistrand behavior makes sense, 
        but not for this specific use case.                            """
    
    i = 65
    for complex in stdOptions.start_state:
        for strand in complex.strand_list:
            strand.id = i
            i += 1

    stdOptions.DNA23Metropolis()
    
    return stdOptions


def genAndPrint(numOfPaths, toggle):
  
    simulation_temperature = 50.0
    
    print("Building the statespace from traces for reaction: " + toggle + " at T = %.1f C " % simulation_temperature)
    
    myBuilder = Builder(doReaction, [numOfPaths, test8mer, simulation_temperature, toggle])
    myBuilder.genAndSavePathsFile()
    myBuilder.genAndSavePathsFile()
    print(myBuilder)
    
    builderRate = BuilderRate(myBuilder) 
     
    if not (toggle == str_dissociation or toggle == str_threeway_strand_displacement):
        time = builderRate.averageTimeFromInitial(bimolecular=False)
        print("\nRate = %.2f /s  \n" % (1.0 / time))
    else:
        time = builderRate.averageTimeFromInitial(bimolecular=True)
        print("\nRate = %.2f /M/s  \n" % (1.0 / time))
     
    """ Fields you could look at: """
    if False:
        print(builderRate.stateIndex)
        print(builderRate.rate_matrix_csr.toarray())
       
        print("Initial States")
        print(builderRate.initial_states)
       
        print("Final States")
        print(builderRate.final_states)
       
        print(str(builderRate))
       
        times = builderRate.averageTime()
        print(str(times))
    

  

# # The actual main method
if __name__ == '__main__':
       
    numOfPaths = 30
    
    genAndPrint(numOfPaths, str_association)
    genAndPrint(numOfPaths, str_hairpin_closing)
    genAndPrint(numOfPaths, str_dissociation)

    numOfPaths = 5000        # put this too low, and there won't be any final states found
    genAndPrint(numOfPaths, str_threeway_strand_displacement)
    






    
#    """ See-saw gate code. In case it is needed. """
#     elif toggle == 5:
#         stdOptions.simulation_mode = Literals.trajectory
#         stdOptions.temperature = arguments[3]
#         stdOptions.join_concentration = arguments[2]
#         
#         # Sequence design from Thubagere, 2017
#         CL_LONG_S18 = "TCTTCTAACAT"
#         CL_LONG_S5 = "CCACCAAACTT"
#         CL_LONG_S29 = "CCAATACTCCT"
#         CL_LONG_S44 = "AAACTCTCTCT"
#         CL_LONG_SEQT = "TCT"
#         CLAMP_SEQ = "CA"
#         
#         CL_LONG_GATE_A_SEQ = [CL_LONG_S44, CL_LONG_S18,
#                               CL_LONG_S5, CL_LONG_S29, CL_LONG_SEQT, CLAMP_SEQ]
#         
#         gateA = ClampedSeesawGate(*CL_LONG_GATE_A_SEQ, sameID=True)
#         
#         seesaw_gate_fuel_leak(stdOptions, gateA, trials=arguments[0], supersample=1, doFirstPassage=True)
#         
#         stdOptions.simulation_mode = Literals.trajectory
#     
