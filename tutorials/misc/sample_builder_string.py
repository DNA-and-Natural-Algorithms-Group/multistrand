import time 

from multistrand.experiment import standardOptions, hybridization, hairpinclosing, dissociation, seesaw_gate_fuel_leak, threewayDisplacement
from multistrand.concurrent import FirstPassageRate
from multistrand.builder import Builder, BuilderRate, hybridizationString, dissociationString, threewaybmString
from multistrand.options import Options, Literals
from multistrand.system import SimSystem
from multistrand.objects import StopCondition

str_association = "association"
str_hairpin_closing = "hairpinclosing"
str_threeway_strand_displacement = "threewaystranddisplacement"
str_dissociation = "dissociation"

test3mer = "TTT"
test6mer = "TTGGTG"
test8mer = "TTGGTGAT"
test10mer = "TTGGTGATCC"
test15mer = "AGATTAGCAGGTTTC"
test20mer = "AGATTAGCAGGTTTCCCACC"

sumTime = 0.0


def getOptions(arguments):
     
    o = standardOptions()
    o.simulation_mode = Literals.trajectory
    o.num_simulations = 80

    o.temperature = 30.0
    o.simulation_time = 0.0000001
    
    endComplex = arguments[0]
    
    stopSuccess = StopCondition(Literals.success, [(endComplex, Literals.exact_macrostate, 0)])
    o.stop_conditions = [stopSuccess]
    
    return o
  

def getString(arguments): 
    
    toggle = arguments[0]

    if toggle == str_association:
        return hybridizationString(arguments[1])    
        
    elif toggle == str_threeway_strand_displacement:
        return threewaybmString("TTT", arguments[1], "")
    
    elif toggle == str_dissociation:
        return dissociationString(arguments[1])


def genAndPrint(numOfPaths, toggle):
    
    global sumTime
    
    print("Building the statespace from traces for reaction: " + toggle) 
    
    startStates = getString([toggle, test15mer])
    endState = startStates[-1]
    
    myBuilder = Builder(getOptions, [endState[0]])
    
    ''' setting the precision to just 2 states will ensure the builder stops after a single iteration. '''
    startTime = time.time()
#     Builder.verbosity = False
    myBuilder.genAndSavePathsFromString(startStates[:(len(startStates) - 1)])
#     Builder.verbosity = True
    myBuilder.fattenStateSpace()
    
    buildTime = time.time() - startTime
    print("Build time was %.2f s" % buildTime)
    sumTime += buildTime
    
    print(myBuilder)
    
    builderRate = BuilderRate(myBuilder) 
     
    if not (toggle == str_dissociation or toggle == str_threeway_strand_displacement):
        compTime = builderRate.averageTimeFromInitial(bimolecular=False)
        print("\nRate = %.2f /s  \n" % (1.0 / compTime))
    else:
        compTime = builderRate.averageTimeFromInitial(bimolecular=True)
        print("\nRate = %.2f /M/s  \n" % (1.0 / compTime))
     
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
       
    numOfPaths = 40
    genAndPrint(numOfPaths, str_association)
#     genAndPrint(numOfPaths, str_hairpin_closing)
#     genAndPrint(numOfPaths, str_dissociation)

#     numOfPaths = 5000        # put this too low, and there won't be any final states found
#     genAndPrint(numOfPaths, str_threeway_strand_displacement)
    
    print("Overal construction time was %.2f seconds" % sumTime)
    
