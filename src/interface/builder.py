# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
The builder code has two classes: Builder and BuilderRate. A Builder will accept
references to a function that returns a Multistrand Options object. The sampled
states and transitions will be collected.

The BuilderRate accepts a Builder object and will compute the mean first passage
time from the starting states until hitting an end state.
"""

import time, copy, os, sys

from multistrand.system import SimSystem
from multistrand.utils import uniqueStateID, seqComplement, GAS_CONSTANT
from multistrand.options import Literals
from multistrand.experiment import standardOptions, makeComplex

from scipy.sparse import csr_matrix, coo_matrix
from scipy.sparse.linalg import spsolve, bicg, bicgstab, cg, gmres, lgmres

import numpy as np

floatT = np.float64


def hybridizationString(seq):
    """
    Constructs a surface of N * ( N -1 ) / 2 points along the reaction frontier.
    This method returns a list of starting states that should be used as initial
    states for the string method.
    A starting state is a list of complexes.
    """
    cutoff = 0.75
    ids = [65, 66]

    N = len(seq)
    output = []

    """ start with the seperated, zero-bp state """
    dotparen1 = "."*N + "+" + "."*N

    complex0 = makeComplex([seq], "."*N, [ids[0]])
    complex1 = makeComplex([seqComplement(seq)], "."*N, [ids[1]])
    seperated = [complex0, complex1]

    output.append(seperated)

    """ 
        bps is the number of basepairs between the strands
        offset is the offset of where the pairing starts
    """

    dotparen0 = [seq, seqComplement(seq)]

    for bps in range(1, N + 1):
        for offset in range(N):
            dotparen1 = list("."*N + "+" + "."*N)
            if (bps + offset < (N + 1)) and bps < np.ceil(cutoff * N):
                #             if (bps + offset < (N + 1)):
                for i in range(bps):
                    dotparen1[offset + i] = "("
                    dotparen1[ 2 * N - offset - i ] = ")"

                dotparen1 = "".join(dotparen1)
                complex = makeComplex(dotparen0, dotparen1, ids)
                output.append([complex])

    """ We always require the final state to be the success state"""
    complex0 = makeComplex(dotparen0, "("*N + "+" + ")"*N, ids)
    output.append([complex0])
    return output


def dissociationString(seq):
    ''' Exactly like the association, but swap the initial and final states '''
    
    myList = hybridizationString(seq)
    
    first = myList[0]
    last = myList[-1]
    
    myList[0] = last
    myList[-1] = first

    return myList

    
''' returns a list of strings pairs that represent the hybridization steps 
    for toeholds / domains during pathway elaboration method    '''
def weave(M):
    output = []
    for bps in range(1, M + 1):
        for offset in range(M):
            leftStr = list("."*M)
            rightStr = list("."*M)
            if (bps + offset < (M + 1)):
                for i in range(bps):
                    leftStr[offset + i] = "("
                    rightStr[ M - offset - i - 1 ] = ")"
    
                leftStr = "".join(leftStr)
                rightStr = "".join(rightStr)
                output.append([leftStr, rightStr])
    return output


def threewaybmString(lefttoe, displace, righttoe):

    ''' toehold switch around in the substrate '''
    N = len(displace)
    rT = len(lefttoe)
    lT = len(righttoe) 
    
    output = list()
    
    invaderSq = lefttoe + displace + righttoe
    incumbentSq = displace
    substrateSq = seqComplement(invaderSq)
    
    invaderID = 65
    incumbentID = 66
    substrateID = 67    
    
    """ start with the separated, zero-bp state """
    complex0 = makeComplex([invaderSq], "."*(lT + N + rT), [invaderID])
    complex1 = makeComplex([incumbentSq, substrateSq], "(" * N + "+" + "."* lT + ")" * N + "." *rT, [incumbentID, substrateID])
    seperated = [complex0, complex1]
    
    output.append(seperated)

    ''' left invasion toehold '''
    seqsL = [invaderSq, substrateSq, incumbentSq]
    idsL = [invaderID, substrateID, incumbentID]
    weaving = weave(lT)
    
    for pair in weaving:
        dotparen = "."*(N + rT) + pair[0] + "+" + pair[1] + "("*N + "."*rT + "+" + ")"*N 
        output.append([makeComplex(seqsL, dotparen, idsL)])
        
    if lT == 0:
        dotparen = "."*(N + rT) + "" + "+" + "" + "("*N + "."*rT + "+" + ")"*N
        
    ''' the dotparen of the weave is the fully hybridized toehold.
        Toggle the basepairs one step at a time.
    '''
    for invasion in range(N - 1):
        parenList = list(dotparen)
        parenList[rT + N - 1 - invasion] = "("
        parenList[lT + N + rT + 1 + lT + invasion ] = ")"
        parenList[-invasion - 1] = "."
        dotparen = "".join(parenList)
        output.append([makeComplex(seqsL, dotparen, idsL)])
        
    ''' right invasion toehold '''
    seqsR = [invaderSq, incumbentSq, substrateSq]
    idsR = [invaderID, incumbentID, substrateID]
    weaving = weave(rT)
    
    for pair in weaving:
        dotparen = pair[0] + "."*(N + lT) + "+" + "("*N + "+" + "."*lT + ")"*N + pair[1] 
        output.append([makeComplex(seqsR, dotparen, idsR)])

    if rT == 0:
        dotparen = "" + "."*(N + lT) + "+" + "("*N + "+" + "."*lT + ")"*N + "" 

    ''' 
        Toggle the basepairs one step at a time.
    '''
    for invasion in range(N - 1):
        parenList = list(dotparen)
        parenList[rT + invasion] = "("
        parenList[lT + N + rT + 1 + invasion ] = "."
        dotparen = "".join(parenList)
        output.append([makeComplex(seqsR, dotparen, idsR)])
        
    ''' Do not forget to set the final state.
        This is just the displaced strand floating freely.
    '''
    output.append([makeComplex([incumbentSq], "."*N, [incumbentID])])
    return output


class ConvergeCrit:

    period = 4  # average out over past X increases.

    def __init__(self):
        self.maximumIterations = 1000
        self.minimumStateIncrement = 4
        
        self.currIteration = 0
        self.currStates = -99

        self.precision = 0.05
        self.array = [-99.0] * self.period

    def converged(self, rateIn=None, statespace_size=None):  # also saves the rateIn
        if self.currIteration > self.maximumIterations:
            return True

        if self.precision < 1.0 :
            averageL = (sum(self.array) / self.period) * (1.0 - self.precision)
            averageH = (sum(self.array) / self.period) * (1.0 + self.precision)

            conv1 = rateIn > averageL
            conv2 = rateIn < averageH

            self.array[self.currIteration % self.period] = rateIn
            self.currIteration += 1;
            return (conv1 and conv2)

        else:
            if statespace_size - self.currStates < self.minimumStateIncrement:
                return True
                        
            self.currStates = statespace_size;
            return statespace_size >= self.precision

    def __str__(self):
        return  str(self.array)


class transitiontype:

    unimolecular = "uni"
    bimolecularIn = "bi-in"
    bimolecularOut = "bi-out"

    array = [unimolecular, bimolecularIn, bimolecularOut ]


class localtype:

    end = "End"
    loop = "Loop"
    stack = "Stack"
    stackstack = "StackStack"
    loopend = "LoopEnd"
    stackend = "StackEnd"
    stackloop = "StackLoop"

    array = [end, loop, stack, stackstack, loopend, stackend, stackloop]


class Energy:

    dH = 0.0;
    dS = 0.0;

    def __init__(self, dG, dH, temp, concentration, n_complexes, n_strands):
        assert(200 < temp < 400)

        RT = GAS_CONSTANT * floatT(temp);
        dG_volume = RT * (n_strands - n_complexes) * np.log(1.0 / concentration)

        self.dH = floatT(dH)
        self.dS = floatT(-(dG - dG_volume - dH) / temp)
        
    def dG(self, temp):
        return self.dH - floatT(temp) * self.dS

    def __str__(self):
        return "dH= " + str(self.dH) + " dS= " + str(self.dS)

    def __eq__(self, that):
        if isinstance(that, Energy):
            return self.dH == that.dH and self.dS == that.dS
        else:
            raise Exception("Not an acceptable argument", "__eq__", that)


def codeToDesc(code):

    output = []
    primes = [3, 5, 7, 11, 13, 17, 19]
    primesq = [9, 25, 49, 121, 169, 289, 361]

    for i in range(len(primes)):

        if code % primes[i] == 0:
            output.append(localtype.array[i])

            if code % primesq[i] == 0:  # code is a square number
                output.append(localtype.array[i])
                return output

    return output


class InitCountFlux:

    def __init__(self):
        self.count = 0;
        self.flux = 0.0;

    def __add__(self, x):
        self.count + x.count

    def __str__(self):
        return "Count: " + str(self.count) + "  flux: " + str(self.flux)

    def  __repr__(self):
        return str(self)


class Builder:

    verbosity = False

    # input function returns the multistrand options object for which to build the statespace
    def __init__(self, inputFunction, arguments):
        # initial argument has to be the number of trials
        self.optionsFunction = inputFunction
        self.optionsArgs = copy.deepcopy(arguments)

        self.doMultiprocessing = False
        self.printTimer = True
        self.numOfThreads = 8

        self.protoSpace = dict()  # key: states. Value: Energy
        self.protoTransitions = dict()  # key: transitions. Value: ArrheniusType (negative if it is a bimolecular transition)
        self.protoInitialStates = dict()  # key: states. Value: a InitCountFlux object that tells how many times the state has been the initial state and the join flux (rate)
        self.protoFinalStates = dict()  # key: states: Value: the result of this final state can be SUCCES or FAILURE
        self.protoSequences = dict()  # key: name of strand. Value: sequence

        self.firstStepMode = True
        self.startTime = time.time()
        
        self.mergingCounter = 0  # counts how many transitions from the merging have been found

        # save a copy for later processing -- note this copy will not have results attached to it,
        # so it won't have a large memory footprint.
        self.options = self.optionsFunction(self.optionsArgs)

        self.the_dir = "p_statespace/"

    def __str__(self):
        output = "states / transitions / initS / finalS / merged     \n "
        output += str(len(self.protoSpace)) 
        output += "   -    " + str(len(self.protoTransitions))
        output += "   -    " + str(len(self.protoInitialStates)) 
        output += "   -    " + str(len(self.protoFinalStates))
        output += "   -    " + str(self.mergingCounter) 
        return output

    def reset(self):
        self.protoSpace.clear()
        self.protoTransitions.clear()
        self.protoInitialStates.clear()
        self.protoFinalStates.clear()
        self.protoSequences.clear()

    def printOverlap(self, other):
        N = len(self.protoSpace)
        overlap = 0
        for state in self.protoSpace:
            if state in other.protoSpace:
                overlap += 1

        print("The overlap is " + str(100.0 * overlap / N) + " percent. ")

    def mergeSet(self, this, that):
        for key, val in that.items():
            if not key in this:
                this[key] = val

    def mergeBuilder(self, other):
        self.mergeSet(self.protoSpace, other.protoSpace)
        self.mergeSet(self.protoTransitions, other.protoTransitions)
        self.mergeSet(self.protoInitialStates, other.protoInitialStates)
        self.mergeSet(self.protoFinalStates, other.protoFinalStates)
        self.mergeSet(self.protoSequences, other.protoSequences)

    ''' Merges if both source and the target exist '''
    def transitionMerge(self, other):
        for key, value in other.protoTransitions.items():
            sFrom = key[0]
            sTo = key[1]
            
            if sFrom in self.protoSpace and sTo in self.protoSpace:
                if not key in self.protoTransitions:
                    self.protoTransitions[key] = value
                    self.mergingCounter += 1

    def parseState(self, line, simulatedTemperature, simulatedConc):
        mywords = line.split()

        n_complexes = int(mywords[0])
        n_strands = 0

        ids = []
        sequences = []
        structs = []

        for i in range(n_complexes):
            ids.append(mywords[1 + i])
            sequences.append(mywords[1 + n_complexes + i])
            structs.append(mywords[1 + 2 * n_complexes + i])
            n_strands += len(mywords[1 + 2 * n_complexes + i].split('+'))

            uniqueID2 = uniqueStateID(ids, structs)
            uniqueID = ""
            uniqueID =  tuple(uniqueID2)

        dG = float(mywords[1 + 3 * n_complexes])
        dH = float(mywords[1 + 3 * n_complexes + 1])

        energyvals = Energy(dG, dH, simulatedTemperature, simulatedConc, n_complexes, n_strands)

        return uniqueID, energyvals, (sequences, ids, structs)

    """ Runs genAndSavePathsFile until convergence is reached"""
    def genUntilConvergence(self, precision):
        crit = ConvergeCrit()
        crit.precision = precision

        currTime = -1.0

        while not crit.converged(currTime, len(self.protoSpace)):
            self.genAndSavePathsFile()

            if self.verbosity:
                print("Size     = %i " % len(self.protoSpace))
            if precision < 1.0:
                    builderRate = BuilderRate(self)
                    currTime = builderRate.averageTimeFromInitial()

        self.fattenStateSpace()
        if self.verbosity:
            print("Size     = %i " % len(self.protoSpace))

    """ Runs genAndSavePathsFile until convergence is reached,
        given a list of initial states"""
    def genUntilConvergenceWithInitialState(self, precision, initialStates, printMeanTime=False):
        crit = ConvergeCrit()
        crit.precision = precision

        currTime = -1.0
        while not crit.converged(currTime, len(self.protoSpace)) :
            self.genAndSavePathsFromString(initialStates, printMeanTime=printMeanTime)
            
            if precision < 1.0:
                builderRate = BuilderRate(self)
                currTime = builderRate.averageTimeFromInitial()
            if printMeanTime:
                print("Mean first passage time = %.2E" % currTime)

        self.fattenStateSpace()
        if self.verbosity:
            print("Size     = %i " % len(self.protoSpace))

    ''' A single iteration of the pathway elaboration method '''
    def genAndSavePathsFromString(self, pathway, printMeanTime=False):
        startTime = time.time()

        """ Only the first state will count towards the set of initial states """
        ignoreInitial = False
        for state in pathway:
            otherBuilder = Builder(self.optionsFunction, self.optionsArgs)
            otherBuilder.genAndSavePathsFile(supplyInitialState=state, ignoreInitialState=ignoreInitial)
            ignoreInitial = True

            self.mergeBuilder(otherBuilder)
            del otherBuilder

        if self.verbosity or printMeanTime:
            print("Size     = %i    ---  bytesize = %i " % (len(self.protoSpace), sys.getsizeof (self.protoSpace)))
            print("Size T   = %i    ---  bytesize = %i " % (len(self.protoTransitions), sys.getsizeof (self.protoTransitions)))
            print("Time = %f" % (time.time() - startTime))

    '''
    Generates all transitions between states in the statespaces and adds missing transitions
    '''
    def fattenStateSpace(self):
        """ Load up the energy model with a near zero-length simulation """
        myOptions = self.optionsFunction(copy.deepcopy(self.optionsArgs))
        myOptions.simulation_mode = Literals.trajectory
        myOptions.activestatespace = False
        myOptions.simulation_time = 0.0000000001
        myOptions.start_state = hybridizationString("AAA")[0] #pathway[0]
#          
        s = SimSystem(myOptions)
        s.start()  
        """ Done setting the energy model        """

        ogVerb = Builder.verbosity
        Builder.verbosity = False
        counter = 0
        
        def inspectionSim(inputs):
            o1 = standardOptions()
            """ Disable regenerating the energymodel """
            o1.reuse_energymodel = True 
            o1.rate_method = self.options.rate_method
            o1.start_state = inputs[0]
            return o1
        
        for key, _ in self.protoSpace.items():
            if ogVerb and ((counter % 100) == 0):
                print("Searching for missing transitions. Progress " + str(counter) + " / " + str(len(self.protoSpace)))

            (seqs, ids, structs) = self.protoSequences[key]
            myState = []
            
            for seq, id, struct in zip(seqs, ids, structs):
                seqs = seq.split('+')
                ids = id.split(',')
                myC = makeComplex(seqs, struct, ids)
                myState.append(myC)
            
            ''' post: myState is the state we want to explore transitions for. '''
    
            myB = Builder(inspectionSim, [myState])
            myB.genAndSavePathsFile(inspecting=True)
            
            self.transitionMerge(myB)
            counter += 1
        
        Builder.verbosity = ogVerb

    """
    Computes the mean first pasasage times, 
    then selects states that are delta-close 
    to the set of final states. 
    Those states are then added to the set of final states. 
    This reduces the size of the matrix that is constructed.
    """
    def deltaPruning(self, delta=0.01, printCount=False):
        builderRate = BuilderRate(self)
        firstpassagetimes = builderRate.averageTime()

        sumTime = 0.0
        sumStart = 0.0

        if printCount:
            beforeN = len(self.protoFinalStates)

        for state in builderRate.initial_states:

            stateindex = builderRate.stateIndex[state]
            sumTime += builderRate.initial_states[state].count * firstpassagetimes[stateindex]
            sumStart += builderRate.initial_states[state].count

        averagedMFPT = sumTime / sumStart

        """Now add states that are delta close to the set of final states"""

        for state in builderRate.statespace:
            if state not in self.protoFinalStates:
                stateindex = builderRate.stateIndex[state]
                if firstpassagetimes[stateindex] < delta * averagedMFPT:
                    self.protoFinalStates[state] = Literals.success

        if printCount:
            print("Number of final states was %i but now is %i" % (beforeN, len(self.protoFinalStates)))

    """
    supplyInitialState : a Complex that serves as initial state for that simulation
    ignoreIntiialState: The initial state is not added to the set of initial states
    """
    def genAndSavePathsFile(self, ignoreInitialState=False, supplyInitialState=None, inspecting=False):
        self.startTime = time.time()

        space = dict()
        transitions = dict()
        initStates = dict()
        finalStates = dict()
        sequences = dict()

        # the first argument is always the number of paths
        inputArgs = copy.deepcopy(self.optionsArgs)

        def runPaths(optionsF, optionsArgs, space, transitions, initStates, finalStates, sequences):
            myOptions = optionsF(optionsArgs)
            myOptions.activestatespace = True
            myOptions.output_interval = 1

            if not supplyInitialState == None:
                myOptions.start_state = supplyInitialState

            """ Set longer searching time for the initial state. """
            if not ignoreInitialState:
                myOptions.simulation_time = myOptions.simulation_time * 10.0

            simTime = time.time()
            s = SimSystem(myOptions)
                
            ''' a specialized routine that ends after taking all transitions'''
            if inspecting == True:
                s.localTransitions() 
            else:
                s.start()  # after this line, the computation is finished.

            if self.verbosity:
                print("Multistrand simulation is now done,      time = %.2f" % (time.time() - simTime))

            """ load the space """
            myFile = open(self.the_dir + str(myOptions.interface.current_seed) + "/protospace.txt", "r")
            for line in myFile:
                uniqueID, energyvals, seqs = self.parseState(line, myOptions._temperature_kelvin, myOptions.join_concentration)
                
                if not uniqueID in sequences:
                    sequences[uniqueID] = seqs
                if not uniqueID in space:
                    space[uniqueID] = energyvals
                elif not space[uniqueID] == energyvals:
                    print("My hashmap contains " + str(uniqueID) + " with Energy " + str(space[uniqueID]) + " but found: " + str(energyvals))
                    print("Line = " + line)

            """ load the transitions """
            myFile = open(self.the_dir + str(myOptions.interface.current_seed) + "/prototransitions.txt", "r")
            index = 0
            go_on = True
            myLines = []
            for line in myFile:
                myLines.append(line)

            while go_on:
                line1 = myLines[index];
                line2 = myLines[index + 1];
                line3 = myLines[index + 2];
                
                index = index + 4  # note the whitespace

                go_on = len(myLines) > index

                uID1, ev1, seq1 = self.parseState(line2, myOptions._temperature_kelvin, myOptions.join_concentration)
                uID2, ev2, seq2 = self.parseState(line3, myOptions._temperature_kelvin, myOptions.join_concentration)

                transitionPair = (uID1, uID2)

                if not transitionPair in transitions:

                    transitionList = list()

                    n_complex1 = int(line2.split()[0])
                    n_complex2 = int(line3.split()[0])

                    if n_complex1 == n_complex2:
                        transitionList.append(transitiontype.unimolecular)

                    if n_complex1 > n_complex2:
                        transitionList.append(transitiontype.bimolecularIn)

                    if n_complex2 > n_complex1:
                        transitionList.append(transitiontype.bimolecularOut)

                    if myOptions.rate_method == Literals.arrhenius:
                        # decode the transition and add it
                        transitionList.extend(codeToDesc(int(float(line1))))

                    transitions[transitionPair] = transitionList

            """ load the initial states """
            myFile = open(self.the_dir + str(myOptions.interface.current_seed) + "/protoinitialstates.txt", "r")
            myLines = []
            for line in myFile:
                myLines.append(line)
            index = 0
            go_on = True
            if len(myLines) == 0:
                print("No initial states found!")

            while go_on:
                line1 = myLines[index];
                line2 = myLines[index + 1];

                index = index + 2  # note the whitespace
                go_on = len(myLines) > index

                uID1, ev1, seq1 = self.parseState(line2, myOptions._temperature_kelvin, myOptions.join_concentration)
                count = int(line1.split()[0])

                if not uID1 in initStates:
                    newEntry = InitCountFlux()
                    newEntry.count = count
                    newEntry.flux = 777777  # arrType is the flux, and is unique to the initial state
                    initStates[uID1] = newEntry

            """ load the final states """
            myFile = open(self.the_dir + str(myOptions.interface.current_seed) + "/protofinalstates.txt", "r")
            myLines = []
            for line in myFile:
                myLines.append(line)
            index = 0
            go_on = True
            if len(myLines) == 0:
                #                 raise ValueError("No succesful final states found -- mean first passage time would be infinite ")
                go_on = False

            while go_on:
                line1 = myLines[index];
                line2 = myLines[index + 1];
                index = index + 2

                go_on = len(myLines) > (index + 1)

                uID1, ev1, seq1 = self.parseState(line1, myOptions._temperature_kelvin, myOptions.join_concentration)
                tag = line2.split()[0]

                if not uID1 in finalStates:
                    finalStates[uID1] = tag

            """ Now delete the files as they can get quite large """
            os.remove(self.the_dir + str(myOptions.interface.current_seed) + "/protospace.txt")
            os.remove(self.the_dir + str(myOptions.interface.current_seed) + "/prototransitions.txt")
            os.remove(self.the_dir + str(myOptions.interface.current_seed) + "/protoinitialstates.txt")
            os.remove(self.the_dir + str(myOptions.interface.current_seed) + "/protofinalstates.txt")
            os.rmdir(self.the_dir + str(myOptions.interface.current_seed)) 

        runPaths(self.optionsFunction, inputArgs, space, transitions, initStates, finalStates, sequences)

        # do not forget to merge the objects back
        self.mergeSet(self.protoSpace, space)
        self.mergeSet(self.protoTransitions, transitions)
        self.mergeSet(self.protoFinalStates, finalStates)
        self.mergeSet(self.protoSequences, sequences)

        if not ignoreInitialState:
            for key, val in initStates.items():
                if not key in self.protoInitialStates:
                    self.protoInitialStates[key] = val;
                else:
                    self.protoInitialStates[key].count += val.count


# FD: This class is in progress.
# this class takes a builder object and computes the average time between
# starting in an initial state and reaching a final state.
class BuilderRate:

    solveToggle = 2

    # input function returns the multistrand options object for which to build the statespace
    def __init__(self, builderIn):

        self.build = builderIn
        self.rateLimit = 1e-5

        if len(self.build.protoFinalStates) == 0 :
            raise ValueError('No final states found.')

        self.processStates()  # prunes statespace and creates objects that can be used to create the rate matrx
        self.setMatrix()  # generates the matrix for the current temperature

    """ Generates the state space by traversing from the final states """
    def findConnectedStates(self, final_states, neighborsdict):

        new_statespace = set()
        unexplored = set()

        for state in final_states:
            unexplored.add(state)

        while unexplored:  # empty dicts evaluate to false in python
            new_unexplored = set()
            for state in unexplored:
                new_statespace.add(state)

                # Make sure all neighbors are explored whenever we add a state.
                neighbors = neighborsdict[state]
                for s in neighbors :
                    if s not in new_statespace :
                        new_unexplored.add(s)
            unexplored = new_unexplored

        return new_statespace

    """ generates a list of neighbors for each state -- each transition goes both ways """
    def genNeighborsTransitive(self, transitions, statespace):

        output = dict()  # key: state, value: a set of states that are neighboring

        # each state starts out with no neighbors
        for state in statespace:
            output[state] = list()

        for transition in transitions:
            output[transition[0]].append(transition[1])
            output[transition[1]].append(transition[0])

        return output

    """ generates a list of neighbors so that for each transition only one direction is added
        This is important when we build the matrix - so that we do not doubly add transitions """
    def genNeighbors(self, transitions, statespace):

        output = dict()  # key: state, value: a set of states that are neighboring

        # each state starts out with no neighbors
        for state in statespace:
            output[state] = list()

        for transition in transitions:
            # ensure the state is connected and add only one direction of a transition
            if (transition[0] in statespace) and not (transition[0] in output[transition[1]]) :
                output[transition[0]].append(transition[1])

        return output

    def processStates(self):

        self.statespace = set()
        self.initial_states = dict()  # key: state, value: number of times started in this state
        self.final_states = set()

        # Storing the final states
        for state in self.build.protoFinalStates:
            if self.build.protoFinalStates[state] == Literals.success:
                self.final_states.add(state)

        if len(self.final_states) == 0:
            raise ValueError("No final states found!")

        startT = time.time()
        # generate a list of neighbors for each state  -- output is placed in self.neighbors
        self.neighbors = self.genNeighborsTransitive(self.build.protoTransitions, self.build.protoSpace)

        startT = time.time()
        # now prune the statespace to only include states that can reach the final state
        self.statespace = self.findConnectedStates(self.final_states, self.neighbors)

        startT = time.time()
        # Now re-generated the neighbors, but only for states in the statespace -- and only "forward" transitions instead of both ways.
        self.neighbors = self.genNeighbors(self.build.protoTransitions, self.build.protoSpace)

        # Update the initial states to only include those connected states
        for state in self.build.protoInitialStates:
            if state in self.statespace:
                self.initial_states[state] = self.build.protoInitialStates[state]

    """
    Uses the Arrhenius kinetic model to calculate transition rates. 
    Returns the transition rate from state1 to state2 
    and then also the reverse rate, from state2 to state1 
    """
    def halfcontext_parameter(self, localContext):

        if localContext == localtype.stack:
            return self.build.options.lnAStack , self.build.options.EStack
        elif localContext == localtype.loop :
            return self.build.options.lnALoop, self.build.options.ELoop
        elif localContext == localtype.end:
            return self.build.options.lnAEnd, self.build.options.EEnd
        elif localContext == localtype.stackloop:
            return self.build.options.lnAStackLoop, self.build.options.EStackLoop
        elif localContext == localtype.stackend :
            return self.build.options.lnAStackEnd, self.build.options.EStackEnd
        elif localContext == localtype.loopend :
            return self.build.options.lnALoopEnd, self.build.options.ELoopEnd
        elif localContext == localtype.stackstack:
            return self.build.options.lnAStackStack, self.build.options.EStackStack
        else:
            raise ValueError('The transition code name is unexpected ')

    """Use the builder options object to determine which rates to compute """
    def get_rate(self, state1, state2):

        transitionlist = self.build.protoTransitions[(state1, state2)]
        if self.build.options.rate_method == Literals.arrhenius:
            return self.arrhenius_rate(state1, state2, transitionlist)
        else:
            return self.metropolis_rate(state1, state2, transitionlist)

    """    Returns the transition rate from state1 to state2 and then also the reverse rate, from state2 to state1 """
    def metropolis_rate(self, state1, state2, transitionlist):

        myT = self.build.options._temperature_kelvin
        RT = GAS_CONSTANT * myT

        dG1 = self.build.protoSpace[state1].dG(myT)
        dG2 = self.build.protoSpace[state2].dG(myT)

        if transitionlist[0] == transitiontype.unimolecular:
            if dG1 > dG2 :  # state2 is more stable (negative), dG1 - dG2 is positive
                rate1 = self.build.options.unimolecular_scaling
                rate2 = self.build.options.unimolecular_scaling * np.exp(-(dG1 - dG2) / RT)
            else:  # state2 is less or equally stable (negative), dG1 - dG2 is negative
                rate1 = self.build.options.unimolecular_scaling * np.exp((dG1 - dG2) / RT)
                rate2 = self.build.options.unimolecular_scaling 
            return rate1, rate2

        else:  # bimolecular rate
            collisionRate = self.build.options.join_concentration * self.build.options.bimolecular_scaling
            if transitionlist[0] == transitiontype.bimolecularIn:
                outR = self.build.options.bimolecular_scaling * np.e ** (-(dG1 - dG2) / RT)
                return collisionRate, outR
            elif transitionlist[0] == transitiontype.bimolecularOut:
                outR = self.build.options.bimolecular_scaling * np.e ** ((dG1 - dG2) / RT)
                return outR, collisionRate
            else:
                raise ValueError('The transition code is unexpected ')

    def arrhenius_rate(self, state1, state2, transitionlist):

        lnA_left, E_left = self.halfcontext_parameter(transitionlist[1])
        lnA_right, E_right = self.halfcontext_parameter(transitionlist[2])

        lnA = lnA_left + lnA_right
        E = E_left + E_right

        bimolecular_scaling = self.build.options.bimolecular_scaling
        concentration = self.build.options.join_concentration

        myT = self.build.options._temperature_kelvin
        RT = GAS_CONSTANT * myT
        dG1 = self.build.protoSpace[state1].dG(myT)
        dG2 = self.build.protoSpace[state2].dG(myT)

        DeltaG = dG2 - dG1
        DeltaG2 = -DeltaG

        if transitionlist[0] == transitiontype.unimolecular:
            if DeltaG > 0.0:
                rate1 = np.e ** (lnA - (DeltaG + E) / RT)
                rate2 = np.e ** (lnA - E / RT)
            else:
                rate1 = np.e ** (lnA - E / RT)
                rate2 = np.e ** (lnA - (DeltaG2 + E) / RT)

            # if rate1 < self.rateLimit:
            #         rate1 = 0.0
            #
            # if rate2 < self.rateLimit:
            #     rate2 = 0.0

        elif transitionlist[0] == transitiontype.bimolecularIn:
            rate1 = bimolecular_scaling * concentration * (np.e ** (lnA - E / RT))
            rate2 = bimolecular_scaling * (np.e ** (lnA - (DeltaG2 + E) / RT))
        elif transitionlist[0] == transitiontype.bimolecularOut:
            rate1 = bimolecular_scaling * (np.e ** (lnA - (DeltaG + E) / RT))
            rate2 = bimolecular_scaling * concentration * (np.e ** (lnA - E / RT))
        else:
            raise ValueError('Exception in transition rate calculations.')
        return rate1, rate2

    """Set the rate matrix for this transition. If the target state is a final state, only subtract the outgoing rate from the diagonal. """
    def addTransition(self, state, neighbor, rate, rates, iArray, jArray, stateIndex):

        if not state in self.final_states:
            # now add the negative outgoing rate on the diagonal
            rates[stateIndex[state]] += -1.0 * rate
            if not neighbor in self.final_states:
                # set the transition
                rates.append(rate)
                iArray.append(stateIndex[state])
                jArray.append(stateIndex[neighbor])

    """ given the statespace, the initial states and final states, and the original builder object,
    build rate matrix for the current temperature. """
    def setMatrix(self):

        # give every state an explicit index
        self.stateIndex = dict()
        N = 0

        # the index works like this: rate[i] is in position iArray[i], jArray[j]
        # This is then interperted as an NxN matrix
        rates = list()
        iArray = list()
        jArray = list()

        for state in self.statespace:
            if not state in self.final_states:
                # give every state an explicit index
                self.stateIndex[state] = N

                # populate the diagonal so we can manipulate this as we go along
                rates.append(0.0)
                iArray.append(N)
                jArray.append(N)

                N += 1

        # post: N is the size of the statespace
        # post: rates[stateIndex[state]] is the diagonal for that state

        # first, set all transitions
        for state in self.statespace:
            for neighbor in self.neighbors[state]:
                myRate, revRate = self.get_rate(state, neighbor)
                if myRate < 0.0 or revRate < 0.0:
                    ValueError('Negative transition rate found.')

                # This handles either state being a final state (in which case, subtract from the non-final state diagonal,
                # but do not add the transition rate.
                if myRate > self.rateLimit:
                    self.addTransition(state, neighbor, floatT(myRate), rates, iArray, jArray, self.stateIndex)
                if revRate > self.rateLimit:
                    self.addTransition(neighbor, state, floatT(revRate), rates, iArray, jArray, self.stateIndex)

        # now actually create the matrix
        rate_matrix_coo = coo_matrix((rates, (iArray, jArray)), shape=(N, N) , dtype=floatT)

        self.rate_matrix_csr = csr_matrix(rate_matrix_coo, dtype=floatT)
        self.b = -1 * np.ones(N, dtype=floatT)

        #         # FD: pre-compute the matrix diagonal for preconditioning
        diagons = [ ((1.0 / x), i, j) for x, i, j in zip(rates, iArray, jArray) if i == j ]
        diagons0 = [x[0] for x in diagons]
        diagons1 = [x[1] for x in diagons]
        diagons2 = [x[2] for x in diagons]

        diagonal_matrix_coo = coo_matrix((diagons0, (diagons1, diagons2)), shape=(N, N), dtype=floatT)
        self.rate_matrix_inverse = csr_matrix(diagonal_matrix_coo, dtype=floatT)

        # save two counts for later interest
        self.n_states = N
        self.n_transitions = len(rates)

    """ 
        Computes the first passage times
    """
    def averageTime(self, x0=None, maxiter=None):

        startTime = time.time()

        if self.solveToggle == 1:
            firstpassagetimes, info = bicg(self.rate_matrix_csr, self.b, x0=x0, maxiter=maxiter)
        elif self.solveToggle == 2:
            firstpassagetimes, info = bicg(self.rate_matrix_csr, self.b, M=self.rate_matrix_inverse, x0=x0, maxiter=maxiter)
        elif self.solveToggle == 3:
            firstpassagetimes, info = bicgstab(self.rate_matrix_csr, self.b, x0=x0, maxiter=maxiter)
        elif self.solveToggle == 4:
            firstpassagetimes, info = bicgstab(self.rate_matrix_csr, self.b, M=self.rate_matrix_inverse, x0=x0, maxiter=maxiter)
        elif self.solveToggle == 5:
            firstpassagetimes, info = gmres(self.rate_matrix_csr, self.b, x0=x0, maxiter=maxiter)
        elif self.solveToggle == 6:
            firstpassagetimes, info = gmres(self.rate_matrix_csr, self.b, M=self.rate_matrix_inverse, x0=x0, maxiter=maxiter)
        elif self.solveToggle == 7:
            firstpassagetimes, info = cg(self.rate_matrix_csr, self.b, x0=x0, maxiter=maxiter)
        elif self.solveToggle == 8:
            firstpassagetimes, info = cg(self.rate_matrix_csr, self.b, M=self.rate_matrix_inverse, x0=x0, maxiter=maxiter)
        elif self.solveToggle == 9:
            firstpassagetimes, info = lgmres(self.rate_matrix_csr, self.b, x0=x0, maxiter=maxiter)
        elif self.solveToggle == 10:
            firstpassagetimes, info = lgmres(self.rate_matrix_csr, self.b, M=self.rate_matrix_inverse, x0=x0, maxiter=maxiter)
        else:
            self.setMatrix()
            firstpassagetimes = spsolve(self.rate_matrix_csr, self.b)

        self.matrixTime = time.time() - startTime
        return firstpassagetimes

    def averageTimeFromInitial(self, bimolecular=False, printMeanTime=False):

        times = self.averageTime()
        if self.build.verbosity or printMeanTime:
            print("Solving matrix took %.2f s" % self.matrixTime)
        mfpt = self.weightedPassageTime(times=times, bimolecular=bimolecular)
        return mfpt

    """ Weights the solution vector by the frequency of the initial states 
    """
    def weightedPassageTime(self, times, bimolecular=False):

        sumTime = 0.0
        sumStart = 0.0

        if len(self.initial_states) == 0:
            raise ValueError("The number of initial states connected to a final state is zero.")

        for state in self.initial_states:
            stateindex = self.stateIndex[state]
            sumTime += self.initial_states[state].count * times[stateindex]
            sumStart += self.initial_states[state].count

        if not bimolecular:
            return sumTime / sumStart
        else:
            return (sumTime / sumStart) * (self.build.options.join_concentration)

    def __str__(self):

        output = "This is a builder-rate object. Printing all states and transitions \n"
        for state in self.statespace:
            output += str(state) + " \n"
        output += "Transitions \n"

        # first, set all transitions
        for state in self.statespace:
            for neighbor in self.neighbors[state]:
                myRate, revRate = self.get_rate(state, neighbor)
                transitionlist = self.build.protoTransitions[(state, neighbor)]
                output += "s1 = " + str(state) + "  s2 = " + str(neighbor) + " for = " + "%.2E" % myRate + " back = " + "%.2E" % revRate + "   tlist = " + str(transitionlist) + "\n"

        output += str(self.rate_matrix_csr.toarray())
        return output

    def numOfTransitions(self):

        count = 0
        # first, count all transitions
        for state in self.statespace:
            for neighbor in self.neighbors[state]:
                transitionlist = self.build.protoTransitions[(state, neighbor)]
                count += len(transitionlist)
        return count

    def statsInfo(self):

        output = str(len(self.statespace)) + "   -    " + str(self.numOfTransitions())
        output += "   -    " + str(len(self.initial_states)) + "   -    " + str(len(self.final_states)) + "   -   (builderRate Object)"
        return output
