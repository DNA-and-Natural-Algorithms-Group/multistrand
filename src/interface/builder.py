'''
Created on Oct 2, 2017

    Frits Dannenberg, June 1th, 2018.

    The builder code has two classes: Builder and BuilderRate.
    A Builder will accept references to a function that returns a multistrand options object. 
    The sampled states and transitions will be collected. 
    
    The BuilderRate accepts a Builder object and will compute the 
    mean first passage time from the starting states until hitting an end state.    
    
'''

import time, copy, os, sys

from multistrand.system import SimSystem
from multistrand.utils import uniqueStateID, seqComplement
from multistrand.options import Options, Literals
from multistrand.experiment import standardOptions, makeComplex

from scipy.sparse import csr_matrix, coo_matrix, csc_matrix
from scipy.sparse.linalg import spsolve, bicg, bicgstab, cg, cgs, gmres, lgmres, qmr, inv

import numpy as np

"""
    Constructs a surface of  N *  ( N -1 ) / 2 points along the reaction frontier.
    This method returns a list of starting states that should be used as initial states for the string method.
    A starting state is a list of complexes.
"""


def hybridizationString(seq):
    
    cutoff = 0.45    
    
    N = len(seq)
    output = []
    
    """ start with the seperated, zero-bp state """
    dotparen1 = "."*N + "+" + "."*N
    
    complex0 = makeComplex([seq], "."*N)
    complex1 = makeComplex([seqComplement(seq)], "."*N)
    seperated = [complex0, complex1]
    
    i = 65
    for complex in seperated:
        for strand in complex.strand_list:
            strand.id = i
            strand.name = "auto_" + str(i)
            i += 1
        
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
                complex = makeComplex(dotparen0, dotparen1)
            
                i = 65
                for strand in complex.strand_list:
                    strand.id = i
                    strand.name = "auto_" + str(i)
                    i += 1
                        
                output.append([complex])
                
    """ We always require the final state to be the success state"""
    complex0 = makeComplex(dotparen0, "("*N + "+" + ")"*N)
         
    i = 65
    for strand in complex0.strand_list:
        strand.id = i
        strand.name = "auto_" + str(i)
        i += 1
         
    output.append([complex0])
#     
                
    return output


class ConvergeCrit(object):
    
    iit = 0
    currRate = -1.0
    
    period = 4  # average out over past X increases.
    array = [-99.0] * period
    precision = 0.05
    
    def converged(self, rateIn):  # also saves the rateIn
        
        averageL = (sum(self.array) / self.period) * (1.0 - self.precision)
        averageH = (sum(self.array) / self.period) * (1.0 + self.precision)
        
        conv1 = rateIn > averageL
        conv2 = rateIn < averageH 
        
        self.array[self.iit % self.period] = rateIn
        self.iit += 1;
        
        return (conv1 and conv2)
    
    def __str__(self):
        
        return  str(self.array)


class transitiontype(object):
    
    unimolecular = "uni"
    bimolecularIn = "bi-in"
    bimolecularOut = "bi-out"
    
    array = [unimolecular, bimolecularIn, bimolecularOut ]
    

class localtype(object):
    
    end = "End"
    loop = "Loop"
    stack = "Stack"
    stackstack = "StackStack"
    loopend = "LoopEnd"
    stackend = "StackEnd"
    stackloop = "StackLoop"
    
    array = [end, loop, stack, stackstack, loopend, stackend, stackloop]
    

class energy(object):
    
    GAS_CONSTANT = 0.0019872036  # kcal / K mol
        
    dH = 0.0;
    dS = 0.0;

    def __init__(self, dG, dH, temp, concentration, n_complexes, n_strands):
        
        assert(200 < temp < 400)
        
        RT = self.GAS_CONSTANT * temp;
        dG_volume = RT * (n_strands - n_complexes) * np.log(1.0 / concentration)
        
        self.dH = dH
        self.dS = -(dG - dG_volume - dH) / temp

    def dG(self, temp):
        
        return self.dH - temp * self.dS
        
    def __str__(self):
        
        return "dH= " + str(self.dH) + " dS= " + str(self.dS)

    def __eq__(self, that):
        if isinstance(that, energy):
            return self.dH == that.dH and self.dS == that.dS
        else:
            raise Exception("Not an acceptible argument", "__eq__", that)


def codeToDesc(code):
    
    output = []
    primes = [5, 7, 11, 13, 17, 19]
    primesq = [25, 49, 121, 169, 289, 361]

    for i in range(len(primes)):

        if code % primes[i] == 0:
            output.append(localtype.array[i])
            
            if code % primesq[i] == 0:  # code is a square number
                output.append(localtype.array[i])
                return output

    return output


class initCountFlux(object):

    def __init__(self):
        self.count = 0;
        self.flux = 0.0;

    def __add__(self, x):
        self.count + x.count

    def __str__(self):
        return "Count: " + str(self.count) + "  flux: " + str(self.flux)

    def  __repr__(self):
        return str(self)


class Builder(object):
    
    verbosity = False
    
    # input function returns the multistrand options object for which to build the statespace
    def __init__(self, inputFunction, arguments):
    
        # initial argument has to be the number of trials
        self.optionsFunction = inputFunction
        self.optionsArgs = copy.deepcopy(arguments)
    
        self.doMultiprocessing = False
        self.printTimer = True
        self.numOfThreads = 8
    
        self.protoSpace = dict()  # key: states. Value: energy
        self.protoTransitions = dict()  # key: transitions. Value: ArrheniusType (negative if it is a bimolecular transition)
        self.protoInitialStates = dict()  # key: states. Value: a initCountFlux object that tells how many times the state has been the initial state and the join flux (rate)
        self.protoFinalStates = dict()  # key: states: Value: the result of this final state can be SUCCES or FAILURE
        
        self.firstStepMode = True
        self.startTime = time.time()
        
        # save a copy for later processing -- note this copy will not have results attached to it,
        # so it won't have a large memory footprint.
        self.options = self.optionsFunction(self.optionsArgs)
        
        self.the_dir = "p_statespace/"
    
    def __str__(self):
        
        output = "states / transitions / initS / finalS      \n "
        output += str(len(self.protoSpace)) + "   -    " + str(len(self.protoTransitions)) 
        output += "   -    " + str(len(self.protoInitialStates)) + "   -    " + str(len(self.protoFinalStates))
    
        return output
    
    def reset(self):
        
        self.protoSpace.clear()
        self.protoTransitions.clear()
        self.protoInitialStates.clear()
        self.protoFinalStates.clear()
    
    def printOverlap(self, other):
        
        N = len(self.protoSpace)
        
        overlap = 0;
        
        for state in self.protoSpace:
            if state in other.protoSpace:
                overlap += 1
        
        print "The overlap is " + str(100.0 * overlap / N) + " percent. "

    def mergeSet(self, this, that):
        
            for key, val in that.iteritems():
                
                if not key in this:
                    this[key] = val;

    def mergeBuilder(self, other):
    
        self.mergeSet(self.protoSpace, other.protoSpace)
        self.mergeSet(self.protoTransitions, other.protoTransitions)
        self.mergeSet(self.protoInitialStates, other.protoInitialStates)
        self.mergeSet(self.protoFinalStates, other.protoFinalStates)
                
    def parseState(self, line, simulatedTemperature, simulatedConc):
        
        mywords = line.split()
        
        n_complexes = int(mywords[0])
        n_strands = 0
        
        ids = []
        structs = []
        
        for i in range(n_complexes):
            
            ids.append(mywords[1 + i])
            structs.append(mywords[1 + 2 * n_complexes + i])
            n_strands += len(mywords[1 + 2 * n_complexes + i].split('+'))
            
            uniqueID2 = uniqueStateID(ids, structs)
            uniqueID = ""
                    
            for iddd in uniqueID2 :
                uniqueID += str(iddd)
        
        dG = float(mywords[1 + 3 * n_complexes])
        dH = float(mywords[1 + 3 * n_complexes + 1])
        
        energyvals = energy(dG, dH, simulatedTemperature, simulatedConc, n_complexes, n_strands)
    
        return uniqueID, energyvals
         
    """ Runs genAndSavePathsFile until convergence is reached"""

    def genUntilConvergence(self, precision):
        
        crit = ConvergeCrit()
        crit.precision = precision
        
        currTime = -1.0
        
        while not crit.converged(currTime) :
            
            self.genAndSavePathsFile()
            
            if self.verbosity:
                print "Size     = %i " % len(self.protoSpace)

            builderRate = BuilderRate(self)
            currTime = builderRate.averageTimeFromInitial()

        if self.verbosity:
            print "Size     = %i " % len(self.protoSpace)
    
    """ Runs genAndSavePathsFile until convergence is reached,
        given a list of initial states"""

    def genUntilConvergenceWithInitialState(self, precision, initialStates, printMeanTime=False):
         
        crit = ConvergeCrit()
        crit.precision = precision
         
        currTime = -1.0
         
        while not crit.converged(currTime) :

            otherBuilder = Builder(self.optionsFunction, self.optionsArgs)

            startTime = time.time()

            """ Only the first state will count towards the set of initial states """
            ignoreInitial = False
            for state in initialStates:
               otherBuilder.genAndSavePathsFile(supplyInitialState=state, ignoreInitialState=ignoreInitial)
               ignoreInitial = True   
            
            self.mergeBuilder(otherBuilder)
             
            if self.verbosity or printMeanTime:
                print "Size     = %i    ---  bytesize = %f mb" % (len(self.protoSpace), sys.getsizeof (self.protoSpace) / (1024 * 1024))
                print "Size T   = %i    ---  bytesize = %f mb" % (len(self.protoTransitions), sys.getsizeof (self.protoTransitions) / (1024 * 1024))
                print "Time = %f" % (time.time() - startTime) 
                
            del otherBuilder
             
            builderRate = BuilderRate(self)
            currTime = builderRate.averageTimeFromInitial(printMeanTime)
            
            if printMeanTime:
                print "Mean first passage time = %.2E" % currTime
 
        if self.verbosity:
            print "Size     = %i " % len(self.protoSpace)
    
    
    
    def deltaPruning(self, delta):
        
        0
        
    
    
    """
    supplyInitialState : this needs to be a list of complexes if used.
    """

    def genAndSavePathsFile(self, ignoreInitialState=False, supplyInitialState=None):
    
        self.startTime = time.time()
    
        space = dict() 
        transitions = dict()
        initStates = dict()
        finalStates = dict()
        
        # the first argument is always the number of paths
        inputArgs = copy.deepcopy(self.optionsArgs)
        
        def runPaths(optionsF, optionsArgs, space, transitions, initStates, finalStates):
    
            myOptions = optionsF(optionsArgs)
            myOptions.activestatespace = True
        
            if not supplyInitialState == None:
                myOptions.start_state = supplyInitialState
                
            """ Set longer searching time for the initial state. """
            if not ignoreInitialState:
                myOptions.simulation_time = myOptions.simulation_time * 10.0 
        
            simTime = time.time()
        
            s = SimSystem(myOptions)
            s.start()  # after this line, the computation is finished.
            
            if self.verbosity:
                 print "Multistrand simulation is now done,      time = %.2f" % (time.time() - simTime)
            
#             loadTime = time.time()
            
            """ load the space """
            myFile = open(self.the_dir + str(myOptions.interface.current_seed) + "/protospace.txt", "r")
            
            for line in myFile:
                
                uniqueID, energyvals = self.parseState(line, myOptions._temperature_kelvin, myOptions.join_concentration)
    
                if not uniqueID in space:
                    
                    space[uniqueID] = energyvals
                
                elif not space[uniqueID] == energyvals:
                    
                    print "My hashmap contains " + str(uniqueID) + " with energy " + str(space[uniqueID]) + " but found: " + str(energyvals) 
                    print "Line = " + line
                    
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
                
                uID1, ev1 = self.parseState(line2, myOptions._temperature_kelvin, myOptions.join_concentration)
                uID2, ev2 = self.parseState(line3, myOptions._temperature_kelvin, myOptions.join_concentration)
                
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
                print "No initial states found!"
            
            while go_on:
                
                line1 = myLines[index];
                line2 = myLines[index + 1];
    
                index = index + 2  # note the whitespace
                go_on = len(myLines) > index
                    
                uID1, ev1 = self.parseState(line2, myOptions._temperature_kelvin, myOptions.join_concentration)
                count = int(line1.split()[0])
                
                if not uID1 in initStates:
                        
                    newEntry = initCountFlux()
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
    
                uID1, ev1 = self.parseState(line1, myOptions._temperature_kelvin, myOptions.join_concentration)
                tag = line2.split()[0]
    
                if not uID1 in finalStates:
                    finalStates[uID1] = tag
    
            """ Now delete the files as they can get quite large """
            os.remove(self.the_dir + str(myOptions.interface.current_seed) + "/protospace.txt")
            os.remove(self.the_dir + str(myOptions.interface.current_seed) + "/prototransitions.txt")
            os.remove(self.the_dir + str(myOptions.interface.current_seed) + "/protoinitialstates.txt")
            os.remove(self.the_dir + str(myOptions.interface.current_seed) + "/protofinalstates.txt")
    
    #             print "Loading is now done,                     time = %.2f" % (time.time() - loadTime)
    
        runPaths(self.optionsFunction, inputArgs, space, transitions, initStates, finalStates)
    
    #         mergeTime = time.time()
    
        # do not forget to merge the objects back
        self.mergeSet(self.protoSpace, space)
        self.mergeSet(self.protoTransitions, transitions) 
        self.mergeSet(self.protoFinalStates, finalStates)
        
    #         print "Done statespace merging!                 time = %.2f " % (time.time() - mergeTime)
    
        if not ignoreInitialState:    
            
            for key, val in initStates.iteritems():
                if not key in self.protoInitialStates:
                    self.protoInitialStates[key] = val;
                else:
                    self.protoInitialStates[key].count += val.count


# FD: This class is in progress.
# this class takes a builder object and computes the average time between
# starting in an initial state and reaching a final state. 
class BuilderRate(object):
    
    debugPrinter = False;
    
    solveToggle = 2

    # input function returns the multistrand options object for which to build the statespace
    def __init__(self, builderIn):
        
        self.build = builderIn

        if len(self.build.protoFinalStates) == 0 :
            raise ValueError('No final states found.')
        
        if self.debugPrinter:
            print str(self.build)
        
#         stdOptions = standardOptions()
#         stdOptions.DNA23Metropolis()
        
#         self.unimolecular_scaling = stdOptions.unimolecular_scaling
#         self.bimolecular_scaling = stdOptions.bimolecular_scaling
        
#         print "Starting to process states"
#         processTime = time.time()
        self.processStates()  # prunes statespace and creates objects that can be used to create the rate matrx
#         print "Prunning states is done,  %.2f" % (time.time() - processTime) 
#         print "Starting to build matrix"
#         matrixTime = time.time()
        self.setMatrix()  # generates the matrix for the current temperature
#         print "Building matrix is now done, %.2f" % (time.time() - matrixTime)
    
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
#         print "done generating list of neighbors             time = %.2f" % (time.time() - startT) 

        startT = time.time()
        # now prune the statespace to only include states that can reach the final state
        self.statespace = self.findConnectedStates(self.final_states, self.neighbors)
#         print "done pruning the statespace via connected     time = %.2f" % (time.time() - startT) 
        
        startT = time.time()
        # Now re-generated the neighbors, but only for states in the statespace -- and only "forward" transitions instead of both ways.
        self.neighbors = self.genNeighbors(self.build.protoTransitions, self.build.protoSpace)
#         print "done re-generating the neighbors               time = %.2f" % (time.time() - startT) 

        # Update the initial states to only include those connected states
        for state in self.build.protoInitialStates:
            if state in self.statespace:
                self.initial_states[state] = self.build.protoInitialStates[state]

    """Use the builder options object to determine which rates to compute """

    def get_rate(self, state1, state2):
        
        transitionlist = self.build.protoTransitions[(state1, state2)]
        
        if self.build.options.rate_method == Literals.arrhenius:
            return self.arrhenius_rate(state1, state2, transitionlist)
        else :
            return self.metropolis_rate(state1, state2, transitionlist)

    """Uses the Arrhenius kinetic model to calculate the transition rate from state1 to state and the transition rate from state2 to state1. 
    Returns the transition rate from state1 to state2 
    and then also the reverse rate, from state2 to state1 """

    def arrhenius_rate(self, state1, state2, transitionlist):
        
            # # TODO -- fix this
        lnA_left, E_left = codeToDesc[transitionlist[1]]
        lnA_right, E_right = codeToDesc[transitionlist[2]]
        
        lnA = lnA_left + lnA_right
        E = E_left + E_right
        
        DeltaG = self.build.protoSpace[state1].dG(self.build.options._temperature_kelvin)
        DeltaG -= self.build.protoSpace[state2].dG(self.build.options._temperature_kelvin)
                                                
        RT = energy.GAS_CONSTANT * (self.build.options._temperature_kelvin)
        
        uphill = np.e ** (lnA - (DeltaG + E) / RT)
        downhill = np.e ** (lnA - E / RT)
        
        if DeltaG > 0 :
            return uphill, downhill 
        else:
            return downhill, uphill 
        
    """    Returns the transition rate from state1 to state2 
    and then also the reverse rate, from state2 to state1 """

    def metropolis_rate(self, state1, state2, transitionlist):

        myT = self.build.options._temperature_kelvin

        dG1 = self.build.protoSpace[state1].dG(myT)
        dG2 = self.build.protoSpace[state2].dG(myT)
        
        RT = energy.GAS_CONSTANT * myT
        
        if transitionlist[0] == transitiontype.unimolecular:
            
            if dG1 > dG2 :  # state2 is more stable (negative), dG1 - dG2 is positive
                return self.build.options.unimolecular_scaling  , self.build.options.unimolecular_scaling * np.e ** (-(dG1 - dG2) / RT)
            else:  # state2 is less or equally stable (negative), dG1 - dG2 is negative
                return self.build.options.unimolecular_scaling * np.e ** ((dG1 - dG2) / RT), self.build.options.unimolecular_scaling   
        
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
                
#                 iArray.append(stateIndex[neighbor])         
#                 jArray.append(stateIndex[state])

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

#         print "Setting all transitions, temperature is now " + str(self.build.options._temperature_kelvin)

        # first, set all transitions
        for state in self.statespace:
                
            for neighbor in self.neighbors[state]:
                
                myRate, revRate = self.get_rate(state, neighbor)

                if myRate < 0.0 or revRate < 0.0:
                    print "FOUND negative rate!!"

                # This handles either state being a final state (in which case, subtract from the non-final state diagonal,
                # but do not add the transition rate.
                self.addTransition(state, neighbor, myRate, rates, iArray, jArray, self.stateIndex)
                self.addTransition(neighbor, state, revRate, rates, iArray, jArray, self.stateIndex)

        # now actually create the matrix
        rate_matrix_coo = coo_matrix((rates, (iArray, jArray)), shape=(N, N) , dtype=np.float64)
        
        self.rate_matrix_csr = csr_matrix(rate_matrix_coo)
        self.b = -1 * np.ones(N)
        
#         # FD: pre-compute the matrix diagonal for preconditioning
        diagons = [ (1.0 / x, i, j) for x, i, j in zip(rates, iArray, jArray) if i == j ]
        diagons0 = [x[0] for x in diagons]
        diagons1 = [x[1] for x in diagons]
        diagons2 = [x[2] for x in diagons]
        
        diagonal_matrix_coo = coo_matrix((diagons0, (diagons1, diagons2)), shape=(N, N), dtype=np.float64)
        self.rate_matrix_inverse = csr_matrix(diagonal_matrix_coo)
        
        # save two counts for later interest
        self.n_states = N
        self.n_transitions = len(rates)

    def averageTime(self):
        
        # now solve for the first passagetimes

        if self.solveToggle == 1:
            firstpassagetimes, info = bicg(self.rate_matrix_csr, self.b)
            return firstpassagetimes
        
        if self.solveToggle == 2:
            firstpassagetimes, info = bicg(self.rate_matrix_csr, self.b, M=self.rate_matrix_inverse)
            return firstpassagetimes
        
        if self.solveToggle == 3:
            firstpassagetimes, info = bicgstab(self.rate_matrix_csr, self.b)
            return firstpassagetimes

        if self.solveToggle == 4:
            firstpassagetimes, info = bicgstab(self.rate_matrix_csr, self.b, M=self.rate_matrix_inverse)
            return firstpassagetimes

        if self.solveToggle == 5:
            firstpassagetimes, info = gmres(self.rate_matrix_csr, self.b)
            return firstpassagetimes
        
        if self.solveToggle == 6:
            firstpassagetimes, info = gmres(self.rate_matrix_csr, self.b, M=self.rate_matrix_inverse)
            return firstpassagetimes
        
        if self.solveToggle == 7:
            firstpassagetimes, info = cg(self.rate_matrix_csr, self.b)
            return firstpassagetimes
        
        if self.solveToggle == 8:
            firstpassagetimes, info = cg(self.rate_matrix_csr, self.b, M=self.rate_matrix_inverse)
            return firstpassagetimes
        
        if self.solveToggle == 9:
            firstpassagetimes, info = lgmres(self.rate_matrix_csr, self.b)
            return firstpassagetimes
        
        if self.solveToggle == 10:
            firstpassagetimes, info = lgmres(self.rate_matrix_csr, self.b, M=self.rate_matrix_inverse)
            return firstpassagetimes

        firstpassagetimes = spsolve(self.rate_matrix_csr, self.b)
        return firstpassagetimes
        
    def averageTimeFromInitial(self, bimolecular=False, cutCrit=0.01, printMeanTime=True):

        startTime = time.time()

        if self.build.verbosity:
            print "Starting to solve the matrix equation.."
            
        if printMeanTime:
            print "The size of the matrix is %f mb " % ( sys.getsizeof(self.rate_matrix_csr) / (1024  * 1024))
            
            
            
        firstpassagetimes = self.averageTime()  # spsolve(self.rate_matrix_csr  , self.b)        

        sumTime = 0.0
        sumStart = 0.0
        
        for state in self.initial_states:
            
            stateindex = self.stateIndex[state]
            sumTime += self.initial_states[state].count * firstpassagetimes[stateindex]
            sumStart += self.initial_states[state].count
            
#             print "Initial: " + str(state) +  "Count = " + str(self.initial_states[state].count)
        
        self.matrixTime = time.time() - startTime
        
        if self.build.verbosity or printMeanTime:
            print "Solving matrix took %.2f s" % self.matrixTime
        
        cuttOff = cutCrit * (sumTime / sumStart)
        
        self.reducedSize = len([ x < cuttOff for x in firstpassagetimes])
        
        if not bimolecular:
                    return sumTime / sumStart
        else: 
        
            if self.build.verbosity:
                print "Computing concentration-adjusted concentration: " + str(self.build.options.join_concentration)
        
            return (sumTime / sumStart) * (self.build.options.join_concentration)
        
#     def testMatrix(self):
        
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
