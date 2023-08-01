# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
This module has three main classes

MergeResult: Given multistrand results, compute reaction rate constants.
MergeSim:    A class for running multistrand concurrently
Bootstrap:   A convenience class to compute bootstrapped confidence bounds.

MergeResult children: FirstStepRate, FirstStepLeakRate, FirstPassageRate

The MergeSim.run() method always returns of the MergeResult sub-classes.
"""

import time
import datetime
import math
import copy
import sys

import numpy as np
import multiprocess

from .options import Literals
from .system import SimSystem
from .utils import printTrajectory
from .__init__ import __version__

MINIMUM_RATE = 1e-36


class MergeResult:

    """ Endstates are not saved for first step leak mode """

    def __init__(self, dataset=None, endStates=None):
        if dataset == None:
            dataset = []
        if endStates == None:
            endStates = []
 
        # save the dataset for re-sampling and merging results. Also needed for certain rates
        self.myBootstrap = None
         
        self.nForward = 0
        self.nReverse = 0
        self.nForwardAlt = 0
        self.nTotal = 0
         
        self.dataset = dataset
        self.endStates = endStates
         
        self.generateCounts()
            
    def generateCounts(self):
        # Pre-computing some metrics
        self.nForward = sum((i.tag == Literals.success for i in self.dataset))
        self.nReverse = sum((i.tag == Literals.failure for i in self.dataset))
        self.nForwardAlt = sum((i.tag == Literals.alt_success for i in self.dataset))
        self.nTotal = len(self.dataset)
        
    def merge(self, that, deepCopy=False):
        # Now merge the existing datastructures with the ones from the new dataset
        if deepCopy:
            for data in that.dataset:
                self.dataset.append(data)
            for state in that.endStates:
                self.endStates.append(copy.deepcopy(state))                
        else:
            self.dataset.append(that.dataset)
            self.endStates.append(that.endStates)
        
        self.generateCounts()
        
    """ Convenience methods  """

    def doBootstrap(self, NIn=1000):
        self.myBootstrap = Bootstrap(self, N=NIn, computek1=True)
        return self.myBootstrap.ninetyFivePercentiles()

    def log10KEff(self, concentration):
        return np.log10(self.kEff(concentration))

    def testForTwoStateness(self, concentration):
        print("Test for two-stateness is not implemented for this object (type: " + type(self).__name__ + ")\n")

    def __str__(self):
        return (f"{self.__class__.__name__}:\n" + (
            "" if self.myBootstrap is None else f"{self.myBootstrap}\n"))


# Migration rates for first step
class FirstStepRate(MergeResult):

    def sumCollisionForward(self):
        return sum((np.float64(i.collision_rate) for i in self.dataset if i.tag == Literals.success))

    def sumCollisionForwardAlt(self):
        return sum((np.float64(i.collision_rate) for i in self.dataset if i.tag == Literals.alt_success))

    def sumCollisionReverse(self):
        return sum((np.float64(i.collision_rate) for i in self.dataset if i.tag == Literals.failure))

    def weightedForwardUni(self):
        mean_collision_forward = np.float64(self.sumCollisionForward()) / np.float64(self.nForward)
        weightedForwardUni = sum((np.float64(i.collision_rate) * np.float64(i.time) for i in self.dataset if i.tag == Literals.success))

        return weightedForwardUni / (mean_collision_forward * np.float64(self.nForward))

    def weightedReverseUni(self):
        if self.nReverse == 0:
            return np.float64(0)

        mean_collision_reverse = np.float64(self.sumCollisionReverse()) / np.float64(self.nReverse)
        weightedReverseUni = sum((np.float64(i.collision_rate) * np.float64(i.time) for i in self.dataset if i.tag == Literals.failure))

        return weightedReverseUni / (mean_collision_reverse * np.float64(self.nReverse))

    def k1(self):
        if self.nForward == 0:
            return MINIMUM_RATE
        else:
            return self.sumCollisionForward() / np.float64(self.nTotal)

    def k1Alt(self):
        if self.nForwardAlt == 0:
            return MINIMUM_RATE
        else:
            return self.sumCollisionForwardAlt() / np.float64(self.nTotal)

    def k1Prime(self):
        if self.nReverse == 0:
            return MINIMUM_RATE
        else:
            return self.sumCollisionReverse() / np.float64(self.nTotal)

    def k2(self):
        if self.nForward == 0:
            return MINIMUM_RATE
        else:
            return np.float64(1.0) / self.weightedForwardUni()

    def k2Prime(self):
        if self.nReverse == 0:
            return MINIMUM_RATE
        else:
            return np.float64(1.0) / self.weightedReverseUni()

    def kEff(self, concentration=None):
        if concentration == None:
            print("Cannot compute k_effective without concentration")
            return MINIMUM_RATE

        concentration = np.float64(concentration)

        if self.nForward == 0:
            return MINIMUM_RATE

        # expected number of failed collisions
        multiple = (self.k1Prime() / self.k1())

        # the expected rate for a collision
        collTime = self.k1() + self.k1Prime()

        dTForward = np.float64(1.0) / self.k2() + np.float64(1.0) / (concentration * collTime)
        dTReverse = np.float64(1.0) / self.k2Prime() + np.float64(1.0) / (concentration * collTime)
        dT = dTReverse * multiple + dTForward

        return (np.float64(1.0) / dT) * (np.float64(1.0) / concentration)

    def testForTwoStateness(self, concentration=None):
        if concentration == None:
            print("Warning! Attempting to test for two-state behaviour but concentration was not given. \n")
            return True

        # Test if the failed trajectory and the success trajectory are dominated ( > 10% of total ) by the unimolecular phase
        tau_bi_succ = np.float64(1) / (self.k1() * concentration)
        tau_bi_fail = np.float64(1) / (self.k1Prime() * concentration)

        testFail = (tau_bi_fail / self.weightedReverseUni()) < 9
        testSucces = (tau_bi_succ / self.weightedForwardUni()) < 9

        if(testFail | testSucces):
            print(''.join(["Warning! At the chosen concentration, ",
                           str(concentration),
                           " M, the reaction might violate two-state secondary order rate kinetics"]))
            print(self)

    def resample(self):
        # returns a new rates object with resampled data
        N = len(self.dataset)
        newDataset = np.random.choice(self.dataset, N, True).tolist()
        return FirstStepRate(newDataset)

    def castToNumpyArray(self):
        self.collision_rates = np.array(self.collision_rates)
        self.was_success = np.array(self.was_success)
        self.was_failure = np.array(self.was_failure)

    # # override toString
    def __str__(self):
        output = super(FirstStepRate, self).__str__()
        output += f"  nForward = {self.nForward}\n"
        output += f"  nReverse = {self.nReverse}\n"

        if self.nForward > 0:
            output += "  k1  = %.3g /M /s\n" % self.k1()
        if self.nReverse > 0:
            output += "  k1' = %.3g /M /s\n" % self.k1Prime()
        if self.nForward > 0:
            output += "  k2  = %.3g /M /s\n" % self.k2()
        if self.nReverse > 0:
            output += "  k2' = %.3g /M /s\n" % self.k2Prime()

        # suc = (x for x in self.dataset if (x.tag==Literals.success or x.tag==Literals.failure))
        # for x in suc:
        #     output+= x.__str__() +"\n\n"
        return output


class FirstStepLeakRate(MergeResult):
    
    def __init__(self, dataset=None, endStates=None):
        super(FirstStepLeakRate, self).__init__(dataset, endStates)
    
        """ Now discard all non-success data """
        self.dataset = [x for x in self.dataset if ((x.tag == Literals.success) or x.tag == Literals.alt_success)]

    def sumCollisionForward(self):
        return sum((np.float64(i.collision_rate) for i in self.dataset if i.tag == Literals.success))

    def sumCollisionForwardAlt(self):
        return sum((np.float64(i.collision_rate) for i in self.dataset if i.tag == Literals.alt_success))

    def k1(self):
        if self.nForward == 0:
            return MINIMUM_RATE
        else:
            return self.sumCollisionForward() / np.float64(self.nTotal)

    def k1Alt(self):
        if self.nForwardAlt == 0:
            return MINIMUM_RATE
        else:
            return self.sumCollisionForwardAlt() / np.float64(self.nTotal)

    def resample(self):
        new_dataset = []
        time_outs = self.nTotal - self.nForward - self.nReverse - self.nForwardAlt
        successful_trials = len(self.dataset)
        p = np.float64(successful_trials) / self.nTotal
        # the number of succesful trials
        success = np.random.binomial(self.nTotal, p)

        # use random.choice rather than random.sample because this samples
        # WITH REPLACEMENT, as required.
        if success > 0:
            new_dataset = np.random.choice(self.dataset, success, True).tolist()
        else:
            new_dataset = []

        resampledRates = FirstStepLeakRate(new_dataset)
        # we have passed in an only successful dataset here. Hence our metrics will NOT
        # be correct, only those for successful reactions
        resampledRates.nTotal = self.nTotal
        resampledRates.nReverse = self.nTotal - success

        return resampledRates

    def merge(self, that, deepCopy=True):
        # that is always a FirstStepLeakRate object
        if deepCopy:
            for data in that.dataset:
                self.dataset.append(copy.deepcopy(data))
        else:
            self.dataset.append(that.dataset)

        self.nForward += that.nForward
        self.nForwardAlt += that.nForwardAlt
        self.nReverse += that.nReverse
        self.nTotal += that.nTotal

    def __str__(self):
        output = super(FirstStepLeakRate, self).__str__()
        output += f"  nForward = {self.nForward}\n"
        output += f"  nReverse = {self.nReverse}\n"

        if self.nForwardAlt > 0:
            output += "  nForwardAlt = " + str(self.nForwardAlt) + "\n"
        if self.nForward > 0:
            output += "  k1    = %.3g /M /s\n" % self.k1()
        if self.nForwardAlt > 0:
            output += "  k1Alt = %.3g /M /s\n" % self.k1Alt()
        return output


""" Uses data from first passage time rather than first step mode"""
class FirstPassageRate(MergeResult):

    def resample(self):
        N = len(self.dataset)
        newDataset = np.random.choice(self.dataset, N, True).tolist()

        return FirstPassageRate(newDataset)

    def k1(self):
        mean = np.mean([i.time for i in self.dataset])
        return np.float64(1.0) / mean

    def kEff(self, concentration):
        mean = np.mean([i.time for i in self.dataset])
        kEff = np.float64(1.0) / (mean * concentration)
        return kEff

    def __str__(self):
        output = super(FirstPassageRate, self).__str__()
        output += "  nForward: %d\n" % self.nForward
        output += "  k1 = %.3g /M /s\n" % self.k1()
        return output


class Bootstrap():

    def __init__(self, myRates, N=1000, concentration=None, computek1=False, computek1Alt=False):

        self.myRates = myRates

        # # Resample the dataset
        # FD: This is more expensive than strictly required.
        # FD: Note that this computes the CI for kEff().
        b_start_time = time.time()
        self.effectiveRates = list()
        self.effectiveAltRates = list()
        self.logEffectiveRates = list()
        self.N = N

        print("Bootstrapping " + type(myRates).__name__ + ", using " + str(self.N) + " samples.", end=' ')

        for i in range(self.N):
            # create a new sample, with replacement
            sample = self.myRates.resample()

            # compute k1
            if(computek1):
                rate = float(sample.k1())
            else:
                rate = float(sample.kEff(concentration))
                
            self.effectiveRates.append(rate)

            if(computek1Alt):
                # print "we did resample here"
                rate = float(sample.k1Alt())
                self.effectiveAltRates.append(rate)
            del sample

        # sort for percentiles
        self.effectiveRates.sort(cmp=None, key=None, reverse=False)
        self.effectiveAltRates.sort(cmp=None, key=None, reverse=False)
        b_finish_time = time.time()
        print("   ..finished in %.2f sec.\n" % (b_finish_time - b_start_time))

        # Yet to generate log alt rates
        for rate in self.effectiveRates:
            self.logEffectiveRates.append(np.log10(rate))

    def ninetyFivePercentiles(self):
        low = self.effectiveRates[int(0.025 * self.N)]
        high = self.effectiveRates[int(0.975 * self.N)]
        return low, high

    def ninetyFivePercentilesAlt(self):
        low = self.effectiveAltRates[int(0.025 * self.N)]
        high = self.effectiveAltRates[int(0.975 * self.N)]
        return low, high

    def standardDev(self):
        return np.std(self.effectiveRates)

    def logStd(self):
        # Returns standard deviation from the mean of log10 of the original rates
        return np.std(self.logEffectiveRates)

    def __str__(self):
        if len(self.effectiveRates) > 0 : 
            low, high = self.ninetyFivePercentiles()    
            output = "95%% Confidence Interval: [%.3g /M /s, %.3g /M /s] with N=%d" % (low, high, self.N)
            return output
        else :
            return "No successful reactions observed "


# # Concurrent classes start here ==============================================


class optionsFactory:

    def __init__(self, funct, put0, put1, put2, put3, put4, put5, put6):
        self.myFunction = funct
        self.input0 = put0
        self.input1 = put1
        self.input2 = put2
        self.input3 = put3
        self.input4 = put4
        self.input5 = put5
        self.input6 = put6

    def new(self, inputSeed):
        # The input0 is always trials.
        assert isinstance(self.input0, int)

        output = None
        if self.input1 == None:
            output = self.myFunction(self.input0)
        elif self.input2 == None:
            output = self.myFunction(self.input0, self.input1)
        elif self.input3 == None:
            output = self.myFunction(self.input0, self.input1, self.input2)
        elif self.input4 == None:
            output = self.myFunction(
                self.input0, self.input1, self.input2, self.input3)
        elif self.input5 == None:
            output = self.myFunction(
                self.input0, self.input1, self.input2, self.input3, self.input4)
        elif self.input6 == None:
            output = self.myFunction(
                self.input0, self.input1, self.input2, self.input3, self.input4, self.input5)
        else:
            output = self.myFunction(
                self.input0, self.input1, self.input2, self.input3, self.input4, self.input5, self.input6)

        if output == None:
            sys.exit("MergeSim error: Did not recieve Options object from the factory function.")
        output.initial_seed = inputSeed
        return output


class MergeSimSettings:

    RESULTTYPE1 = "FirstStepRate"
    RESULTTYPE2 = "FirstStepLeakRate"
    RESULTTYPE3 = "FirstPassageRate"
    
    DEFAULT_OUTPUT_FILE = "mergesim_output.log"

    def __init__(self):
        self.debug = False
        self.resultsType = self.RESULTTYPE1
        self.terminationCount = None
        self.max_trials = 250000000
        self.timeOut = 24 * 60 * 60
        
        self.saveInterval = 500000
        self.saveIncrement = self.saveInterval
    
        self.outputFile = self.DEFAULT_OUTPUT_FILE
    
        self.bootstrap = False
        self.bootstrapN = 0

    def rateFactory(self, dataset=None, endStates=None):
        if self.resultsType == self.RESULTTYPE1:
            return FirstStepRate(dataset, endStates)
        if self.resultsType == self.RESULTTYPE2:
            return FirstStepLeakRate(dataset)
        if self.resultsType == self.RESULTTYPE3:
            return FirstPassageRate(dataset, endStates)

    def shouldTerminate(self, printFlag, nForwardIn, nReverseIn, timeStart):
        if self.terminationCount == None:
            return True
        else:
            if printFlag:
                print("nForward = %i " % nForwardIn.value)
                print("nReverse = %i \n" % nReverseIn.value)
            if nForwardIn.value >= self.terminationCount:
                print("Found " + str(nForwardIn.value) + " successful trials, terminating.")
                return True
            elif nForwardIn.value + nReverseIn.value > self.max_trials:
                print("Simulated " + str(nForwardIn.value + nReverseIn.value) + " trials, terminating.")
                return True
            elif time.time() - timeStart > self.timeOut:
                print("Runtime exeeded " + str(self.timeOut) + " seconds, terminating.")
                return True
            else:
                return False

    def setOutputFile(self, filetitle):
            # sets flags and updates the file title but not more than this.
            # Here, we will wish to append to the file in the future but not normally.
            self.outputFile = filetitle + time.strftime("-%d-%m-%y+%H-%M-%S") + ".log"

    def setBootstrap(self, doBootstrap, numIn):
        self.bootstrap = doBootstrap
        self.bootstrapN = numIn


def timeStamp(inTime=None):

    if inTime == None:
        inTime = time.time()
    return str(datetime.datetime.fromtimestamp(inTime).strftime('%Y-%m-%d %H:%M:%S'))


class MergeSim:
    """
    This class has two modus operandi:

    1) `self.settings.resultsType == self.settings.RESULTTYPE1`  (default)

          In this mode, `SimSystem.results` and `SimSystem.end_states` are made
          available as `myMR = MergeSim.run()` via `myMR.dataset` and
          `myMR.endStates` (renaming is deliberate)

    2) `self.settings.resultsType != self.settings.RESULTTYPE1`

          In this mode, the system generates `MergeResults` objects that compute
          `k1(), kEff()` (possibly more) in a multi-processing fashion.

    In addition, users can `setAnaylsisFactory()` for custom in-process analysis
    (for example, computing metrics on traces).
    """

    numOfThreads = 2
    seed = 7713147777
    ctx = multiprocess.get_context()

    def __init__(self, settings=None):
        self.initializationTime = time.time()
        if multiprocess.current_process().name == "MainProcess":
            print(f"{timeStamp()}   Starting Multistrand {__version__}  "
                  f"(c) 2008-2023 Caltech")

        self.factory = optionsFactory
        self.aFactory = None
        if settings == None:
            self.settings = MergeSimSettings()

    # The argument is the count of successfull trials before stopping the simulation
    def setTerminationCriteria(self, terminationCount=25):
        self.settings.terminationCount = terminationCount

    def setFirstStepMode(self):
        self.settings.resultsType = self.settings.RESULTTYPE1

    def setLeakMode(self):  # leak mode is a  type of First Step mode
        self.settings.resultsType = self.settings.RESULTTYPE2

    def setPassageMode(self):
        self.settings.resultsType = self.settings.RESULTTYPE3

    def setOutputFile(self, title):
        self.settings.setOutputFile(title)

    def setBootstrap(self, bootstrapIn, num):
        self.settings.setBootstrap(bootstrapIn, num)

    def timeSinceStart(self):
        print("Time since creating object %.5f seconds" %
              (time.time() - self.initializationTime))

    def setDebug(self, value):
        self.settings.debug = value

    def setNumOfThreads(self, numOfThreads):
        self.numOfThreads = numOfThreads

    # this can be re-done using an args[] obj.
    def setOptionsFactory(self, optionsFactory):
        self.factory = optionsFactory

    def setOptionsFactory1(self, myFun, put0):
        self.factory = optionsFactory(
            myFun, put0, None, None, None, None, None, None)

    def setOptionsFactory2(self, myFun, put0, put1):
        self.factory = optionsFactory(
            myFun, put0, put1, None, None, None, None, None)

    def setOptionsFactory3(self, myFun, put0, put1, put2):
        self.factory = optionsFactory(
            myFun, put0, put1, put2, None, None, None, None)

    def setOptionsFactory4(self, myFun, put0, put1, put2, put3):
        self.factory = optionsFactory(
            myFun, put0, put1, put2, put3, None, None, None)

    def setOptionsFactory5(self, myFun, put0, put1, put2, put3, put4):
        self.factory = optionsFactory(
            myFun, put0, put1, put2, put3, put4, None, None)

    def setOptionsFactory6(self, myFun, put0, put1, put2, put3, put4, put5):
        self.factory = optionsFactory(
            myFun, put0, put1, put2, put3, put4, put5, None)

    def setOptionsFactory7(self, myFun, put0, put1, put2, put3, put4, put5, put6):
        self.factory = optionsFactory(
            myFun, put0, put1, put2, put3, put4, put5, put6)

    def setAnaylsisFactory(self, aFactoryIn):
        """
        If the analysis factory is set,
        then perform an threaded analysis of the returned data.
        E.g. the analysis factory receives a set of locks and
        the simulation will call aFactory.doAnalysis(options)
        once for each thread.
        This is used really only for one case study, but could be reused in the future
        """
        lockArray = list()
        for i in range(16):
            lockArray.append(multiprocess.Lock())

        self.aFactory = aFactoryIn
        self.aFactory.lockArray = lockArray

    def clearAnalysisFactory(self):
        """
        Reset the multithreading objects.
        """
        if not self.aFactory == None:
            self.aFactory.clear()

    def printTrajectories(self, myOptions):
        """
        Print all the trajectories we can find.
        Debug function primairly.
        """
        print("Printing trajectory and times: \n")

        trajs = myOptions.full_trajectory
        times = myOptions.full_trajectory_times

        for t, time in zip(trajs, times):
            print(t, "  t=", time, "\n")

    def printTrajectory(self):
        instanceSeed = self.seed + (time.time() * 10000) % (math.pow(2, 32) - 1)
        o1 = self.factory.new(instanceSeed)

        o1.num_simulations = 1
        o1.output_interval = 1

        s = SimSystem(o1)
        s.start()

        printTrajectory(o1)

    def initialInfo(self):
        myOptions = self.factory.new(777)
        simSys = SimSystem(myOptions)
        simSys.initialInfo()

    def startSimMessage(self):
        welcomeMessage = "Using Results Type: " + self.settings.resultsType + "\n"
        welcomeMessage += ''.join(["Computing ", str(
            self.numOfThreads * self.trialsPerThread), " trials, using ", str(self.numOfThreads), " threads .. \n"])

        if not self.settings.terminationCount == None:
            welcomeMessage += " .. and rolling " + str(self.trialsPerThread)
            welcomeMessage += " trajectories per thread until " + str(self.settings.terminationCount) + " successful trials occur. \n"
        return welcomeMessage

    def createOutputMessage(self):
        outputString = "Trials:              " + str(self.factory.input0) + "\n"
        outputString += "Max Trials:          " + str(self.settings.max_trials) + "\n"
        outputString += "Termination  Count:  " + str(self.settings.terminationCount) + "\n"
        outputString += "Num. of Threads:     " + str(self.numOfThreads) + "\n"
        outputString += "Run Time:            " + str(self.runTime) + "s \n"
        outputString += "Using Results Type:  " + self.settings.resultsType + "\n\n"

        def hiddenPrint(myString):
            options = self.factory.new(0)
            myString += "Temperature: " + str(options.temperature) + "K \n\n"
            myString += "Start States\n"
            for x in options.start_state:
                myString += x.__str__() + "\n\n"
            for x in options.stop_conditions:
                myString += x.__str__() + "\n"

        myProc = self.ctx.Process(target=hiddenPrint, args=[outputString])
        myProc.start()
        myProc.join()
        myProc.terminate()

        outputString += "\n" + self.results.__str__() + "\n"
        return outputString
    
    def writeToFile(self):
        outputString = self.createOutputMessage()
        
        if(self.settings.outputFile == self.settings.DEFAULT_OUTPUT_FILE):
            # if a file has been specified, open is default mode.
            f = open(self.settings.outputFile, 'a+')
        else:
            # i.e. the user has not passed in any data at all in terms of file logging.
            # Open in write, not append mode!!
            f = open(self.settings.outputFile, 'w')
            
        f.write(outputString)
        f.close()
          
    def printStates(self):
        def actualPrint():
            print("Start states:")
            for i in self.factory.new(0).start_state:
                print(i)
                print()
                
            print("Stop conditions: ")
            for i in self.factory.new(0).stop_conditions:
                print(i)
                print()
             
        myProc = self.ctx.Process(target=actualPrint, args=[])
        myProc.start()
        myProc.join()
        myProc.terminate()

    def run(self):
        self.trialsPerThread = int(
            math.ceil(float(self.factory.input0) / float(self.numOfThreads)))
        startTime = time.time()

        assert(self.numOfThreads > 0)

        manager = multiprocess.Manager()

        self.exceptionFlag = manager.Value('b', True)
        self.managed_result = manager.list()
        self.managed_endStates = manager.list()
        self.nForward = manager.Value('i', 0)
        self.nReverse = manager.Value('i', 0)

        self.results = self.settings.rateFactory()
        self.endStates = []

        def doSim(myFactory, aFactory, list0, list1, instanceSeed, nForwardIn, nReverseIn):
            try:
                myOptions = myFactory.new(instanceSeed)
                myOptions.num_simulations = self.trialsPerThread
            except:
                self.exceptionFlag.value = False
                return
            
            # Overwrite the result factory method if we are not using First Step Mode.
            # By default, the results object is a First Step object.
            if myOptions.simulation_mode != Literals.first_step:
                self.setPassageMode()

            try:
                s = SimSystem(myOptions)
                s.start()
            except:
                self.exceptionFlag.value = False
                return

            myFSR = self.settings.rateFactory(myOptions.interface.results)
            nForwardIn.value += myFSR.nForward + myFSR.nForwardAlt
            nReverseIn.value += myFSR.nReverse

            for result in myOptions.interface.results:
                list0.append(result)
            for endState in myOptions.interface.end_states:
                list1.append(endState)
            if self.settings.debug:
                self.printTrajectories(myOptions)
            if not(aFactory == None):
                aFactory.doAnalysis(myOptions)

        def getSimulation(input):
            instanceSeed = self.seed + input * 3 * 5 * 19 + (time.time() * 10000) % (math.pow(2, 32) - 1)
            return multiprocess.Process(target=doSim, args=(
                self.factory, self.aFactory, self.managed_result,
                self.managed_endStates, instanceSeed, self.nForward, self.nReverse))

        # this saves the results generated so far as regular Python objects,
        # and clears the concurrent result lists.
        def saveResults():
            if self.settings.terminationCount == None:
                # just let the threads join peacefully
               for i in range(self.numOfThreads):
                   procs[i].join()
            else :
                # join all running threads -- the process has 999 seconds to close or it will be terminated.
                for i in range(self.numOfThreads):
                    procs[i].join(timeout=999)
                    procs[i].terminate()

            self.runTime = (time.time() - startTime)
            print("Done.  %.5f seconds -- now processing results \n" % (time.time() - startTime))

            # Leak - the below is a leak rates object
            # NB: Initialize with a dataset, but we merge with
            # a differrent rates object.
            myFSR = self.settings.rateFactory(self.managed_result, self.managed_endStates)

            self.results.merge(myFSR, deepCopy=True)

            # save the terminal states if we are not in leak mode
            if self.settings.resultsType == self.settings.RESULTTYPE2:
                print(self.results)

            # reset the multiprocessing results lists.
            self.managed_result = manager.list()
            self.managed_endStates = manager.list()
            # this should also reset the 
            self.settings.saveInterval += self.settings.saveIncrement 

        # give a print of the initial states and stopping conditions
        self.printStates()

        # start the initial bulk
        print(self.startSimMessage())
        
        procs = []
        for i in range(self.numOfThreads):
            p = getSimulation(i)
            procs.append(p)
            p.start()
                
        printFlag = False

        # check for stop conditions, restart sims if needed
        while self.exceptionFlag.value:
            if self.settings.shouldTerminate(printFlag, self.nForward, self.nReverse, startTime):
                break
            printFlag = False
            # find and re-start finished threads
            for i in range(self.numOfThreads):
                if not procs[i].is_alive():
                    procs[i].join()
                    procs[i].terminate()
                    procs[i] = getSimulation(i)
                    procs[i].start()
                    printFlag = True
            time.sleep(0.2)

            # if >500 000 results have been generated, then store
            if (self.nForward.value + self.nReverse.value) > self.settings.saveInterval:
                saveResults()

        if not self.exceptionFlag.value:
            raise Exception("MergeSim: exception found in child process.")

        saveResults()

        if not self.settings.resultsType == MergeSimSettings.RESULTTYPE2:
            self.results.generateCounts()
        
        if self.settings.bootstrap == True:
            self.results.doBootstrap(self.settings.bootstrapN)

        self.writeToFile()    
        return 0
