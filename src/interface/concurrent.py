"""
Frits Dannenberg, May 17th, 2017.
This module has two parts.
1) Given multistrand results, compute reaction rate constants.
2) A module for running multistrand concurrently. Results are merged automatically
"""
import operator
import os
import time
import datetime
import math
import traceback
import copy
import sys
import random

from collections import Counter
from multistrand.options import Options
import multiprocessing
import numpy as np

from multistrand.system import SimSystem


MINIMUM_RATE = 1e-36

# # Rate-computation classes start here


# Shared class methods
class MergeResult(object):

    myBootstrap = None
    nForward = 0
    nReverse = 0
    nForwardAlt = 0
    nTotal = 0

    dataset = []

    def log10KEff(self, concentration):
        return np.log10(self.kEff(concentration))

    def testForTwoStateness(self, concentration):
        print "Test for two-stateness is not implemented for this object (type: " + type(self).__name__ + ")\n"

    # Convenience.
    def doBootstrap(self, NIn=10000):
        self.myBootstrap = Bootstrap(self, N=NIn,  computek1=True)
        return self.myBootstrap.ninetyFivePercentiles()

    def __str__(self):
        if self.myBootstrap != None:
            return self.myBootstrap.__str__() + "\n"
        else:
            return ""

# # Migration rates for first step
class FirstStepRate(MergeResult):

    def __init__(self, dataset=None):

        if dataset == None:
            dataset = []

        # save the dataset for re-sampling and merging results. Also needed for certain rates
        self.dataset = dataset
        self.generateRates()

    def generateRates(self):

        # Pre-computing some metrics
        self.nForward = sum(
            [i.tag == Options.STR_SUCCESS for i in self.dataset])
        self.nReverse = sum(
            [i.tag == Options.STR_FAILURE for i in self.dataset])
        self.nForwardAlt = sum(
            [i.tag == Options.STR_ALT_SUCCESS for i in self.dataset])

        self.nTotal = len(self.dataset)

    def sumCollisionForward(self):
        return sum([np.float(i.collision_rate) for i in self.dataset if i.tag == Options.STR_SUCCESS])

    def sumCollisionForwardAlt(self):
        return sum([np.float(i.collision_rate) for i in self.dataset if i.tag == Options.STR_ALT_SUCCESS])

    def sumCollisionReverse(self):
        return sum([np.float(i.collision_rate) for i in self.dataset if i.tag == Options.STR_FAILURE])

    def weightedForwardUni(self):

        mean_collision_forward = np.float(
            self.sumCollisionForward()) / np.float(self.nForward)
        weightedForwardUni = sum([np.float(i.collision_rate) * np.float(i.time)
                                  for i in self.dataset if i.tag == Options.STR_SUCCESS])

        return weightedForwardUni / (mean_collision_forward * np.float(self.nForward))

    def weightedReverseUni(self):

        if self.nReverse == 0:
            return np.float(0)

        mean_collision_reverse = np.float(
            self.sumCollisionReverse()) / np.float(self.nReverse)
        weightedReverseUni = sum([np.float(i.collision_rate) * np.float(i.time)
                                  for i in self.dataset if i.tag == Options.STR_FAILURE])

        return weightedReverseUni / (mean_collision_reverse * np.float(self.nReverse))

    def k1(self):

        if self.nForward == 0:
            return MINIMUM_RATE
        else:
            return self.sumCollisionForward() / np.float(self.nTotal)

    def k1Alt(self):
        if self.nForwardAlt == 0:
            return MINIMUM_RATE
        else:
            return self.sumCollisionForwardAlt() / np.float(self.nTotal)

    def k1Prime(self):

        if self.nReverse == 0:
            return MINIMUM_RATE
        else:
            return self.sumCollisionReverse() / np.float(self.nTotal)

    def k2(self):

        if self.nForward == 0:
            return MINIMUM_RATE
        else:
            return np.float(1.0) / self.weightedForwardUni()

    def k2Prime(self):

        if self.nReverse == 0:
            return MINIMUM_RATE
        else:
            return np.float(1.0) / self.weightedReverseUni()

    def kEff(self, concentration=None):

        if concentration == None:
            print "Cannot compute k_effective without concentration"
            return MINIMUM_RATE

        concentration = np.float(concentration)

        if self.nForward == 0:
            return MINIMUM_RATE

        # expected number of failed collisions
        multiple = (self.k1Prime() / self.k1())

        # the expected rate for a collision
        collTime = self.k1() + self.k1Prime()

        dTForward = np.float(1.0) / self.k2() + \
            np.float(1.0) / (concentration * collTime)
        dTReverse = np.float(1.0) / self.k2Prime() + \
            np.float(1.0) / (concentration * collTime)

        dT = dTReverse * multiple + dTForward

        return (np.float(1.0) / dT) * (np.float(1.0) / concentration)

    def testForTwoStateness(self, concentration=None):

        if concentration == None:
            print "Warning! Attempting to test for two-state behaviour but concentration was not given. \n"
            return True

        # Test if the failed trajectory and the success trajectory are dominated ( > 10% of total ) by the unimolecular phase
        tau_bi_succ = np.float(1) / (self.k1() * concentration)
        tau_bi_fail = np.float(1) / (self.k1Prime() * concentration)

        testFail = (tau_bi_fail / self.weightedReverseUni()) < 9
        testSucces = (tau_bi_succ / self.weightedForwardUni()) < 9

        if(testFail | testSucces):

            print(''.join(["Warning! At the chosen concentration, ", str(
                concentration), " M, the reaction might violate two-state secondary order rate kinetics"]))
            print(self)

    def resample(self):
        # returns a new rates object with resampled data

        newDataset = []
        N = len(self.dataset)

        for i in range(N):
            # generate random between 0 and N-1
            index = int(np.floor(np.random.uniform(high=N)))
            newDataset.append(self.dataset[index])

        return FirstStepRate(newDataset)

    def castToNumpyArray(self):

        self.collision_rates = np.array(self.collision_rates)
        self.was_success = np.array(self.was_success)
        self.was_failure = np.array(self.was_failure)

    def merge(self, that, deepCopy=False):
        
        # Now merge the existing datastructures with the ones from the new dataset
        if deepCopy:
            for result in that.dataset:
                self.dataset.append(copy.deepcopy(result))
        else:
            self.dataset.append(that.dataset)

    # # override toString
    def __str__(self):

        output = "nForward = " + str(self.nForward) + " \n"
        output += "nReverse = " + str(self.nReverse) + " \n \n"

        if(self.nForward > 0):
            output += "k1       = %.2e  /M /s  \n" % self.k1()

        if(self.nReverse > 0):
            output += "k1'      = %.2e /M /s \n" % self.k1Prime()

        if(self.nForward > 0):
            output += "k2       = %.2e  /M /s  \n" % self.k2()

        if(self.nReverse > 0):
            output += "k2'      = %.2e /M /s \n" % self.k2Prime()

        output+= super(FirstStepRate, self).__str__()

        # suc = (x for x in self.dataset if (x.tag==Options.STR_SUCCESS or x.tag==Options.STR_FAILURE))
        # for x in suc:
        #     output+= x.__str__() +"\n\n"
        return output

    def shortString(self):

        output = "nForward = " + str(self.nForward) + " \n"
        output += "nReverse = " + str(self.nReverse) + " \n \n"

        if self.nForwardAlt > 0:
            output += "nForwardAlt = " + str(self.nForwardAlt)

        if(self.nForward > 0):
            output += "k1       = %.2e  /M /s  \n" % self.k1()

        return output

    def typeName(self):
        return "First Step Rate"


class FirstStepLeakRate(MergeResult):

    # take a dataset with failed trajectories, save only the important information.
    def __init__(self, dataset=None, generate_rates=True):

        if dataset == None:
            dataset = []

        c = Counter(x.tag for x in dataset)

        self.nForward = c[Options.STR_SUCCESS]
        self.nForwardAlt = c[Options.STR_ALT_SUCCESS]
        self.nReverse = c[Options.STR_FAILURE]
        self.nTotal = len(dataset)

        self.dataset = [x for x in dataset if (
            (x.tag == Options.STR_SUCCESS) or x.tag == Options.STR_ALT_SUCCESS)]

    def generateRates(self, dataset=None):
        pass

    def sumCollisionForward(self):
        return sum([np.float(i.collision_rate) for i in self.dataset if i.tag == Options.STR_SUCCESS])

    def sumCollisionForwardAlt(self):
        return sum([np.float(i.collision_rate) for i in self.dataset if i.tag == Options.STR_ALT_SUCCESS])

    def k1(self):
        if self.nForward == 0:
            return MINIMUM_RATE
        else:
            return self.sumCollisionForward() / np.float(self.nTotal)

    def k1Alt(self):
        if self.nForwardAlt == 0:
            return MINIMUM_RATE
        else:
            return self.sumCollisionForwardAlt() / np.float(self.nTotal)

    def resample(self):

        new_dataset = []
        time_outs = self.nTotal - self.nForward - self.nReverse - self.nForwardAlt
        successful_trials = len(self.dataset)
        p = np.float(successful_trials) / self.nTotal
        # the number of succesful trials
        success = np.random.binomial(self.nTotal, p)

        # use random.choice rather than random.sample because this samples
        # WITH REPLACEMENT, as required.
        if success > 0:
            new_dataset = np.random.choice(
                self.dataset, success, True).tolist()
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

        output = "nForward = " + str(self.nForward) + " \n"
        output += "nReverse = " + str(self.nReverse) + " \n \n"

        # for x in self.dataset:
        #     output+= x.__str__() +"\n\n"
        
        if self.nForwardAlt > 0:
            output += "nForwardAlt = " + str(self.nForwardAlt) + "\n"

        if(self.nForward > 0):
            output += "k1       = %.2e  /M /s  \n\n" % self.k1()


        if(self.nForwardAlt > 0):
            output += "k1Alt       = %.2e  /M /s  \n" % self.k1Alt()
            
        
            

        output+= super(FirstStepLeakRate, self).__str__() + "\n"
        
        
        return output

    def shortString(self):
        output = "nForward = " + str(self.nForward) + " \n"
        output += "nReverse = " + str(self.nReverse) + " \n \n"

        if self.nForwardAlt > 0:
            output += "nForwardAlt = " + str(self.nForwardAlt) + "\n"

        if(self.nForward > 0):
            output += "k1       = %.2e  /M /s  \n\n" % self.k1()

        if(self.nForwardAlt > 0):
            output += "k1Alt       = %.2e  /M /s  \n" % self.k1Alt()

       
        return output

# Like migrationrate, but uses data from first passage time rather than first step mode
class FirstPassageRate(MergeResult):

    def __init__(self, dataset=None):

        if dataset == None:
            # set up for resampling
            self.dataset = []
            self.times = []
            self.timeouts = []

        else:
            self.process(dataset)
            self.generateRates()

    def process(self, dataset):

        self.times = np.array([i.time for i in dataset])
        self.timeouts = [
            i for i in dataset if not i.tag == Options.STR_SUCCESS]

    def castToNumpyArray(self):

        self.times = np.array(self.times)

    def resample(self):

        # avoid resampling time-outs.
        newRates = FirstPassageRate()

        N = len(self.times)

        # generate random between 0 and N-1
        for i in range(N):

            index = int(np.floor(np.random.uniform(high=N)))
            newRates.times.append(self.times[index])

        newRates.castToNumpyArray()
        newRates.generateRates()

        return newRates

    def merge(self, that, deepCopy=False):

        # Now merge the existing datastructures with the ones from the new dataset
        if deepCopy:
            for time in that.times:
                self.times.append(copy.deepcopy(time))
            for timeout in that.timeouts:
                self.timeouts.append(copy.deepcopy(timeout))
        else:
            self.times.append(that.times)
            self.timeouts.append(that.timesouts)
        
        self.generateRates()
            
    def generateRates(self):
            
        self.nForward = len(self.times) - len(self.timeouts) 

    def k1(self):
        
        mean = np.mean(self.times)
        return np.float(1.0) / (mean)
        

    def kEff(self, concentration):

        if len(self.timeouts) > 0:

            print("# association trajectories did not finish =",
                  str(len(self.timeouts)))

            for timeout in self.timeouts:
                print timeout

            for i in self.timeouts:
                assert (i.type_name == 'Time')
                assert (i.tag == None)
                assert (i.time >= 10.0)

        mean = np.mean(self.times)
        kEff = np.float(1.0) / (mean * concentration)

        return kEff

    def __str__(self):
        output = "nForward: %d \n" % self.nForward
        output+= "k1 = %.3g \n" % self.k1()
        output+= super(FirstPassageRate, self).__str__() 
        return output


# Like FirstStepRate, but the collected data is from trajectory mode. 
# Currently, this only saves the final state.
class TrajectoryMerge(MergeResult):

    def __init__(self, dataset=None, endStates = None):

        if dataset == None:
            # set up for resampling
            self.dataset = []
            self.endStates = []
        else :
            self.dataset = dataset
            self.endStates = endStates

    def resample(self):

        # avoid resampling time-outs.
        newRates = TrajectoryMerge()

        return newRates

    def merge(self, that, deepCopy=False):
        
        # Now merge the existing datastructures with the ones from the new dataset
        if deepCopy:
            for result in that.dataset:
                self.dataset.append(copy.deepcopy(result))
            for state in that.endStates:
                self.endStates.append(copy.deepcopy(result))
                
        else:
            self.dataset.append(that.dataset)
            self.endStates.append(that.endStates)
            
    def generateRates(self):
            
        self.nForward = 0;

    def k1(self):
        
        return 0;
        

    def kEff(self, concentration):

        return 0;

    def __str__(self):

        nSucc = len([i for i in self.dataset if i.tag == Options.STR_SUCCESS])
        nSuccAlt = len([i for i in self.dataset if i.tag == Options.STR_ALT_SUCCESS])
        nFail = len([i for i in self.dataset if i.tag == Options.STR_FAILURE])

        output = "Trajectory Rate, size = %d \n" % len(self.dataset)
        output = output + "nSucc = " + str(nSucc) + "\n"
        if nSuccAlt > 0 :
            output = output + "nAltS = " + str(nSuccAlt)     + "\n"
        output = output +  "  nFail " + str(nFail)  + "\n"
        return output




class Bootstrap():

    def __init__(self, myRates, N=10000, concentration=None, computek1=False, computek1Alt=False):

        self.myRates = myRates

        # # Resample the dataset
        # FD: This is more expensive than strictly required.
        # FD: Note that this computes the CI for kEff().
        b_start_time = time.time()
        self.effectiveRates = list()
        self.effectiveAltRates = list()
        self.logEffectiveRates = list()
        self.N = N

        print "Bootstrapping " + type(myRates).__name__ + ", using " + str(self.N) + " samples.",

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
        print "   ..finished in %.2f sec.\n" % (b_finish_time - b_start_time)

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
            output= "95%% Confidence Interval: [%.3g /M /s, %.3g /M /s] with N=%d" % (low, high, self.N)
            return output
        else :
            return "No successful reactions observed "

# # Concurrent classes start here


class optionsFactory(object):

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


class MergeSimSettings(object):

    RESULTTYPE1 = "FirstStepRate"
    RESULTTYPE2 = "FirstStepLeakRate"
    RESULTTYPE3 = "FirstPassageRate"
    RESULTTYPE4 = "TrajectoryMerge"
    
    DEFAULT_OUTPUT_FILE = "mergesim_output.log"

    debug = False
    resultsType = RESULTTYPE1
    terminationCount = None
    max_trials = 250000000
    timeOut = 24*60*60
    
    saveInterval = 500000
    saveIncrement = saveInterval

    outputFile = DEFAULT_OUTPUT_FILE

    bootstrap = False
    bootstrapN = 0

    def rateFactory(self, dataset=None, endStates= None):
        
        if self.resultsType == self.RESULTTYPE1:
            return FirstStepRate(dataset=dataset)
        if self.resultsType == self.RESULTTYPE2:
            return FirstStepLeakRate(dataset=dataset)
        if self.resultsType == self.RESULTTYPE3:
            return FirstPassageRate(dataset=dataset)
        if self.resultsType == self.RESULTTYPE4:
            return TrajectoryMerge(dataset, endStates)

    def shouldTerminate(self, printFlag, nForwardIn, nReverseIn, timeStart):

        if self.terminationCount == None:
            return True
        else:
            if printFlag:
                print "nForward = %i " % nForwardIn.value
                print "nReverse = %i \n" % nReverseIn.value

       
            if(nForwardIn.value >= self.terminationCount):
                print "Found " + str(nForwardIn.value) + " successful trials, terminating."
                return True

            elif((nForwardIn.value + nReverseIn.value) > self.max_trials):
                print "Simulated " + str(nForwardIn.value + nReverseIn.value) +  " trials, terminating."
                return True
            
            elif( time.time() - timeStart  > self.timeOut):
                print "Runtime exeeded " + str(self.timeOut) + " seconds, terminating."
                return True


            else:
                return False


    def setOutputFile(self,filetitle):
            # sets flags and updates the file title but not more than this.
            # Here, we will wish to append to the file in the future but not normally.
            self.outputFile = filetitle + time.strftime("-%d-%m-%y+%H-%M-%S")+".log"

    def setBootstrap(self, doBootstrap, numIn):
        self.bootstrap = doBootstrap
        self.bootstrapN = numIn

def timeStamp(inTime=None):

    if inTime == None:
        inTime = time.time()

    return str(datetime.datetime.fromtimestamp(inTime).strftime('%Y-%m-%d %H:%M:%S'))



# This class has two modus operandi:
# 1) self.settings.resultsType == self.settings.RESULTTYPE4  (default)
#   In this mode, SimSystem.results and SimSystem.end_states are made available as
#   myMR = MergeSim.run ()  via myMR.datasetand and myMR.endStates (renaming is deliberate)
# 2) self.settings.resultsType != self.settings.RESULTTYPE0:
#    In this mode, the system generates MergeResults objects that compute k1(), kEff() (possibly more)
#    in a multi-processing fashion.
# In addition, users can populate the analysisFactory with an objec to custom
# in-process analysis (for example, computing metrics on traces) 


class MergeSim(object):

    numOfThreads = 2
    seed = 7713147777

    def __init__(self, settings=None):

        self.initializationTime = time.time()
        print("%s%s" % (timeStamp(),
                        "  Starting Multistrand 2.1      (c) 2008-2017 Caltech      "))

        self.factory = optionsFactory
        self.aFactory = None

        if settings == None:
            self.settings = MergeSimSettings()

    # The argument is the count of successfull trials before stopping the simulation
    def setTerminationCriteria(self, terminationCount=25):
        self.settings.terminationCount = terminationCount

    def setFirstStepMode(self):
        self.settings.resultsType = self.settings.RESULTTYPE1

    def setLeakMode(self): # leak mode is a  type of First Step mode
        self.settings.resultsType = self.settings.RESULTTYPE2

    def setPassageMode(self):
        self.settings.resultsType = self.settings.RESULTTYPE3

    def setTrajectoryMode(self):
        self.settings.resultsType = self.settings.RESULTTYPE4

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

    # If the analysis factory is set,
    # then perform an threaded analysis of the returned data.
    # E.g. the analysis factory receives a set of locks and
    # the simulation will call aFactory.doAnalysis(options)
    # once for each thread.
    # This is used really only for one case study, but could be reused in the future
    def setAnaylsisFactory(self, aFactoryIn):

        lockArray = list()
        for i in range(16):
            lockArray.append(multiprocessing.Lock())

        self.aFactory = aFactoryIn
        self.aFactory.lockArray = lockArray

    # reset the multithreading objects
    def clearAnalysisFactory(self):

        if not self.aFactory == None:
            self.aFactory.clear()

    def printTrajectories(self, myOptions):
        # # Print all the trajectories we can find.
        # # Debug function primairly.
        print("Printing trajectory and times: \n")

        trajs = myOptions.full_trajectory
        times = myOptions.full_trajectory_times

        for t, time in zip(trajs, times):

            print (t, "  t=", time, "\n")

        0

    def initialInfo(self):

        myOptions = self.factory.new(777)
        simSys = SimSystem(myOptions)
        simSys.initialInfo()

    def startSimMessage(self):
        
        welcomeMessage = "Using Results Type: " + self.settings.resultsType + "\n"
        # python2 style printing because of an compilation issue with web.py
        welcomeMessage += ''.join(["Computing ", str(
            self.numOfThreads * self.trialsPerThread), " trials, using ", str(self.numOfThreads), " threads .. \n"])

        if not self.settings.terminationCount == None:

            welcomeMessage += " .. and rolling " + str(self.trialsPerThread)
            welcomeMessage += " trajectories per thread until " + str(self.settings.terminationCount) + " successful trials occur. \n"

        return welcomeMessage


    def createOutputMessage(self):
        
        outputString="Trials:              " + str(self.factory.input0) + "\n"
        outputString+="Max Trials:          " + str(self.settings.max_trials) + "\n"
        outputString+="Termination  Count:  " + str(self.settings.terminationCount) + "\n"
        outputString+="Num. of Threads:     " + str(self.numOfThreads) + "\n"
        outputString+="Run Time:            " + str(self.runTime) + "s \n"
        outputString+="Using Results Type:  " + self.settings.resultsType + "\n\n"

        def hiddenPrint(myString):
            
            options = self.factory.new(0)
            myString+="Temperature: " + str(options.temperature) +  "K \n\n"
    
            myString+="Start States\n"
            for x in options.start_state:
                myString+=x.__str__()+"\n\n"
    
            for x in options.stop_conditions:
                myString+=x.__str__()+"\n"


        myProc = multiprocessing.Process(target=hiddenPrint,  args= [outputString])
        myProc.start()
        myProc.join()
        myProc.terminate()

        outputString+="\n"+self.results.__str__() + "\n"
        
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
            print  "Start states:"
            for i in self.factory.new(0).start_state:
                print i
                print "\n"
                
            print "Stop conditions: "
            for i in self.factory.new(0).stop_conditions:
                print i
                print "\n"

             
        myProc = multiprocessing.Process(target=actualPrint,  args= [])
        myProc.start()
        myProc.join()
        myProc.terminate()
                

    def run(self):

        # The input0 is always trials.
        self.trialsPerThread = int(
            math.ceil(float(self.factory.input0) / float(self.numOfThreads)))
        startTime = time.time()

        assert(self.numOfThreads > 0)

        manager = multiprocessing.Manager()

        self.managed_result = manager.list()
        self.managed_endStates = manager.list()
        self.nForward = manager.Value('i', 0)
        self.nReverse = manager.Value('i', 0)

        self.results = self.settings.rateFactory()
        self.endStates = []

        def doSim(myFactory, aFactory, list0, list1, instanceSeed, nForwardIn, nReverseIn):

            myOptions = myFactory.new(instanceSeed)
            myOptions.num_simulations = self.trialsPerThread

            s = SimSystem(myOptions)
            s.start()

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
            return multiprocessing.Process(target=doSim, args=(self.factory, self.aFactory, self.managed_result, self.managed_endStates, instanceSeed, self.nForward, self.nReverse))          


        # this saves the results generated so far as regular Python objects,
        # and clears the concurrent result lists.
        def saveResults():
            
            if self.settings.terminationCount == None:
                # just let the threads join peacefully
               for i in range(self.numOfThreads):
                   procs[i].join()
            else :
                # join all running threads -- the process has 2 seconds to close or it will be terminated.
                for i in range(self.numOfThreads):
                    procs[i].join(timeout = 2)
                    procs[i].terminate()

            self.runTime = (time.time() - startTime)
            print("Done.  %.5f seconds -- now processing results \n" % (time.time() - startTime))


            # Leak - the below is a leak rates object
            # NB: Initialize with a dataset, but we merge with
            # a differrent rates object.
            myFSR = self.settings.rateFactory(self.managed_result, self.managed_endStates)

            self.results.merge(myFSR, deepCopy=True)

            # save the terminal states if we are not in leak mode
            if self.settings.resultsType != self.settings.RESULTTYPE2:
                for endState in self.managed_endStates:
                    self.endStates.append(copy.deepcopy(endState))
            else: #in case of studying leak reaction, expect few hits, and print out info more often.
                print self.results

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
        while True:

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

            time.sleep(0.25)

            # if >500 000 results have been generated, then store
            if (self.nForward.value + self.nReverse.value) > self.settings.saveInterval:
                saveResults()


        saveResults()

        # print final results to the user
        self.results.generateRates()
        
        print self.results


        if self.settings.bootstrap == True:
            self.results.doBootstrap(self.settings.bootstrapN)

        self.writeToFile()    
        return 0





# # The default multistrand object
myMultistrand = MergeSim()
