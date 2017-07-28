"""
Frits Dannenberg, May 17th, 2017.
This module has two parts.
1) Given multistrand results, compute reaction rate constants. 
2) A module for running multistrand concurrently. Results are merged automatically
"""
import operator, os
import time, datetime, math
import traceback

from multistrand.options import Options
import multiprocessing
import numpy as np

from multistrand.system import SimSystem 


MINIMUM_RATE = 1e-36


# # Rate-computation classes start here


# Shared class methods
class basicRate(object):
    
    def log10KEff(self):
        
        return np.log10(self.kEff())


# # Migration rates for first step
class FirstStepRate(basicRate):
    
    
    def __init__(self, dataset=None, concentration=None):
        
        if not concentration == None:        
            self.z = np.float(concentration)
        
        if dataset == None:            
            dataset = []
            print "Setting dataset to empty list"
        

        self.dataset = dataset  # save the dataset for re-sampling and merging results. Also needed for certain rates           
        self.z = concentration
        self.generateRates()


    def generateRates(self):
        
        # Pre-computing some metrics
        self.nForward = sum([i.tag == Options.STR_SUCCESS  for i in self.dataset])
        self.nReverse = sum([i.tag == Options.STR_FAILURE  for i in self.dataset])
        self.nForwardAlt = sum([i.tag == Options.STR_ALT_SUCCESS  for i in self.dataset])
                         
        self.nTotal = len(self.dataset)
        
    
    # return TRUE if enough correct trajectories have been simulated
    # Second argument controls printing during the testing
    def checkTermination(self, datasetIn, printFlag, nForwardIn = None):

        if not nForwardIn == None:
            if printFlag:
                print "nForward = %i \n" % nForwardIn.value
            
            return nForwardIn.value >= self.terminationCount
        
        newRates = FirstStepRate(dataset=datasetIn, concentration=1e-99)

        if printFlag:
            print newRates.shortString()
        
        if newRates.nForward + newRates.nForwardAlt < self.terminationCount:
            return False
        else: 
            print str(newRates)
            return True
        
    
    def sumCollisionForward(self):
        return sum([np.float(i.collision_rate) for i in self.dataset if i.tag == Options.STR_SUCCESS])

    def sumCollisionForwardAlt(self):
        return sum([np.float(i.collision_rate) for i in self.dataset if i.tag == Options.STR_ALT_SUCCESS])
    
    def sumCollisionReverse(self):
        return sum([np.float(i.collision_rate) for i in self.dataset if i.tag == Options.STR_FAILURE])
    
    def weightedForwardUni(self):
        
        mean_collision_forward = np.float(self.sumCollisionForward()) / np.float(self.nForward) 
        weightedForwardUni = sum([ np.float(i.collision_rate) * np.float(i.time) for i in self.dataset if i.tag == Options.STR_SUCCESS]) 
        
        return weightedForwardUni / (mean_collision_forward * np.float(self.nForward))
    
    def weightedReverseUni(self):
        
        if self.nReverse == 0:
            return np.float(0)
        
        mean_collision_reverse = np.float(self.sumCollisionReverse()) / np.float(self.nReverse) 
        weightedReverseUni = sum([np.float(i.collision_rate) * np.float(i.time) for i in self.dataset if i.tag == Options.STR_FAILURE]) 
        
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
            return  np.float(1.0) / self.weightedForwardUni()  
        
            
    def k2Prime(self):
        
        if self.nReverse == 0:
            return MINIMUM_RATE
        else:
            return np.float(1.0) / self.weightedReverseUni()
        
            
    def kEff(self, concentration=None):
         
        if concentration == None:
            z = self.z
        else:
            z = np.float(concentration)
   
        if(self.nForward == 0):
            return MINIMUM_RATE
         
        # expected number of failed collisions
        multiple = (self.k1Prime() / self.k1())   
        
        # the expected collision rates for success and failed collisions respectively
        # this is not equal to k1 k1' because those are weighted by the probability of success 
        # netto, k1 k1' are lower
        
        if self.nForward == 0:
            return  MINIMUM_RATE
        else:
            collForward = self.sumCollisionForward() / np.float(self.nForward)  
        
        if self.nReverse == 0:
            collReverse = MINIMUM_RATE
        else:
            collReverse = self.sumCollisionReverse() / np.float(self.nReverse) 
        
        dTForward = self.weightedForwardUni() + np.float(1.0) / (z * collForward) 
        dTReverse = self.weightedReverseUni() + np.float(1.0) / (z * collReverse) 
         
        dT = dTReverse * multiple + dTForward
        
        return (np.float(1.0) / dT) * (np.float(1.0) / np.float(z))



    def testForTwoStateness(self):
        
        if(self.z == None):
            print("Warning! Attempting to test for two-state behaviour but concentration was not set! \n")
        
        # Test if the failed trajectory and the success trajectory are dominated ( > 10% of total ) by the unimolecular phase
        
        tau_bi_succ = np.float(1) / (self.k1() * self.z)
        tau_bi_fail = np.float(1) / (self.k1Prime() * self.z)
        
        testFail = (tau_bi_fail / self.weightedReverseUni()) < 9
        testSucces = (tau_bi_succ / self.weightedForwardUni()) < 9
        
        
        if(testFail | testSucces):
            
            print(''.join(["Warning! At the chosen concentration, ", str(self.z), " M, the reaction might violate two-state secondary order rate kinetics"]))                         
            print(self)


    def resample(self):
        # returns a new rates object with resampled data
        
        newDataset = []
        N = len(self.dataset)
        
        for i in range(N):
            # generate random between 0 and N-1
            index = int(np.floor(np.random.uniform(high=N)))
            newDataset.append(self.dataset[index])
        
        return FirstStepRate(newDataset, concentration=self.z)
    

    def castToNumpyArray(self):
        

        self.collision_rates = np.array(self.collision_rates)
        self.was_success = np.array(self.was_success)
        self.was_failure = np.array(self.was_failure)
        

    def merge(self, that):
        
        if not self.z == that.z:
            print("Error! Cannot combine results from different concentrations")
            
        # Now merge the existing datastructures with the ones from the new dataset
        self.dataset.append(that.dataset)

        # re-generate basic metrics
        self.generateRates()


        
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
        
        if(self.nForward > 0):
            output += "k_eff    = %.2e  /M /s  \n   " % self.kEff()
        
        
        return output

    def shortString(self):
        
        output = "nForward = " + str(self.nForward) + " \n"
        output += "nReverse = " + str(self.nReverse) + " \n \n"
        
        if self.nForwardAlt > 0 :
            output += "nForwardAlt = " + str(self.nForwardAlt)
        
        if(self.nForward > 0):
            output += "k1       = %.2e  /M /s  \n" % self.k1()
        
        return output


# Like migrationrate, but uses data from first passage time rather than first step mode
class FirstPassageRate(basicRate):
    
    def __init__(self, dataset=None, concentration=None):
        
        if dataset == None:
            
            # set up for resampling 
            self.concentration = concentration
            self.times = []
            self.timeouts = []
            
        else:
            self.process(dataset, concentration)
            self.generateRates()
        
    def process(self, dataset, concentration):
               
        self.concentration = concentration
        self.times = np.array([i.time for i in dataset])
        self.timeouts = [i for i in dataset if not i.tag == Options.STR_SUCCESS]
        
    def castToNumpyArray(self):
        
        self.times = np.array(self.times)
    
    def resample(self):
        
        # avoid resampling time-outs.
        newRates = FirstPassageRate(concentration=self.concentration)
        
        N = len(self.times)
        
        # generate random between 0 and N-1        
        for i in range(N):
        
            index = int(np.floor(np.random.uniform(high=N)))
            newRates.times.append(self.times[index])
        
        newRates.castToNumpyArray()       
        newRates.generateRates()
        
        return newRates
                
    
    def generateRates(self):
        0


    def kEff(self):
        
        if len(self.timeouts) > 0 :
            
            print("# association trajectories did not finish =", str(len(self.timeouts)))
            for i in self.timeouts :
                assert (i.type_name == 'Time')
                assert (i.tag == None)
                assert (i.time >= 10.0)
        
#         print "average completion time = %g seconds at %s" % (np.mean(self.times), concentration_string(self.concentration))
        
        
        mean = np.mean(self.times)
        kEff = np.float(1.0) / (mean * self.concentration)
        
        
        return kEff


class Bootstrap():
    
    
    def __init__(self, myRates, concentration=None, computek1=False, computek1Alt=False):
        
        self.myRates = myRates
        
        self.process(concentration, computek1, computek1Alt)
        
    def process(self, concentration, computek1, computek1Alt):
        
        # # Resample the dataset
        # FD: This is more expensive than strictly required. 
        # FD: Note that this computes the CI for kEff().
        
        self.effectiveRates = list()
        self.effectiveAltRates = list()
        self.logEffectiveRates = list()
        self.N = 400
        
        if not concentration == None:
            self.myRates.z = concentration
        
        for i in range(self.N):
            
            # create a new sample, with replacement
            sample = self.myRates.resample()
            
            # compute k1
            if(computek1):
                rate = float(sample.k1())
            else: 
                rate = float(sample.kEff())

            self.effectiveRates.append(rate)

            if(computek1Alt):
                # print "we did resample here"
                rate = float(sample.k1Alt())
                self.effectiveAltRates.append(rate)
                
            del sample
                
        # sort for percentiles
        self.effectiveRates.sort(cmp=None, key=None, reverse=False)
        self.effectiveAltRates.sort(cmp=None, key=None, reverse=False)

        # Yet to generate log alt rates
        for rate in self.effectiveRates:
            self.logEffectiveRates.append(np.log10(rate))
        
        
    
    def ninetyFivePercentiles(self):
        
        low = self.effectiveRates[int(0.025 * self.N)]
        high = self.effectiveRates[int(0.975 * self.N)]
        
        return low, high

    def ninetyFivePercentilesAlt(self):
        
        try:
            low = self.effectiveAltRates[int(0.025 * self.N)]
            high = self.effectiveAltRates[int(0.975 * self.N)]
            return low, high
        except Exception:
            return 0,0
        
        
    
    def standardDev(self):
        
        return np.std(self.effectiveRates)
    
    def logStd(self):
        
        # Returns standard deviation from the mean of log10 of the original rates
        return np.std(self.logEffectiveRates)

        


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
            output = self.myFunction(self.input0, self.input1, self.input2, self.input3)
        elif self.input5 == None:
            output = self.myFunction(self.input0, self.input1, self.input2, self.input3, self.input4)
        elif self.input6 == None:
            output = self.myFunction(self.input0, self.input1, self.input2, self.input3, self.input4, self.input5)
        else: 
            output = self.myFunction(self.input0, self.input1, self.input2, self.input3, self.input4, self.input5, self.input6)
            
        output.initial_seed = inputSeed
            
        return output





class msSettings(object):
    
    debug = False
    
def timeStamp(inTime=None):
    
    if inTime == None:
        inTime = time.time()
    
    return str(datetime.datetime.fromtimestamp(inTime).strftime('%Y-%m-%d %H:%M:%S'))

class MergeSim(object):

    numOfThreads = 2
    seed = 7713147777;
    

    def __init__(self, settings=None):

        
        self.initializationTime = time.time()
        print("%s%s" % (timeStamp() , "  Starting Multistrand 2.1      (c) 2008-2017 Caltech      "))
                

        self.factory = optionsFactory
        self.aFactory = None
        
        if(settings == None):
            self.settings = msSettings()
            
        # override to enable dynamic simulation length
        self.terminationCriteria = None


    # First argument is typically a FirstStepRate(),
    # Second argument is the count of successfull trials before stopping the simulation 
    def setTerminationCriteria(self, input, terminationCount = 25):
        
        self.terminationCriteria = input
        self.terminationCriteria.terminationCount = terminationCount
        
    def timeSinceStart(self):
        
        print("Time since creating object %.5f seconds" % (time.time() - self.initializationTime))
#                 print("Done.  %.5f seconds" % (time.time() - startTime))
        
    
    def setDebug(self, value):
        
        self.settings.debug = value
        
    
    def setNumOfThreads(self, numOfThreads):
        
        self.numOfThreads = numOfThreads

    def setOptionsFactory(self, optionsFactory):  # this can be re-done using an args[] obj.
        
        self.factory = optionsFactory

    def setOptionsFactory1(self, myFun, put0):
        
        self.factory = optionsFactory(myFun, put0, None, None, None, None, None, None) 
 
    def setOptionsFactory2(self, myFun, put0, put1):
        
        self.factory = optionsFactory(myFun, put0, put1, None, None, None, None, None) 

    def setOptionsFactory3(self, myFun, put0, put1, put2):
        
        self.factory = optionsFactory(myFun, put0, put1, put2, None, None, None, None) 

    def setOptionsFactory4(self, myFun, put0, put1, put2, put3):
        
        self.factory = optionsFactory(myFun, put0, put1, put2, put3, None, None, None)
    
    def setOptionsFactory5(self, myFun, put0, put1, put2, put3, put4):
        
        self.factory = optionsFactory(myFun, put0, put1, put2, put3, put4, None, None)  
    
    def setOptionsFactory6(self, myFun, put0, put1, put2, put3, put4, put5):
        
        self.factory = optionsFactory(myFun, put0, put1, put2, put3, put4, put5, None)

    def setOptionsFactory7(self, myFun, put0, put1, put2, put3, put4, put5, put6):
        
        self.factory = optionsFactory(myFun, put0, put1, put2, put3, put4, put5, put6)

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
        

    def numOfSuccesTrials(self):
        
        count = 0
        
        for result in self.results:
            if result.tag == Options.STR_SUCCESS:
                count += 1

        return count

    # reset the multithreading objects
    def clear(self):
        
        if not self.aFactory == None:
            self.aFactory.clear()
    
    def printTrajectories(self, myOptions):
        # # Print all the trajectories we can find. 
        # # Debug function primairly. 
        
        print("Printing trajectory and times: \n")
        
        
        trajs = myOptions.full_trajectory 
        times = myOptions.full_trajectory_times
        
        for t, time in zip (trajs, times):
            
            print (t, "  t=", time, "\n")
        
        0
    
    def initialInfo(self):
        
        myOptions = self.factory.new(777)
        simSys = SimSystem(myOptions)
        simSys.initialInfo()


    def startSimMessage(self):
        
        # python2 style printing because of an compilation issue with web.py        
        welcomeMessage = ''.join(["Computing ", str(self.numOfThreads * self.trialsPerThread) , " trials, using " , str(self.numOfThreads), " threads .. \n" ])
        
        if not self.terminationCriteria == None:
            
            welcomeMessage += " .. and rolling " + str(self.trialsPerThread)
            welcomeMessage += " trajectories per thread until " + str(self.terminationCriteria.terminationCount) + " successful trials occur. \n"                    
        return welcomeMessage


    def run(self):
        
        
        # The input0 is always trails.
        self.trialsPerThread = int(math.ceil(float(self.factory.input0) / float(self.numOfThreads)))
        startTime = time.time()

        
        try:
            for i in self.factory.new(0).stop_conditions:
                print i
        except Exception:
                print "No stop conditions defined"
        
        

        def doSim(myFactory, aFactory, list0, list1, instanceSeed, nForwardIn):
            
            
            myOptions = myFactory.new(instanceSeed)
            myOptions.num_simulations = self.trialsPerThread

            s = SimSystem(myOptions)
            s.start()
            
            
            if myOptions.simulation_mode == Options.firstStep :
                
                myFSR = FirstStepRate(myOptions.interface.results)
                nForwardIn.value +=  myFSR.nForward + myFSR.nForwardAlt
                
            for result in myOptions.interface.results:
                
                list0.append(result)   
                 
            for endState in myOptions.interface.end_states:
                
                list1.append(endState)
                            
            if self.settings.debug:
                
                self.printTrajectories(myOptions)
                
            if not(aFactory == None):
                 
                aFactory.doAnalysis(myOptions)

        
        assert(self.numOfThreads > 0)

        manager = multiprocessing.Manager()
        
        self.results = manager.list()
        self.endStates = manager.list()
        self.nForward = manager.Value('i',0)
        

        def getSimulation():
            
            instanceSeed =  self.seed + i * 3 * 5 * 19 + (time.time() * 10000) % (math.pow(2, 32) - 1)
            return multiprocessing.Process(target=doSim, args=(self.factory, self.aFactory, self.results, self.endStates, instanceSeed, self.nForward))
        
        def shouldTerminate(printFlag):
            return (self.terminationCriteria == None) or self.terminationCriteria.checkTermination(self.results, printFlag, self.nForward)
        

        # start the initial bulk
        print(self.startSimMessage()) 
        procs = []
         
        for i in range(self.numOfThreads):
                       
            p = getSimulation()
            procs.append(p)
            p.start()

        sleepTime = 0.3
        printFlag = False

        # check for stop conditions, restart sims if needed
        while True:
            
            if shouldTerminate(printFlag):
                break 

            printFlag = False
            
            # find and re-start finished threads             
            for i in range(self.numOfThreads):
                 
                if not procs[i].is_alive():
                     
                    procs[i].join()
                    procs[i].terminate()
                     
                    procs[i] = getSimulation()
                    procs[i].start()
                    
                    printFlag = True
                    
             
            time.sleep(sleepTime)
            if(sleepTime < 3.0):
                sleepTime = sleepTime * 1.3
    
             
        # join all running threads
        for i in range(self.numOfThreads):
            procs[i].join()
        
        self.results = list(self.results)
        self.endStates = list(self.endStates)
        
        myOptions = self.factory.new(77) #only to get concentration
        if myOptions.join_concentration == None:
            newRates = FirstStepRate(dataset=self.results, concentration=50e-9)
        else:
            newRates = FirstStepRate(dataset=self.results, concentration=myOptions.join_concentration)
            
        print str(newRates)
        
        self.runTime = (time.time() - startTime)
        print("Done.  %.5f seconds \n" % (time.time() - startTime))
        
        return 0



# # The default multistrand object
myMultistrand = MergeSim()
        
        
          
          
          
          
          
          
          
          
          
          
          
          
        
