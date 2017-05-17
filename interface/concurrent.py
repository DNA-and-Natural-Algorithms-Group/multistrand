"""
Frits Dannenberg, May 17th, 2017.
This module has two parts.
1) Given multistrand results, compute reaction rate constants. 
2) A module for running multistrand concurrently.
"""

import numpy as np
from multistrand.options import Options



MINIMUM_RATE = 1e-18;



# Shared class methods
class basicRate(object):
    
    def log10KEff(self):
        
        return np.log10(self.kEff())


## Migration rates for first step
class migrationRate(basicRate):
    
    
    def __init__(self, dataset=None, concentration=None):
        
        
        if not concentration == None:        
            # empty constructor for resampling code
            self.z = np.float(concentration)

        
        if dataset == None:
        
            
            self.collision_rates = []
            self.was_success = []
            self.was_failure = []
            self.time = []
            
            self.generateRates()

        else:
        
            # process the result on object creation
            self.process(dataset, concentration)
            self.generateRates()
      
    
    def terminate(self, dataset):
        
        newRates = migrationRate(dataset=dataset, concentration=1e-99)
        
        print(str(newRates))
        
        if(newRates.nForward <= 25):
            return False
        else: 
            return True
        
        
    
    def process(self, dataset, concentration):
        
        # concentration not actually used in generateRates
        self.z = concentration
        
        # Now determine the reaction model parameters from the simulation results.  (Simplified from hybridization_first_step_mode.py.)    
        
        # if there is no collision_rate set, then either something is wrong or we're not using first step mode. Default to 
        # a standard collisionrate so we can at least run the code TODO
         
        self.collision_rates = np.array([ i.collision_rate for i in dataset])
        
        self.was_success = np.array([1 if i.tag == Options.STR_SUCCESS else 0 for i in dataset])
        self.was_failure = np.array([1 if i.tag == Options.STR_FAILURE else 0 for i in dataset])
        
        # save unprocessed times for re-sample
        self.time = np.array([ i.time for i in dataset])
        

    def generateRates(self):
        
        
        self.forward_times = np.array([self.time[i] for i in range(len(self.was_success)) if self.was_success[i] == 1])
        self.reverse_times = np.array([self.time[i] for i in range(len(self.was_failure)) if self.was_failure[i] == 1])
       
        self.collision_forward = np.array([self.collision_rates[i] for i in range(len(self.was_success)) if self.was_success[i] == 1])
        self.collision_reverse = np.array([self.collision_rates[i] for i in range(len(self.was_failure)) if self.was_failure[i] == 1])
        
        self.nForward = sum(self.was_success)
        self.nReverse = sum(self.was_failure)
        self.nTotal = self.nForward + self.nReverse

    # FD: this is the concentration-independent rate
    def k1(self):

        # Erik's notes:   
        # Basic calculations from Joseph's PhD thesis.  (See also threewaybm_first_step_mode.py.)
        # The equations there were derived for runs stating with exact microstates.
        # Here, as we are using Boltzmann sampling, which requires that the mean collision rate 
        # be calculated separately for forward and reverse cases, since conditioning on success
        # might alter the distribution of secondary structures, and thus the collision rates.
     
        forward_kcoll = np.mean(self.collision_forward)
        
        if not self.nTotal==0:
            prob = np.float(self.nForward) / np.float(self.nTotal)
        else:
            return MINIMUM_RATE
#             raise Warning("Could not find the results, nTotal = 0")

        k1 = prob * forward_kcoll  
        
        return k1
    
    # FD: this is the concentration-dependent rate
    def kEff(self, concentration=None):
        
        if not concentration == None:
            z = np.float(concentration)
        else:
            z = self.z

        if(self.nForward == 0):
            return MINIMUM_RATE
    
        dTsuccess_uni = np.mean(self.forward_times)
        dTfail_uni = np.mean(self.reverse_times)
        
        kcollision = np.mean(self.collision_rates)
        
        kcollision_foward = np.mean(self.collision_forward)
        kcollision_reverse = np.mean(self.collision_reverse)        
        
        
        multiple = float(self.nReverse) / float(self.nForward)

        dTcoll = np.float(1.0) / ((kcollision) * z)

        if(True):  # toggle if use specific collision rates or not
            dTcoll_forward = np.float(1.0) / ((kcollision_foward) * z)
            dTcoll_reverse = np.float(1.0) / ((kcollision_reverse) * z)
        else:
            dTcoll_forward = dTcoll
            dTcoll_reverse = dTcoll
       
        
        dTfail = dTcoll_reverse + dTfail_uni
        dTforward = dTcoll_forward + dTsuccess_uni
        
        dTcorrect = dTfail * multiple + dTforward
        
        kEff = (1 / dTcorrect) * (1 / z)
        
        return kEff



    def resample(self):
        
        # returns a new rates object with resampled data
        
        newRates = migrationRate(concentration=self.z)
        
        N = len(self.collision_rates)
        
        # # generate random between 0 and N-1
        
        for i in range(N):
        
            index = int(np.floor(np.random.uniform(high=N)))
            
            newRates.collision_rates.append(self.collision_rates[index])
            newRates.was_success.append(self.was_success[index])
            newRates.was_failure.append(self.was_failure[index])
            newRates.time.append(self.time[index])

        
        newRates.castToNumpyArray()       
        newRates.generateRates()
        
        return newRates
    

    def castToNumpyArray(self):
        

        self.collision_rates = np.array(self.collision_rates)
        self.was_success = np.array(self.was_success)
        self.was_failure = np.array(self.was_failure)
        self.time = np.array(self.time)
        

    def merge(self, that):
        
        if not self.z == that.z:
            print("Error! Cannot combine results from different concentrations")
            
        # Now merge the existing datastructures with the ones from the new dataset
        self.collision_rates = np.concatenate([self.collision_rates, that.collision_rates])
        self.was_success = np.concatenate([self.was_success, that.was_success])
        self.was_failure = np.concatenate([self.was_failure, that.was_failure])
        self.time = np.concatenate([self.time, that.time])

        # re-generate the rates
        self.generateRates()


        
    # # override toString
    def __str__(self):
        
        output = "nForward = " + str(self.nForward) + " \n"
        output += "nReverse = " + str(self.nReverse) + " \n"
        output += "k1 = " + str(self.k1()) + " /M/s   \n"
        output += "k_eff = " + str(self.kEff()) + " /M/s   \n"
        output += ", concentration = " + str(self.z) + " M  \n "
        
        return output



# Like migrationrate, but uses data from first passage time rather than first step mode
class migrationRatePassage(basicRate):
    
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
        newRates = migrationRatePassage(concentration=self.concentration)
        
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
            print "# association trajectories did not finish =", len(self.timeouts)
            for i in self.timeouts :
                assert (i.type_name == 'Time')
                assert (i.tag == None)
                assert (i.time >= 10.0)
        
#         print "average completion time = %g seconds at %s" % (np.mean(self.times), concentration_string(self.concentration))
        
        
        mean = np.mean(self.times)
        kEff = np.float(1.0) / (mean * self.concentration)
        
        
        return kEff


class bootstrap():
    
    
    def __init__(self, myRates, concentration=None):
        
        self.myRates = myRates
        
        self.process(concentration)
        
    def process(self, concentration):
        
        # # Resample the dataset
        
        self.effectiveRates = list()
        self.logEffectiveRates = list()
        self.N = 1000
        
        if not concentration == None:
            self.myRates.z = concentration
        
        for i in range(self.N):
            
            sample = self.myRates.resample()
            rate = float(sample.kEff())
            self.effectiveRates.append(rate)
            
            del sample
                
        self.effectiveRates.sort(cmp=None, key=None, reverse=False)
        
        for rate in self.effectiveRates:
            self.logEffectiveRates.append(np.log10(rate))
        
        
    
    def ninetyFivePercentiles(self):
        
        low = self.effectiveRates[int(0.025 * self.N)]
        high = self.effectiveRates[int(0.975 * self.N)]
        
        return low, high
    
    def standardDev(self):
        
        return np.std(self.effectiveRates)
    
    def logStd(self):
        
        # Returns standard deviation from the mean of log10 of the original rates
        return np.std(self.logEffectiveRates)

        
         
        
        
          