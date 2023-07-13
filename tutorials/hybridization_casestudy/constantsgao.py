import operator
import multiprocessing
import numpy as np

from multistrand.options import Options, Literals

colors = ['blue', 'red', 'cyan', 'magenta', 'green', 'k', 'darkblue', 'darkred', 'darkcyan', 'darkmagenta', 'darkgreen']


goa2006_P0 = 'GTTGTCAAGATGCTACCGTTCAGAG'
goa2006_P3 = 'AGATCAGTGCGTCTGTACTAGCAGT'
goa2006_P4 = 'AGATCAGTGCGTCTGTACTAGCACA'

STR_ALL = "ALL"

def setSaltGao2006(o):
    
        o.sodium = 0.5;
        o.magnesium = 0.0;
 
 
 
 
 
 
 
 
 
''' The following are required for the case3 case study file '''
    




# hard-coded locks
# # this is a manager.dict, that is dict, lock is multithreading.lock
def mergeDictAdd (this, that, lock):
    
    with lock:    
        for key, value in that.items():
            
            if(key in this):
                this[key] += value
            else:
                this[key] = value
            
def mergeDictBinary(this, that, lock):
    
    # # If an entry in that is non-zero, add +1 to this using the same key
    with lock:
        for key, value in that.items():
            
            if(value > 0):
            
                if(key in this):
                    this[key] += 1
                else:
                    this[key] = 1


def mergeDict (this, that, lock):
    
    with lock:    
        for key, value in that.items():
            
            if not(key in this):
        
                this[key] = value

def mergeList (this, that, lock):
    
    with lock:    
        for entry in that:
            
            this.append(entry)
        

def mergeCount(mValue, counter, lock):
        
    with lock:
        mValue.value = mValue.value + counter
        



class position(object):
    
    def __init__(self, left, right):
        
        self.posX = left
        self.posY = right

    def __hash__(self):       
        return (self.posX + 1000 * self.posY)  # works for complexes < 1000 nt
    
    def __eq__(self, other):
        return (self.posX == other.posX) and (self.posY == other.posY)
    
    def __str__(self):
        return (str(self.posX) + ' ' + str(self.posY))
    
    
    def toString(self):
        return (str(self.posX) + ' ' + str(self.posY))
    

class analysisResult(object):

    
    def  __init__(self, selector):
  
        self.clause = selector
    

        manager = multiprocessing.Manager()
    
        self.pathCount = multiprocessing.Value('i', 0)
        self.posDict = manager.dict()  # This collects pure time
        self.countDict = manager.dict()
        self.structDict2 = [ manager.dict() for x in range(30 * 30)] 
        self.binaryDict = manager.dict()  # This lists by how many paths each location is visited.
         
        self.pathProps = manager.list()  # Collects if the strands initial binding is aligned or not.
        
        self.initialTrajectory = manager.list()  # This is the first trajectory we encouter, save it to plot it later.
    
            
    def clear(self):
        
        self.pathCount.value = 0 
        self.posDict.clear()
        self.countDict.clear()
        del self.pathProps[:]
        
        manager = multiprocessing.Manager()
        self.structDict2 = [ manager.dict() for x in range(30 * 30)] 
    
    
    def setInitialTrajectory(self, traj, lock):
        
        with lock:
        
            if not traj == None:
            
                self.initialTrajectory.append(traj)



def getTubeStruct(traj, index):
    
    states = traj[index]
    structs = []
    for state in states: structs += [ state[4] ]  # similarly extract the secondary structures for each complex
    tubeStruct = ' '.join(structs)  # give the dot-paren secondary structure for the whole test tube
    
    return tubeStruct



class pathProperties(object):
    
    def __init__(self):
        self.aligned = None
        self.tag = None



class analysisFactory(object):
       
       
    def __init__(self, inSeq, inCutOff):
        
        self.mySeq = inSeq
        self.cutOff = inCutOff

        self.selectors = [STR_ALL, Literals.success, Literals.failure]

        self.result0 = analysisResult(self.selectors[0])
        self.result1 = analysisResult(self.selectors[1])
        self.result2 = analysisResult(self.selectors[2])


    
    
    def clear(self):
        
        self.result0.clear()
        self.result1.clear()
        self.result2.clear()
             
    def processStructs(self, array, index):
        
        # the byte array is of form LLL+RRR where length(LLL) = index
        left1 = 0;
        right1 = 0;
        left2 = 0;
        right2 = 0;
        
        
        for i in range(index):
            if (array[i] == 40):
                left1 = left1 + 1
            if (array[i] == 41):
                right1 = right1 + 1
        
        offset = index + 1 
        for i in range(index):
            if (array[offset + i] == 40):
                left2 = left2 + 1
            if (array[offset + i] == 41):
                right2 = right2 + 1
    
    
        return position(left1 - right1, right1 + left2)


    def processStructsString(self, string, index):
        
        # the byte array is of form LLL+RRR where length(LLL) = index
        left1 = 0;
        right1 = 0;
        left2 = 0;
        right2 = 0;
        
        
        for i in range(index):
            if (string[i] == '('):
                left1 = left1 + 1
            if (string[i] == ')'):
                right1 = right1 + 1
        
        offset = index + 1 
        for i in range(index):
            if (string[offset + i] == '('):
                left2 = left2 + 1
            if (string[offset + i] == ')'):
                right2 = right2 + 1
    
    
        return position(left1 - right1, right1 + left2)

    
    def selectTrajectories(self, clause, results, trajectories, times):
        
        # # The trajectories and end times are not seperated by the results.tag. This is a helper function to return only those trajectories and times for which
        # # the results.tag == clause
        outTraj = list()
        outTime = list()
        
        # # Refactoring: outTraj, outTime are now lists of lists, each entry is a seperate trajectory.
        pathCounter = -1
        pathSelected = False  
        selectedCounter = 0
        
        myTraj = None
        myTime = None
        
        def savePaths(trajs, times):
            if len(trajs) > 0:
                outTraj.append(trajs)
                outTime.append(times)
        
        for traj, time in zip(trajectories, times):
            
            if(time == 0.0):
                
                if not myTraj == None:
                    savePaths(myTraj, myTime)
                
                # # a new path is triggered
                myTraj = list()
                myTime = list()
                
                pathCounter += 1
                pathSelected = (results[pathCounter].tag == clause)
                
                if(clause == "ALL"):
                    pathSelected = True
                
                if(pathSelected):
                        selectedCounter += 1


            if(pathSelected):
                
                myTraj.append(traj)
                myTime.append(time)
                
        savePaths(myTraj, myTime)
                        
        return outTraj, outTime, selectedCounter
    

    
    def checkTag(self, trajectory, prop):

        length = len(self.mySeq)
        succesStruct = length * '(' + "+" + length * ')'
        tubeStruct = getTubeStruct(trajectory, -1)  # Python way of accessing last element

        
        if tubeStruct == succesStruct:
            prop.tag = Literals.success
        else:
            prop.tag = Literals.failure
                
        
    
    def checkAligned(self, trajectory, prop):
        
        tubeStruct = getTubeStruct(trajectory, 0)
        length = len(self.mySeq)

        
        def findAbnormal(string):
            # # Input dot-paren structure, with one unclosed bracket. 
            # # Return the location of the unclosed bracket.
            
            myStack = list()
            
            for i in range(len(string)):
            
                if(string[i] == '('):
                    myStack.append(i)
                
                if(string[i] == ')'):
                    if(len(myStack) == 0):
                        return i
                    else:
                        myStack.pop()
            
#             assert len(myStack) == 1
            return myStack.pop()
            
                    
        left = findAbnormal(tubeStruct[0:length])
        right = findAbnormal(tubeStruct[(length + 1):(2 * length + 1)])
        
        prop.aligned = (left + right) == (length - 1)
        
        
        

    def processSingleTraj(self, trajectory, times, posDict, countDict, structDict2):
        
        myLength = len(self.mySeq)
        
        oldTime = 0.0

     
        for i in range(len(trajectory) - 1):  # go through each output microstate of the trajectory
            
            time = times[i + 1]  # time at which this microstate is exited (does not exist for final state)
            tubeStruct = getTubeStruct(trajectory, i)
                        
            # Re-doing this without bytearray
            position = self.processStructsString(tubeStruct, myLength)            
            timeInState = np.max([time - oldTime, 0.0])
            oldTime = time

            # # Save the time of the position       
            if(position.posY < self.cutOff):
                
                if (position in posDict):
                    posDict[position] += timeInState 
                    countDict[position] += 1
                else:
                    posDict[position] = timeInState
                    countDict[position] = 1
                    
                    
                # Save the count of the position-struct:
                index = position.posX * 30 + position.posY
        
                myDict = structDict2[index];
                if not(tubeStruct in myDict):
                    myDict[tubeStruct] = 1
                else:
                    myDict[tubeStruct] += 1
                    
        

        
                                
    def analyzeTrajectory(self, trajectories, times):
    
        # # Refactoring: trajs and times aren now lists-of-lists.
        posDict = dict()
        countDict = dict()  # counts the number of hits for each position.
        structDict2 = [dict() for x in range(30 * 30)]  # counts the number of structs per position.
        
        pathProps = list()
        
        
        for traj, time in zip(trajectories, times):
            
            self.processSingleTraj(traj, time, posDict, countDict, structDict2)
            
            prop = pathProperties()             
            self.checkAligned(traj, prop)
            self.checkTag(traj, prop)
            pathProps.append(prop)
            
            
                
        # # remove all but the top 100 entries per position
        for i in range(len(structDict2)):
            myDict = structDict2[i]
            structDict2[i] = dict(sorted(iter(myDict.items()), key=operator.itemgetter(1), reverse=True)[:100])

                
        return posDict, countDict, pathProps, structDict2

    def parseFirstTrajectory(self, trajectories, times):
        
        # simply selects the first trajectory, parses it into positions, returns 
        # the list of positions
    
        if len(trajectories) > 0 :
    
            positions = list()
     
            traj = trajectories[0]
            timeList = times[0]
           
            for i in range(len(traj) - 2):  # go through each output microstate of the trajectory
                
                time = timeList[i + 1]  # time at which this microstate is exited (does not exist for final state)
                tubeStruct = getTubeStruct(traj, i)
                            
                # Re-doing this without bytearray
                position = self.processStructsString(tubeStruct, len(self.mySeq))            
                positions.append(position)
#                 print("First trajectory position is ", str(position))
    
            return positions

        else:
            
            return None




    #FD: Analysis factory is only used for the likelihood plots in the Multistrand 2.0 casestudy.
    def doAnalysis(self, myOptions):
    
    
        myAnalysis = analysisFactory(self.mySeq, self.cutOff)
    
        def processAndMergeDicts(clause, myOptions, myAnalysis, analysisResult, lockArray):
        
            
            traj, times, pathCount = myAnalysis.selectTrajectories(clause, myOptions.interface.results, myOptions.full_trajectory, myOptions.full_trajectory_times)
        
            posDict, countDict, pathProps, structDict2 = myAnalysis.analyzeTrajectory(traj, times)   
            firstTraj = myAnalysis.parseFirstTrajectory(traj, times)
            analysisResult.setInitialTrajectory(firstTraj, lockArray[1])
            
            mergeDictAdd(analysisResult.posDict, posDict, lockArray[0])
            mergeDict(analysisResult.countDict, countDict, lockArray[2])
            mergeDictBinary(analysisResult.binaryDict, countDict, lockArray[3])
            
            
            mergeList(analysisResult.pathProps, pathProps, lockArray[4])
            mergeCount(analysisResult.pathCount, pathCount, lockArray[5])
            
            for i in range(len(structDict2)):
                if not structDict2[i] == None:
                    mergeDictAdd(analysisResult.structDict2[i], structDict2[i], lockArray[6])
                
            
            
    
        processAndMergeDicts(myAnalysis.selectors[0], myOptions, myAnalysis, self.result0, self.lockArray)
        processAndMergeDicts(myAnalysis.selectors[1], myOptions, myAnalysis, self.result1, self.lockArray)
        processAndMergeDicts(myAnalysis.selectors[2], myOptions, myAnalysis, self.result2, self.lockArray)
        




 
 