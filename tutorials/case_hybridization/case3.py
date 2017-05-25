# Frits Dannenberg, Caltech, 2017.

"""
This script computes density plots of the hybridization reaction. Use this as 
<Number of Threads> <Number of Trajectories> <Sequence: P0, P3, P4>
Example:

14 80 P0, or:
2 20 P3

Some routines are ommited for testing purposes, use:
2 24 test2
"""

from multistrand.options import Options

from multistrand.concurrent import myMultistrand, MsMulti, position
from multistrand.utils import standardFileName
from multistrand.experiment import colors, setSaltGao2006,  hybridization, standardOptions,  colors, goa2006_P0, goa2006_P3, goa2006_P4



from matplotlib.collections import LineCollection
import matplotlib.pylab as plt
import numpy as np
import operator, sys

SCRIPT_DIR = "Hybridization_F3"
TEMPERATURE = 20.0
ATIME_OUT = 0.0004


class hybridizationSimSettings(object):
    
    def __init__(self, mySeq, cutOff, trials):
    
        self.mySeq = mySeq
        self.cutOff = cutOff
        self.trials = trials
        
        
    def __str__(self):
        
        output = "mySeq=" + self.mySeq
        output += output + str(self.cutOff)
        
        return output

    
def getOptions(trials , settings): #start_complex_top, start_complex_bot, success_stop_condition, failed_stop_condition
      
    options = standardOptions("First Step", TEMPERATURE, settings.trials, ATIME_OUT)

    hybridization(options, settings.mySeq, settings.trials)
    setSaltGao2006(options)

    options.output_interval = 1  # print every state, ever
   
    return options


def first_step_simulation(multistrandObject, settings):
    
    print ("Running first step mode simulations for %s (with Boltzmann sampling)..." % (settings.mySeq))
    
        
    multistrandObject.setOptionsFactory2(getOptions, settings.trials, settings)
    multistrandObject.setAnaylsisFactory(settings.mySeq, settings.cutOff)        
    multistrandObject.run()  
    
    return 0
    
  

def doPosPlotsPrimer(analysisResult, settings, extraTitle):
    
    posDict = analysisResult.posDict
    selectedCount = analysisResult.pathCount.value
    trajectoryIn = analysisResult.initialTrajectory
    
        
    doPosPlots(posDict, settings, extraTitle, selectedCount, "M1", analysisResult.structDict2, trajectory=trajectoryIn)
    

def standardTitle(mySeq, extraTitle, selectedCount, runs):
    
    title = "Hybridization alignment -- " + mySeq + extraTitle + "\n Found " 
    title += str(selectedCount) + "/"      
    if not runs == None:
        title += str(runs) 
    title += " trials." + " T=" + str(TEMPERATURE) + "C, " + str(ATIME_OUT) + "s"
    
    return title


def estimateSuccessProbability(popularStructure, settings):  
    
    # Do a quick first step simulation and see how many end up succeeding.
    # Use a new multistrand object so we don't lose the old results 
    
    newMultistrand = msMulti()
    newMultistrand.setNumOfThreads(myMultistrand.numOfThreads)
 
    oldTrials = settings.trials
 
    settings.trials = int(settings.trials / 3.0)
    settings.initialStructure = popularStructure
  
 
    first_step_simulation(newMultistrand, settings)
    settings.trials = oldTrials
    
    output = np.float(newMultistrand.numOfSuccesTrials()) / np.float(newMultistrand.trialsPerThread * newMultistrand.numOfThreads)
    
    newMultistrand.clear()
    
    return output
    
    

def computeWinProb(f2, pos, structDict2, settings):
    
    # # grab the most frequent structure in structDict2, 
    # # simulate it a few times, and return the expected probabiltiy to hybridize
    
    structs = dict(structDict2[pos.posX * 30 + pos.posY])
#     
    mostPopular = sorted(structs.iteritems(), key=operator.itemgetter(0), reverse=True)[:1]
       
    popularStructure = mostPopular[0][0]
    
    f2.write(str(pos) + " " + popularStructure + "\n")

    return estimateSuccessProbability(popularStructure, settings)
    
def plotMostFrequentStructure(posDict, length):
    
    
            # # Add a spline interpolation of the most popular state per column
        mostFreq = list()
        
        for i in range(length):
            mostFreq.append((-1, -99));
            
        for key, val in posDict.iteritems():
            
            x, y = key.posX, key.posY
            currMax = mostFreq[x][1]            
            
            if val > currMax:
                mostFreq[x] = [y, val]
        
            
        # # determine highest x value for which there are observations
        maxX = max([pos.posX for pos in posDict])
                
        xVals = range(1, maxX + 1)
        yVals = [o[0] for o in mostFreq]
        yVals = yVals[1:maxX + 1]
       
        # Somehow, a line segment doesn't disrupt the axis.       
        segments = list()
        
        for i in range(maxX - 1):
            
            point = [ (xVals[i], yVals[i ]), (xVals[i + 1], yVals[i + 1])]
            segments.append(point)
        
        lc = LineCollection(segments, color="orange", linewidths=6, alpha=0.75)
        plt.gca().add_collection(lc)
        
        return mostFreq, maxX

def plotFirstTrajectory(filename, trajectories, length):
# Add an interpolation of the first four trajectories
    
    f = open(filename + "-trajectories.txt", "w")
    
    upTo = 3
    
    if len(trajectories) < upTo:
        upTo = len(trajectories)
        
        

    for j in range(upTo):

        trajectory = trajectories[j]
        
        trajLength = len(trajectory)
                
        # Somehow, a line segment doesn't disrupt the axis.               
        f.write(str(colors[j]) + " " + "TrajectoryLength= " + str(trajLength) + " \n") 
        
        segments = list()
        
        for i in range(trajLength - 1):
            
            curr = trajectory[i]
            following = trajectory[i + 1]
            
            point = [ (curr.posX, curr.posY), (following.posX, following.posY)]
            segments.append(point)
        
        lc = LineCollection(segments, color=colors[j], alpha=0.65 - 0.12 * j, linewidths=2 + 0.35 * j)
        
        plt.gca().add_collection(lc)
    
    f.close()

def doPosPlots(posDict, settings, extraTitle, selectedCount, extraSettings, structDict2=None, trajectory=None):
    
    fileName = standardFileName(SCRIPT_DIR, settings.mySeq, extraTitle, settings.trials) 
    
    
    length = len(settings.mySeq) + 1
    
    # Make a grid...
    nrows, ncols = settings.cutOff, length 
    image = np.zeros(nrows * ncols)

    for i in range(nrows):
        for j in range(ncols):
            image[i + nrows * j] = -990.0
            

    myMax = -5.0
    myMin = myMax - 7.0        
    corrector = 1.0;
    
    
    goodPosDict = dict(posDict)
    
    if ("winprob" in extraSettings):
        myMax = np.log10(np.float(2.0))
        myMin = myMax - 1.5        
    
    if ("binaryProb" in extraSettings):
        myMax = 0.1
        myMin = -2.0        

    if ("M1" in extraSettings):
        extraTitle = "-ModeOne-" + extraTitle    
        corrector = settings.trials

        
    for pos, val in goodPosDict.iteritems():
 
        value = 0.0
 
        if(val > 0.0):
            value = np.log10(val / corrector)
        else:
            value = -99.9
            
        image[pos.posX + ncols * pos.posY] = np.power(value, 1.0)
   
    # Reshape things into a grid
    image = image.reshape((nrows, ncols))
    
    row_labels = range(nrows)
    col_labels = range(ncols)
    plt.matshow(image, cmap=plt.cm.gray_r, vmin=myMin, vmax=myMax)   
    
    if("M1" in extraSettings):
        
#         if not "test" in extraTitle:
        mostFreq, maxX = plotMostFrequentStructure(goodPosDict, length)
        
        if not (trajectory == None):
            plotFirstTrajectory(fileName, trajectory, length)
    
    plt.title(standardTitle(settings.mySeq, extraTitle, selectedCount, settings.trials))
      
    ax = plt.gca()

    ax.xaxis.set_ticks_position('bottom')
    ax.invert_yaxis()
    
    ax.set_ylabel('Basepairs within strands')    
    ax.set_xlabel('Basepairs between strands')
    
    plt.xticks(range(ncols), col_labels)
    plt.yticks(range(nrows), row_labels)
    plt.colorbar(shrink=0.75)  # orientation='Horizontal'
    
    plt.savefig(fileName + '.pdf')
    plt.close()
    
    # Now generate a new plot for probability of success of the most common structures
    
    if("M1" in extraSettings and (Options.STR_SUCCESS in extraTitle or Options.STR_FAILURE in extraTitle)  and  not "test" in extraTitle) :           
            
        fileName = standardFileName(SCRIPT_DIR, settings.mySeq, extraTitle + "", settings.trials)

        winProb = list()
        plotRange = range(1, maxX)

        
        f2 = open(fileName + "-mostPopularStructs.txt", 'w')
        
        for x in plotRange:
            
            y = mostFreq[x][0]
            myPos = position(x, y)
            
            winProb.append(100.0 * computeWinProb(f2, myPos, structDict2, settings))

        f2.close()
        
        fig = plt.figure()

        plt.fill_between(plotRange, 0.0, winProb, facecolor='orange', edgecolor='orange', alpha=0.5)
        plt.savefig(fileName + '-structSample.pdf')

    
        
        


def writeStructFile(analysisResult, settings, extraTitle):
    
    fileName = standardFileName(SCRIPT_DIR, settings.mySeq, extraTitle, settings.trials)
    f = open(fileName + "-struct.txt", 'w')    
     
    goodDict = dict(analysisResult.posDict)
     
    for pos, val in goodDict.iteritems():
        
        output = "Pos = " + pos.toString() + " Freq= " + str(val) + "\n"     
        f.write(output)
  
    f.close()
    
    # also write the position-struct file
    f = open(fileName + "-posStruct.txt", 'w')    

    for i in range(len(analysisResult.structDict2)):
        
        goodDict = dict(analysisResult.structDict2[i])
        
        # only print the top 20 of structures found
        goodDict = dict(sorted(goodDict.iteritems(), key=operator.itemgetter(1), reverse=True)[:20])
         
        pX = np.int((np.floor(i / 30)))
        pY = np.int(i % 30)
         
        for key, val in goodDict.iteritems():
            
            if(val > 2):
                output = str(pX) + " " + str(pY) + " " + str(key) + " " + str(val) + "\n"     
                f.write(output)
  
    f.close()
    
    

def doProbabilitySuccesPlot(settings, extraTitle):
    
        
    winPosDict = dict()
     
    goodDict = dict(myMultistrand.aFactory.result1.countDict)
    goodDictOther = dict(myMultistrand.aFactory.result2.countDict)

    print("Dict size is ", len(goodDict)) 

    for key, value in goodDict.iteritems():
        
        valueOther = 0
         
        if(goodDictOther.has_key(key)):
            valueOther = goodDictOther[key]
         
        
        denom = np.float(goodDict[key]) + np.float(valueOther)         
        newVal = np.float(value) / np.float(denom)
        winPosDict[key] = newVal

    extraTitle += "-winProb"
    
    doPosPlots(winPosDict, settings, extraTitle, settings.trials, "winprob")
    
    
def doBinaryProbabilityPlot(settings, extraTitle):
     
    def genPlots(settings, extraTitle, selectedCounts, result, extrastr): 
     
        goodDict = dict(result.binaryDict)
        plottingDict = dict()
        
        
        for key, value in goodDict.iteritems():
                        
            plottingDict[key] = np.float(value) / np.float(selectedCounts)
    
        extraTitle += "-binaryProb" + "-" + extrastr 
                
        doPosPlots(plottingDict, settings, extraTitle, selectedCounts, "binaryProb")

    myFact = myMultistrand.aFactory

    genPlots(settings, extraTitle, myFact.result0.pathCount.value, myFact.result0, "")
    genPlots(settings, extraTitle, myFact.result1.pathCount.value, myFact.result1, Options.STR_SUCCESS)
    genPlots(settings, extraTitle, myFact.result2.pathCount.value, myFact.result2, Options.STR_FAILURE)    



def simulationTimeBarplot(settings, extraTitle):
    
    fname = standardFileName(SCRIPT_DIR, settings.mySeq, extraTitle, settings.trials)
    
    all_times = np.array([i.time for i in myMultistrand.results])
    forward_times = np.array([i.time for i in myMultistrand.results if i.tag == Options.STR_SUCCESS])
    reverse_times = np.array([i.time for i in myMultistrand.results if i.tag == Options.STR_FAILURE  or i.tag == None])
    
   
    def makeFig(selector, times):
        
        fig = plt.figure()
        ax = fig.gca()

        ax.hist(times, 20, alpha=0.75)
        ax.set_title(("Trajectory times " + selector + " count=" + str(len(times))))      
         
        ax = plt.gca()
        ax.set_ylabel('Trajectory counts')    
        ax.set_xlabel('Trajectory time')
        
        plt.xticks(rotation=-40)
                
        plt.tight_layout()
        plt.savefig(fname + "-bar" + "-" + selector + '.pdf')
        plt.close()
        
        
        
    makeFig(Options.STR_ALL, all_times) 
    makeFig(Options.STR_SUCCESS, forward_times)
    makeFig(Options.STR_FAILURE, reverse_times)

    
    plt.figure()
    
    plt.style.use('seaborn-deep')

    plt.title("")
    plt.hist([forward_times, reverse_times], bins=12, alpha=0.7, label=['Success', 'Failure'])
    plt.legend(loc='upper right')   
    
    plt.xticks(rotation=-40)
            
    plt.tight_layout()
    plt.savefig(fname + "-bar" + "-" + "combined" + '.pdf')
    plt.close()
    
    
def makeAlignmentSuccessTable(analysis, settings, extraTitle):
    
    fname = standardFileName(SCRIPT_DIR, settings.mySeq, extraTitle, settings.trials)
    
    success = int(analysis.result1.pathCount.value)
    failed = int(analysis.result2.pathCount.value)
    
    def getAligned(result):
        return sum(props.aligned for props in result.pathProps)
         
    alignedSuccess = getAligned(analysis.result1)
    alignedFailed = getAligned(analysis.result2)
    
    nonalignedSuccess = success - alignedSuccess
    nonalignedFailed = failed - alignedFailed
    
    aligned = alignedSuccess + alignedFailed
    nonaligned = nonalignedSuccess + nonalignedFailed
    
    f = open(fname + "-AlignmentTable.txt", 'w')
    
    output = "Success-Aligned " + str(alignedSuccess) + "\n"
    output += "Success-NonAligned " + str(nonalignedSuccess) + "\n"
    output += "Failed-Aligned " + str(alignedFailed) + "\n"
    output += "Failed-NonAligned " + str(nonalignedFailed) + "\n"

    output += "\n"
    
    output += "Success " + str(success) + "\n"
    output += "Failed " + str(failed) + "\n"
    
    output += "\n"
    
    
    output += "aligned " + str(aligned) + "\n"
    output += "non-aligned " + str(nonaligned) + "\n"
    
    output += "\n"
    

    if aligned > 0 :
        output += "Aligned succes prob (%) =  " + str(100.0 * np.float(alignedSuccess) / np.float(aligned)) + "\n"
    if nonaligned > 0 :
        output += "Non aligned succes prob (%) =  " + str(100.0 * np.float(nonalignedSuccess) / np.float(nonaligned)) + "\n"


    f.write(output)
    
    



def doInference(mySeq, extraTitle, cutOff, runs):
    
    
    settings = hybridizationSimSettings(mySeq, cutOff, runs)
    
    first_step_simulation(myMultistrand, settings)

    # The actual number of simulations.
    settings.trials = myMultistrand.numOfThreads * myMultistrand.trialsPerThread
    
    # Uncomment this for the barplots.  #     simulationTimeBarplot(settings, extraTitle)
            
    doPosPlotsPrimer(myMultistrand.aFactory.result0, settings, extraTitle) 
    writeStructFile(myMultistrand.aFactory.result0, settings, extraTitle)
    
    doPosPlotsPrimer(myMultistrand.aFactory.result1, settings, extraTitle + "-SUCCESS") 
    writeStructFile(myMultistrand.aFactory.result1, settings, extraTitle + "-SUCCESS")

    doPosPlotsPrimer(myMultistrand.aFactory.result2, settings, extraTitle + "-FAILURE") 
    writeStructFile(myMultistrand.aFactory.result2, settings, extraTitle + "-FAILURE")
    
    # # Make the probability-success plot   
    doProbabilitySuccesPlot(settings, extraTitle)
    doBinaryProbabilityPlot(settings, extraTitle)

    # # Write the alignment - success table
    makeAlignmentSuccessTable(myMultistrand.aFactory, settings, extraTitle)
    myMultistrand.clear()
    



if __name__ == '__main__':
    
    if len(sys.argv) < 1:
        print """Usage:
              python hybridization_F3 <numOfThreads> <numOfPaths>  P0/P3/P4      \n
              Example: python hybridization_F3 2 100 P3
              """
        sys.exit()
        
    print sys.argv

    numOfThreads = np.int(sys.argv[1])
    numOfPaths = np.int(sys.argv[2])
    toggle = str(sys.argv[3])

    myMultistrand.setNumOfThreads(numOfThreads)
  
    
    if toggle == "test":
        doInference('TACCGT', "P0-test", 10, numOfPaths)  # P0
    if toggle == "test2":
        doInference(goa2006_P0, "P0-test", 10, numOfPaths)  # P0
        
    if toggle == "P0":
        doInference(goa2006_P0, toggle, 14, numOfPaths)  # P0
    if toggle == "P3":
        doInference(goa2006_P3, toggle, 14, numOfPaths)  # P3
    if toggle == "P4":
        doInference(goa2006_P4, toggle, 14, numOfPaths)  # P4



        


        
