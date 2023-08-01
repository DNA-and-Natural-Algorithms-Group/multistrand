# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
Simulation of Figure 2d in
Programmable energy landscapes for kinetic control of DNA strand displacement
RRF Machinek, TE Ouldridge, NEC Haley, J Bath & AJ Turberfield. Nature communications
"""

import sys, os

import matplotlib
matplotlib.use('Agg')

import xlrd
import matplotlib.pyplot as plt

import numpy as np

from matplotlib.ticker import ScalarFormatter

from multistrand.concurrent import MergeSim
from multistrand.experiment import standardOptions, setBoltzmann
from multistrand.objects import StopCondition, Complex, Strand
from multistrand.options import Literals

myMultistrand = MergeSim()
myMultistrand.setNumOfThreads(2)

"""
 Figure 2d has 3x12 = 36 rates plotted. 
 Input: 0 <= selector < 36  

 range 0 -11:  6 nt toehold
 range 12-23:  7 nt toehold
 range 24-36: 10 nt toehold

 order of mismatch position:
 perfect - 2 - 3 - 4 - 5 - 6 - 7 - 8 - 9 -10 - 12 - 14
"""

positionSelector = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14]   

KINETICMETHOD = 1
NUMOFTHREADS = 2
NUMOFPATHS = 400
TERMINATIONCRIT = 2

DUMMYRUN = False


def machinek2014(options, selector, trialsIn):
    # we only allow first step mode at this point.
    
    # these are the sequences we need to build the dot-parens
    incumbent = ""
    target = ""
    invader = ""
    
    toeholdSelect = selector / 12 
    mismatchSelect = selector % 12
    
    positionSelector = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14]    
    mismatchSelect = positionSelector[mismatchSelect]
    
    # decide on toehold sequence
    toeholdSeq = "ATGTGG"  # 6 nt toehold option
    
    if toeholdSelect == 1 :
        toeholdSeq = "ATGTGGA"  # 7 nt toehold option
    if toeholdSelect == 2 :
        toeholdSeq = "ATGTGGAGGG"  # 10 nt toehold option

    # determine the incumbent, target and invader sequences
    # FD: copy-pasting supplementary Table 6 directly

    if mismatchSelect == 0 or mismatchSelect == 2 or mismatchSelect == 12 or mismatchSelect == 14  :
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTTTGAGGTTGA"
        target = "CCCTCCACATTCAACCTCAAACTCACC"
        
        if mismatchSelect == 0:  # perfect
            invader = "GGTGAGTTTGAGGTTGA"

        if mismatchSelect == 2:
            invader = "GGTGAGTTTGAGGTTCA"
        
        if mismatchSelect == 12:
            invader = "GGTGACTTTGAGGTTGA"
            
        if mismatchSelect == 14:
            invader = "GGTCAGTTTGAGGTTGA"
    
    if mismatchSelect == 3:
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTTTGAGGTGAT"
        target = "CCCTCCACATATCACCTCAAACTCACC"
        invader = "GGTGAGTTTGAGGTCAT"
        
    if mismatchSelect == 4:
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTTTGAGTGAGT"
        target = "CCCTCCACATACTCACTCAAACTCACC"
        invader = "GGTGAGTTTGAGTCAGT"

    if mismatchSelect == 5:
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTTTGATGAGGT"
        target = "CCCTCCACATACCTCATCAAACTCACC"
        invader = "GGTGAGTTTGATCAGGT"

    if mismatchSelect == 6:
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTTTGTGAAGGT"
        target = "CCCTCCACATACCTTCACAAACTCACC"
        invader = "GGTGAGTTTGTCAAGGT"
        
    if mismatchSelect == 7:
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTTTTGAGAGGT"
        target = "CCCTCCACATACCTCTCAAAACTCACC"
        invader = "GGTGAGTTTTCAGAGGT"

    if mismatchSelect == 8:
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTTTGATGAGGT"
        target = "CCCTCCACATACCTCATCAAACTCACC"
        invader = "GGTGAGTTTCATGAGGT"
        
    if mismatchSelect == 9:
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTTGATTGAGGT"
        target = "CCCTCCACATACCTCAATCAACTCACC"
        invader = "GGTGAGTTCATTGAGGT"
        
    if mismatchSelect == 10:
        incumbent = "TGGTGTTTGTGGGTGTGGTGAGTGATTTGAGGT"
        target = "CCCTCCACATACCTCAAATCACTCACC"
        invader = "GGTGAGTCATTTGAGGT"

    invader = invader + toeholdSeq
    
    # set up the actual complexes
    strandIncumbent = Strand(name="incumbent", sequence=incumbent)
    strandTarget = Strand(name="target", sequence=target)
    strandInvader = Strand(name="invader", sequence=invader)
    
    intialDotParen = '.' * 16 + '(' * 17 + "+" + '.' * 10 + ')' * 17  
    intialInvaderDotParen = '.' * len(invader)
    successDotParen = '.' * 33
    
    initialComplex = Complex(strands=[strandIncumbent, strandTarget], structure=intialDotParen)
    initialInvader = Complex(strands=[strandInvader], structure=intialInvaderDotParen)
    successComplex = Complex(strands=[strandIncumbent], structure=successDotParen)
    
    # Turns Boltzmann sampling on for this complex and also does sampling more efficiently by sampling 'trials' states.
    if(trialsIn > 0):
        setBoltzmann(initialComplex, trialsIn)
        setBoltzmann(initialInvader, trialsIn)

    stopSuccess = StopCondition(Literals.success, [(successComplex, Literals.dissoc_macrostate, 0)])
    stopFailed = StopCondition(Literals.failure, [(initialComplex, Literals.dissoc_macrostate, 0)])
    
    # actually set the intial and stopping states    
    options.start_state = [initialComplex, initialInvader]
    options.stop_conditions = [stopSuccess, stopFailed]         
    
        
def openDocument(document): 
    reader = xlrd.open_workbook(document, 'rb')
    row = reader.sheet_by_index(0)
    return row            


# trials has to be the first argument.
def genOptions(trialsIn, select):
    stdOptions = standardOptions(Literals.first_step, tempIn=23.0, trials=trialsIn, timeOut=100.0)
    stdOptions.sodium = expSodium(select)
    stdOptions.magnesium = expMagnesium(select)
    stdOptions.join_concentration = 5e-9
    
    # set the DNA23 parameters
    if KINETICMETHOD == 1:
        stdOptions.DNA23Metropolis()
    if KINETICMETHOD == 2:
        stdOptions.DNA23Arrhenius()   
    if KINETICMETHOD == 3:
        stdOptions.JSDefault()
    machinek2014(stdOptions, select, trialsIn)
    return stdOptions


def computeRate(select):
    myMultistrand.setOptionsFactory2(genOptions, NUMOFPATHS, select)
    myMultistrand.setNumOfThreads(NUMOFTHREADS)
    myMultistrand.setTerminationCriteria(TERMINATIONCRIT)
    myMultistrand.setLeakMode()  # the new leak object -- faster bootstrapping.
    myMultistrand.run()
    
    myFSR = myMultistrand.results  
    low, high = myFSR.doBootstrap()
    return myFSR.k1(), low, high


def excelSelect(select):
    # Flip the row -- simply how the excell is formatted
    toeholdSelect = select / 12 
    mismatchSelect = select % 12
    
    # offset one for header
    return 1 + 12 * (2 - toeholdSelect) + mismatchSelect


def excelFind(row, col):
    dir = os.path.dirname(__file__)
    document = os.path.join(dir, 'data/machinek-figure2.xlsx')      
    myDoc = openDocument(document)
    return myDoc.cell(row, col).value


def expSodium(select):
    return excelFind(excelSelect(select), 13)


def expMagnesium(select):
    return excelFind(excelSelect(select), 14)


def measuredRate(select):
    return excelFind(excelSelect(select), 9)


def generateGraph():
    fig , ax = plt.subplots()
    plt.title  ("Threeway strand displacement invasion mismatch \n Machinek et al." , fontsize=18) 
    plt.ylabel (r"$\mathregular{K(M^{-1}s^{-1})}$", fontsize=19)

    ax.xaxis.set_major_formatter(ScalarFormatter())
    plt.xticks([1, 2, 4, 6, 8, 10, 12, 14, 16])
    xticks = ax.xaxis.get_major_ticks()
    xticks[0].label1.set_visible(False)
    
    axes = plt.gca()
    axes.set_ylim([2, 8])
    plt.xlabel('Mismatch position (nt)', fontsize=19)

    colors = ['green', 'blue', 'black']
    
    for i in range(3):

        simRates = []
        simLow = []
        simHigh = []
        realRates = []

        for select in range(12 * i + 0, 12 * i + 12):
            if select % 12 < 11 and DUMMYRUN:
                rate = float(i+1) * (10 ** 5)
                low = float(i+1) * (10 ** 4)
                high = float(i+1) * (10 ** 6)
            else:
                rate, low, high = computeRate(select)
            
            simRates.append(np.log10(rate))
            simLow.append(np.log10(rate) - np.log10(low))
            simHigh.append(np.log10(high) - np.log10(rate))  
            
            realRates.append(np.log10(measuredRate(select)))
            
        print(simRates)
        print(realRates)
        
        ax.plot(positionSelector, realRates, linewidth=2, color=colors[i])
        ax.plot(positionSelector, simRates, linewidth=2, linestyle='--', color=colors[i])
        ax.errorbar(positionSelector, y=simRates, yerr=[simLow, simHigh], color=colors[i])
    # plt.show()
    
    filename = "plot" 
    if not os.path.exists(filename):
        os.makedirs(filename)
        
    plt.savefig(filename + '/machinek2014-' + str(KINETICMETHOD) + "-" + str(TERMINATIONCRIT) + "-" +str(DUMMYRUN) + '.pdf', bbox_inches='tight')
    plt.close(fig) 
    
    
# The actual main method
if __name__ == '__main__':

    if len(sys.argv) < 4:
        print("Please specify the kinetic method, number of successful trials, and number of threads and number of trials.")
        print("Example: python machinek.py metropolis 2 2 80")
        exit()

    if sys.argv[1] == "metropolis":
        KINETICMETHOD = 1
        
    if sys.argv[1] == "arrhenius":
        KINETICMETHOD = 2

    if sys.argv[1] == "jsdefault":
        KINETICMETHOD = 3

    TERMINATIONCRIT = int(sys.argv[2])
    NUMOFTHREADS = int(sys.argv[3])
    NUMOFPATHS = int(sys.argv[4])
    
    if len(sys.argv) > 5 and sys.argv[5] == "DUMMY":
        DUMMYRUN = True
        
    generateGraph()
