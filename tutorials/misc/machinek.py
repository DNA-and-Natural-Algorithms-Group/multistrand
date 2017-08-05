# Frits Dannenberg, July 2017

# Simulation of Figure 2d in 
# Programmable energy landscapes for kinetic control of DNA strand displacement
# RRF Machinek, TE Ouldridge, NEC Haley, J Bath & AJ Turberfield. Nature communications

import sys, os

import xlrd
import matplotlib.pyplot as plt
import matplotlib.lines as lines

import numpy as np

from matplotlib.ticker import ScalarFormatter

from multistrand.concurrent import myMultistrand, FirstStepRate, Bootstrap
from multistrand.experiment import standardOptions, setBoltzmann
from multistrand.objects import StopCondition, Domain, Complex, Strand
from multistrand.options import Options
from msArrhenius import setArrheniusConstantsDNA23


myMultistrand.setNumOfThreads(8)

# Figure 2d has 3x12 = 36 rates plotted. 
# Input: 0 <= selector < 36  

# range 0 -11:  6 nt toehold
# range 12-23:  7 nt toehold
# range 24-36: 10 nt toehold

# order of mismatch position:
# perfect - 2 - 3 - 4 - 5 - 6 - 7 - 8 - 9 -10 - 12 - 14

positionSelector = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14]   


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
        initialComplex.boltzmann_supersample = 25
        initialInvader.boltzmann_supersample = 25

    stopSuccess = StopCondition(Options.STR_SUCCESS, [(successComplex, Options.dissocMacrostate, 0)])
    stopFailed = StopCondition(Options.STR_FAILURE, [(initialComplex, Options.dissocMacrostate, 0)])
    
    # actually set the intial and stopping states    
    options.start_state = [initialComplex, initialInvader]
    options.stop_conditions = [stopSuccess, stopFailed]         
    
        
def openDocument(document): 
    
    reader = xlrd.open_workbook(document, 'rb')
    row = reader.sheet_by_index(0)
    return row            

# trials has to be the first argument.
def genOptions(trialsIn, select):

    stdOptions = standardOptions(Options.firstStep, tempIn=23.0, trials=trialsIn, timeOut=0.1)
    stdOptions.sodium = expSodium(select)
    stdOptions.magnesium = expMagnesium(select)
    stdOptions.join_concentration = 5e-9
    
    # set the DNA23 parameters
    setArrheniusConstantsDNA23(stdOptions)
    
    machinek2014(stdOptions, select, trialsIn)
    return stdOptions


def computeRate(select, trials):
    
    myMultistrand.setOptionsFactory2(genOptions, trials, select) 
    myMultistrand.setTerminationCriteria(5)
    myMultistrand.run()
    
    myFSR = myMultistrand.results  # 5 nM concentration
    
    print myFSR
    
    
    return myFSR.k1()

def excelSelect(select):
    # Flip the row -- simply how the excell is formatted
    toeholdSelect = select / 12 
    mismatchSelect = select % 12
    
    #offset one for header
    return 1 + 12 * (2 - toeholdSelect) + mismatchSelect

def excelFind(row, col):
    dir = os.path.dirname(__file__)
    document = os.path.join(dir, 'data/machinek-figure2.xlsx')      
    myDoc = openDocument(document)
    
    return myDoc.cell(row, col).value

def expSodium(select):
    return excelFind(excelSelect(select),13)

def expMagnesium(select):
    return excelFind(excelSelect(select),14)

def measuredRate(select):
    return excelFind(excelSelect(select),9)


def generateGraph(trials=15):
    
    
    fig , ax= plt.subplots() 
    plt.title  ("Machinek Plot" , fontsize = 24) 
    plt.ylabel ( r"$\mathregular{K(M^{-1}s^{-1})}$",  fontsize = 19)

#            plt.legend( loc=4, borderaxespad=0., prop={'size':17})

    ax.xaxis.set_major_formatter(ScalarFormatter())
    plt.xticks([1,2,4,6,8,10,12,14,16])
    xticks = ax.xaxis.get_major_ticks()
    xticks[0].label1.set_visible(False)
    
    axes = plt.gca()
#    axes.set_ylim([10**2,10**8])
    axes.set_ylim([2,8])
    #plt.yticks(np.arange(100, 10**7 ,  2))
    plt.xlabel('Mismatch position (nt)',fontsize = 19)
    
    simRates = []
    realRates = []
    
    for select in range(0,12):
        simRates.append( np.log10(computeRate(select, trials)))
        realRates.append( np.log10(measuredRate(select)))
        
    print simRates
    print realRates
    
    ax.plot(positionSelector, realRates, linewidth=2)
    ax.plot(positionSelector, simRates, linewidth=2, linestyle = '--')
        
    #plt.show()
    
    filename = "plot" 
    if not os.path.exists(filename):
        os.makedirs(filename)
        
    plt.savefig(filename +'/machinek2014'+'.pdf', bbox_inches='tight')
    plt.close(fig) 
    
    
# The actual main method
if __name__ == '__main__':

    if len(sys.argv) > 1:
        print("No commandline arguments required. Aborting program.")
        exit()
    
    generateGraph(80)



























