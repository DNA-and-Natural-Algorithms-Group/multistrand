# Mrinank Sharma, July 2017
# Frits Dannenberg, Aug 2017
# Simulates Leak Reactions from a DSD oscillator and calculates the rate with bootstraping

import sys, time, os
from os.path import expanduser

dirs = ["~/workspace/multistrand", "~/workspace/multistrandPy",
        "~/multistrand", "~/multistrandPy"]
for x in dirs:
    i = expanduser(x)
    sys.path.append(i)

import xlrd

import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.ticker import ScalarFormatter

from multistrand.concurrent import MeregSim, FirstStepRate, Bootstrap
from multistrand.experiment import standardOptions, setBoltzmann
from multistrand.objects import StopCondition, Domain, Complex, Strand
from multistrand.options import Options, Literals

import numpy as np



REL_DATA_DIR = 'data/table-s5.xlsx'

NORMAL = 0
WITHOUT_GG = 1
WITHOUT_G = 2
HELPER_WITHOUT_CC = 3

SODIUM_COL = 11
MAGNESIUM_COL = 12
TEMP_COL = 5
MEASURED_RATE_COL = 7

myMultistrand = MergeSim()
myMultistrand.setNumOfThreads(8)

def changeComplex(options, expirement_type=NORMAL, trials=500):
    # Easiest to use full dot paren here
    # See SI Document for these sequences and lengths!
    toehold_length = 7
    h_length = 15

    helper = Strand(name="Helper_AAq", sequence="TTTCCTAATCCCAATCAACACCTTTCCTA")
    produce_bot = Strand(name="Produce_BOT_CApAq", sequence="GTAAAGACCAGTGGTGTGAAGATAGGAAAGGTGTTGATTGGGATTAGGAAACC")
    ap = Strand(name="ap", sequence="CATCACTATCAATCATACATGGTTTCCTATCTTCACACCACTGG")
    aq = Strand(name="aq", sequence="CATCACTATCAATCATACATGGTTTCCTAATCCCAATCAACACC")

    # Offset of two to account for clamp domains
    produce_struct = "." * toehold_length + "(" * (len(produce_bot.sequence) - toehold_length) + "+" + '.' * (toehold_length + h_length - 2) + ')' * (
        len(ap.sequence) - toehold_length - h_length + 2 ) + "+" + '.' * (toehold_length + h_length) + ')' * (len(ap.sequence) - toehold_length - h_length)

    if expirement_type == WITHOUT_GG:
        ap = Strand(name="ap", sequence="CATCACTATCAATCATACATTTTCCTATCTTCACACCACTGG")
 
        # bot, aq, ap
        # Offsets due to a) clamp domains b) two b.p. removal
        produce_struct = "." * toehold_length + "(" * (len(produce_bot.sequence) - toehold_length) + "+" + '.' * (toehold_length + h_length - 2) + ')' * (
            len(ap.sequence) - toehold_length - h_length + 4) + "+" + '.' * (toehold_length + h_length - 2) + ')' * (len(ap.sequence) - toehold_length - h_length + 2)

    elif expirement_type == WITHOUT_G:
        ap = Strand(
            name="ap", sequence="CATCACTATCAATCATACATGTTTCCTATCTTCACACCACTGG")
 
        # bot, aq, ap
        # Offsets due to a) clamp domains b) one b.p. removal
        produce_struct = "." * toehold_length + "(" * (len(produce_bot.sequence) - toehold_length) + "+" + '.' * (toehold_length + h_length - 2) + ')' * (
            len(ap.sequence) - toehold_length - h_length + 3) + "+" + '.' * (toehold_length + h_length - 1) + ')' * (len(ap.sequence) - toehold_length - h_length + 1)
            
    elif expirement_type == HELPER_WITHOUT_CC:
        # only modify helper sequence here - remove the two 3' most 'C'.
        helper = Strand(name="Helper_AAq",
                        sequence="TTTCCTAATCCCAATCAACACCTTTTA")
 
        # No change required elsewhere -  we check for the release of strands rather than the complicated
        # leak complex formed for simplicity. We should really check for ANY free strands here i.e. Ap OR Aq
        # but it is hard to imagine a mechanism which results in the release of Ap in this simulation
 
    produce_complex = Complex(name="produce", strands=[produce_bot, aq, ap], structure=produce_struct)
    helper_complex = Complex(name="helper",strands=[helper], structure='.' * len(helper.sequence))
    leak_complex = Complex(name="leak",strands=[aq], structure='.' * len(aq.sequence))
    
    if trials > 0:
        setBoltzmann(produce_complex, trials, 75)
        setBoltzmann(helper_complex, trials, 75)

    success_stop_cond = StopCondition(Literals.success, [(leak_complex, Options.dissoc_macrostate, 0)])
    # the leak has failed if we end up with our initial complexes again.
    # check if we end up with a free helper complex
    failure_stop_cond = StopCondition(Literals.failure, [(helper_complex, Options.dissoc_macrostate, 0)])

    options.start_state = [produce_complex, helper_complex]
    options.stop_conditions = [success_stop_cond, failure_stop_cond]


def openDocument(document):
    reader = xlrd.open_workbook(document)
    sheet = reader.sheet_by_index(0)
    return sheet


def genOptions(trialsIn, experiment_type=NORMAL):
    # NB: Time out MUST be a float
    stdOptions = standardOptions(Options.firstStep, expTemp(experiment_type), trials=trialsIn, timeOut=0.1)
    changeComplex(stdOptions, experiment_type)
    stdOptions.temperature = expTemp(experiment_type)
    stdOptions.sodium = expSodium(experiment_type)
    stdOptions.magnesium = expMagnesium(experiment_type)
    stdOptions.gt_enable = 1
    
    stdOptions.DNA23Metropolis()
    
    

    return stdOptions


def computeRate(trialsIn, experiment_type=NORMAL):
    
    myMultistrand.setOptionsFactory2(genOptions, trialsIn, experiment_type)

    # use the new leak rates class for memory efficiency
    myMultistrand.setLeakMode()
    
    myMultistrand.run()
 
    results = myMultistrand.results
    print(results)
    # see above - no alternative success conditions defined
    confidence = Bootstrap(results, computek1=True)
    print(confidence)
# 
    return results, confidence


def excelSelect(experiment_type):
    # excel sheet is formatted correctly, just need to plus 1 here.
    return experiment_type + 1


def excelFind(row, col):
    # i.e. the path of this file
    dir = os.path.dirname(__file__)
    document = os.path.join(dir, REL_DATA_DIR)
    #document = os.path.join(dir, 'data/table-s5.xlsx')
    doc = openDocument(document)

    return doc.cell(row, col).value


def expSodium(experiment_type):
    return excelFind(excelSelect(experiment_type), SODIUM_COL)


def expMagnesium(experiment_type):
    return excelFind(excelSelect(experiment_type), MAGNESIUM_COL)


def expTemp(experiment_type):
    return excelFind(excelSelect(experiment_type), TEMP_COL)


def measuredRate(experiment_type):
    return excelFind(excelSelect(experiment_type), MEASURED_RATE_COL)


def generateGraph(min_success, increment_trials):
    myMultistrand.setTerminationCriteria(min_success)
    plt.rcdefaults()
    width = 0.35
    measuredRates = [82, 11, 3, 28]
    simRates = []
    low_error = []
    high_error = []
    
#     for x in [WITHOUT_G]:
    for x in [NORMAL, WITHOUT_GG, WITHOUT_G, HELPER_WITHOUT_CC]:
        results = computeRate(increment_trials, x)
        k1 = results[0].k1()
        lower_bound = results[1].ninetyFivePercentiles()[0]
        upper_bound = results[1].ninetyFivePercentiles()[1]
        # expRate = measuredRate[x]
        # measuredRates.append(k1)
        simRates.append(k1)
        low_error.append(k1 - lower_bound)
        high_error.append(upper_bound - k1)
        # measuredRates.append(expRate)

    fig, ax = plt.subplots()
    plt.ylabel(r"$\mathregular{Leak\ Rate\ (M^{-1}s^{-1})}$",  fontsize=12)
    experiment_types = ('Normal', 'W/O GG', 'W/O G', 'HELPER W/O CC')
    n = 4
    X = np.arange(n)
    plt.xticks(X + width / 2, experiment_types)
    plt.bar(X, simRates, width, facecolor="#9999ff",
            yerr=[low_error, high_error], label="Simulated")
    plt.bar(X + width, measuredRates, width, facecolor="#ff9999", label="Measured")
    plt.legend()
    plt.title("Leak Rates in Oscillator Circuit")
    plt.gca().set_ylim(bottom=0)
    plt.show()



# the main method
if __name__ == '__main__':
    generateGraph(6, 10000)
 