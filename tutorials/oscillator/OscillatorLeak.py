# Mrinank Sharma, July 2017
# Frits Dannenberg, Aug 2017
# Simulates Leak Reactions from a DSD oscillator and calculates the rate with bootstraping
import sys, time, os
from os.path import expanduser

import xlrd

import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.ticker import ScalarFormatter

from multistrand.concurrent import myMultistrand, FirstStepRate, Bootstrap
from multistrand.experiment import standardOptions, setBoltzmann
from multistrand.objects import StopCondition, Domain, Complex, Strand
from multistrand.options import Options
from msArrhenius import setArrheniusConstantsDNA23

import numpy as np


REL_DATA_DIR = 'data/table-s5.xlsx'

NORMAL = 0
WITHOUT_GG = 1
WITHOUT_G = 2
HELPER_WITHOUT_CC = 3

# to make - the excel spreadsheet

# To correct - see excel spreadsheet
SODIUM_COL = 11
MAGNESIUM_COL = 12
TEMP_COL = 5
MEASURED_RATE_COL = 7

myMultistrand.setNumOfThreads(4)

# i.e. setup options obiect, given the type of experiment and
# number of trials

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
        len(ap.sequence) - toehold_length - h_length + 2) + "+" + '.' * (toehold_length + h_length) + ')' * (len(ap.sequence) - toehold_length - h_length)

    if expirement_type == WITHOUT_GG:
        
        helper = Strand(name="Helper_AAq", sequence="TTTCCTAATCCCAATCAACACCTTTCCTA")
        produce_bot = Strand(name="Produce_BOT_CApAq", sequence="GTAAAGACCAGTGGTGTGAAGATAGGAAAGGTGTTGATTGGGATTAGGAAACC")
        ap = Strand(name="ap", sequence="CATCACTATCAATCATACATTTTCCTATCTTCACACCACTGG")
        aq = Strand(name="aq", sequence="CATCACTATCAATCATACATGGTTTCCTAATCCCAATCAACACC")

        # bot, aq, ap
        # Offsets due to a) clamp domains b) two b.p. removal
        produce_struct = "." * toehold_length + "(" * (len(produce_bot.sequence) - toehold_length) + "+" + '.' * (toehold_length + h_length - 2) + ')' * (
            len(ap.sequence) - toehold_length - h_length + 2) + "+" + '.' * (toehold_length + h_length - 2) + ')' * (len(ap.sequence) - toehold_length - h_length + 2)

#     elif expirement_type == WITHOUT_G:
#         helper = Strand(name="Helper_AAq",
#                         sequence="TTTCCTAATCCCAATCAACACCTTTCCTA")
#         produce_bot = Strand(name="Produce_BOT_CApAq",
#                              sequence="GTAAAGACCAGTGGTGTGAAGATAGGAAAGGTGTTGATTGGGATTAGGAAACC")
#         ap = Strand(
#             name="ap", sequence="CATCACTATCAATCATACATGTTTCCTATCTTCACACCACTGG")
#         aq = Strand(
#             name="aq", sequence="CATCACTATCAATCATACATGGTTTCCTAATCCCAATCAACACC")
# 
#         # bot, aq, ap
#         # Offsets due to a) clamp domains b) one b.p. removal
#         produce_struct = "." * toehold_length + "(" * (len(produce_bot.sequence) - toehold_length) + "+" + '.' * (toehold_length + h_length - 2) + ')' * (
#             len(ap.sequence) - toehold_length - h_length + 2) + "+" + '.' * (toehold_length + h_length - 1) + ')' * (len(ap.sequence) - toehold_length - h_length + 1)
#     elif expirement_type == HELPER_WITHOUT_CC:
#         # only modify helper sequence here - remove the two 3' most 'C'.
#         helper = Strand(name="Helper_AAq",
#                         sequence="TTTCCTAATCCCAATCAACACCTTTTA")
# 
#         # No change required elsewhere -  we check for the release of strands rather than the complicated
#         # leak complex formed for simplicity. We should really check for ANY free strands here i.e. Ap OR Aq
#         # but it is hard to imagine a mechanism which results in the release of Ap in this simulation
#         produce_bot = Strand(name="Produce_BOT_CApAq",
#                              sequence="GTAAAGACCAGTGGTGTGAAGATAGGAAAGGTGTTGATTGGGATTAGGAAACC")
#         ap = Strand(
#             name="ap", sequence="CATCACTATCAATCATACATGGTTTCCTATCTTCACACCACTGG")
#         aq = Strand(
#             name="aq", sequence="CATCACTATCAATCATACATGGTTTCCTAATCCCAATCAACACC")
# 
#         # Offset of two to account for clamp domains
#         produce_struct = "." * toehold_length + "(" * (len(produce_bot.sequence) - toehold_length) + "+" + '.' * (toehold_length + h_length - 2) + ')' * (
#             len(ap.sequence) - toehold_length - h_length + 2) + "+" + '.' * (toehold_length + h_length) + ')' * (len(ap.sequence) - toehold_length - h_length)

    produce_complex = Complex(strands=[produce_bot, aq, ap], structure=produce_struct)
    produce_struct = Complex(strands=[ap], structure='.' * len(ap.sequence))
    helper_complex = Complex(strands=[helper], structure='.' * len(helper.sequence))
    leak_complex = Complex(strands=[aq], structure='.' * len(aq.sequence))

    if trials > 0:
        setBoltzmann(produce_complex, trials)
        setBoltzmann(helper_complex, trials)

    success_stop_cond = StopCondition(
        Options.STR_SUCCESS, [leak_complex, Options.dissocMacrostate, 0])
    # the leak has failed if we end up with our initial complexes again.
    # check if we end up with a free helper complex
    failure_stop_cond = StopCondition(
        Options.STR_FAILURE, [helper_complex, Options.dissocMacrostate, 0])

    options.start_state = [produce_complex, helper_complex]
    options.stop_conditions = [success_stop_cond, failure_stop_cond]


def openDocument(document):
    reader = xlrd.open_workbook(document)
    sheet = reader.sheet_by_index(0)
    return sheet


def genOptions(trialsIn, experiment_type=NORMAL):
    # NB: Time out MUST be a float
    stdOptions = standardOptions(
        Options.firstStep, expTemp(experiment_type), trials=trialsIn, timeOut=0.00000001)
    stdOptions.temperature = expTemp(experiment_type)
    stdOptions.sodium = expSodium(experiment_type)
    stdOptions.magnesium = expMagnesium(experiment_type)
    #stdOptions.DNA23Metropolis()
    setArrheniusConstantsDNA23(stdOptions)
    changeComplex(stdOptions, experiment_type)

    return stdOptions


def computeRate(trialsIn, experiment_type=NORMAL):
    myMultistrand.setOptionsFactory2(genOptions, trialsIn, experiment_type)

    # use the new leak rates class for memory efficiency
    myMultistrand.setLeakMode()
#     myMultistrand.printStates()
#     myMultistrand.initialInfo()

    
    myMultistrand.run()
 
    results = myMultistrand.results
    print results
    # see above - no alternative success conditions defined
    confidence = Bootstrap(results, computek1=True)
    print confidence
# 
    return results.k1(), confidence


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


def generateGraph(trials=50):
    
    plt.rcdefaults()
    width = 0.35
    measuredRates = []
    simRates = []
    low_error = []
    high_error = []
    
#     for x in [NORMAL, WITHOUT_GG, WITHOUT_G, HELPER_WITHOUT_CC]:
    for x in [NORMAL ]:
        results = computeRate(trials, x)
        k1 = results[0].k1()
        lower_bound = results[1].ninetyFivePercentiles()[0]
        upper_bound = results[1].ninetyFivePercentiles()[1]
        expRate = measuredRate(x)
        measuredRates.append(k1)
        low_error.append(k1 - lower_bound)
        high_error.append(upper_bound - k1)
        measuredRates.append(expRate)

    fig, ax = plt.subplots()
    plt.ylabel(r"$\mathregular{Leak\ Rate\ (M^{-1}s^{-1})}$",  fontsize=12)
    experiment_types = ('Normal', 'W/O GG', 'W/O G', 'HELPER W/O CC')
    n = 4
    X = np.arange(n)
    plt.xticks(X + width / 2, experiment_types)
    plt.bar(X, measuredRates, width, facecolor="#9999ff",
            yerr=[low_error, high_error])
    plt.bar(X + width, simRates, width, facecolor="#ff9999")
    plt.title("Leak Rates in Oscillator Circuit")
    plt.gca().set_ylim(bottom=0)
    plt.show()



# the main method
if __name__ == '__main__':
    generateGraph(5)
 