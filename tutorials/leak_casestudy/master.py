import sys
import time
from os.path import expanduser
from datetime import datetime

# Tad dodgy, but using in non-root environments. Should work
dirs = ["~/workspace/multistrand", "~/workspace/multistrandPy",
        "~/multistrand", "~/multistrandPy"]
for x in dirs:
    i = expanduser(x)
    sys.path.append(i)

from LeakToolkit import setMinimumSuccess, calculateBaseOutputRate, calculateGateGateLeak, calculateBaseFuelRate, calculateGateFuelLeak, calculateBaseThresholdRate, calculateReverseFuelRate, calculateReverseOutputRate
from SeesawGate import NormalSeesawGate, MismatchedSeesawGate
from multistrand.experiment import ClampedSeesawGate
from HiddenSeesawGate import AntiLeakSeesawGate

USE_SHORT_DOMAINS = True

SHORT_SEQ1 = "ACCTCT"
SHORT_SEQ2 = "TCTTTA"
SHORT_SEQ7 = "ACATCC"
SHORT_SEQ5 = "TACTAC"
SHORT_SEQ6 = "ACCATT"
SHORT_SEQT = "CTCT"

LONG_SEQ2 = "CCAAACAAAACCTAT"
LONG_SEQ5 = "AACCACCAAACTTAT"
LONG_SEQ6 = "CCTAACACAATCACT"
# some ive made up, but these shouldn't make much difference
LONG_SEQ7 = "CCACAAAACAAAACT"
LONG_SEQ1 = "CATCCATTCAACTAT"
LONG_SEQT = SHORT_SEQT

CL_LONG_S18 = "TCTTCTAACAT"
CL_LONG_S5 = "CCACCAAACTT"
CL_LONG_S6 = "TAACACAATCA"
CL_LONG_S29 = "CCAATACTCCT"
CL_LONG_S53 = "TATCTAATCTC"
CL_LONG_S44 = "AAACTCTCTCT"
CL_LONG_SEQT = "TCT"
CLAMP_SEQ = "CA"

SHORT_GATE_A_SEQ = [SHORT_SEQ1, SHORT_SEQ2, SHORT_SEQ5, SHORT_SEQ7, SHORT_SEQT]
SHORT_GATE_B_SEQ = [SHORT_SEQ2, SHORT_SEQ5, SHORT_SEQ6, SHORT_SEQ7, SHORT_SEQT]
LONG_GATE_A_SEQ = [LONG_SEQ1, LONG_SEQ2, LONG_SEQ5, LONG_SEQ7, LONG_SEQT]
LONG_GATE_B_SEQ = [LONG_SEQ2, LONG_SEQ5, LONG_SEQ6, LONG_SEQ7, LONG_SEQT]
CL_LONG_GATE_A_SEQ = [CL_LONG_S44, CL_LONG_S18,
                      CL_LONG_S5, CL_LONG_S29, CL_LONG_SEQT, CLAMP_SEQ]
CL_LONG_GATE_B_SEQ = [CL_LONG_S53, CL_LONG_S5,
                      CL_LONG_S6, CL_LONG_S29, CL_LONG_SEQT, CLAMP_SEQ]

start_time = time.time()


def calcMetrics(gateA, gateB):
    rates = []
    print "\n **** Base Output Rates **** \n"
    rates.append(calculateBaseOutputRate(gateA))
    printTimeElapsed()
    rates.append(calculateBaseOutputRate(gateB))
    printTimeElapsed()

    print "\n **** Reverse Output Rates **** \n"
    rates.append(calculateReverseOutputRate(gateA))
    printTimeElapsed()
    rates.append(calculateReverseOutputRate(gateB))
    printTimeElapsed()

    print "\n **** Base Fuel Rates **** \n"
    rates.append(calculateBaseFuelRate(gateA))
    printTimeElapsed()
    rates.append(calculateBaseFuelRate(gateB))
    printTimeElapsed()

    print "\n **** Reverse Fuel Rates **** \n"
    rates.append(calculateReverseFuelRate(gateA))
    printTimeElapsed()
    rates.append(calculateReverseFuelRate(gateB))
    printTimeElapsed()

    print "\n **** Threshold Rates **** \n"
    # rates.append(calculateBaseThresholdRate(gateA))
    printTimeElapsed()
    # rates.append(calculateBaseThresholdRate(gateB))
    printTimeElapsed()

    return rates


def runMismatchSimulations(domainListA, domainListB):
    # 'recognition domain' sequence length
    recog_len = len(domainListA[1])
    # increment by twos
    mismatched_rates = []
    gateA = NormalSeesawGate(*domainListA)
    gateB = NormalSeesawGate(*domainListB)
    print " >>>>>> Normal Gate Rates <<<<<<<"
    mismatched_rates.append(calcMetrics(gateA, gateB))

    for i in range(2, recog_len, 2):
            # Not terribly efficient...
        rates = []
        gateA = MismatchedSeesawGate(*domainListA)
        gateB = MismatchedSeesawGate(*domainListB)

        # Assume that it is required for a mismatch in the output of the first
        # gate and the input of the second gate
        print " >>>>>> Mismatch in position %d <<<<<<<" % i

        # NB: If an upstream gate has a mismatch in its output, this implies
        # that the relevant downstream gate should have a mismatch in its output
        # The numbering remains the same since we are always 5' to 3' and the
        # orientation is not changing
        gateA.placeMismatchInOutput(i)
        gateB.placeMismatchInInputWire(i)

        rates = calcMetrics(gateA, gateB)

        mismatched_rates.append(rates)

    return mismatched_rates


def runClampedSimulations(domainListA, domainListB):
    # 'recognition domain' sequence length
    recog_len = len(domainListA[1])
    # increment by twos
    rates = []
    gateA = ClampedSeesawGate(*domainListA)
    gateB = ClampedSeesawGate(*domainListB)
    setMinimumSuccess(25)
    rates.append(calcMetrics(gateA, gateB))
    setMinimumSuccess(2)
    rates.append(calcLeakMetrics(gateA, gateB))
    return rates


def runAntiLeakSimulations(domainListA, domainListB):
    rates = []
    gateA = AntiLeakSeesawGate(*domainListA)
    gateB = AntiLeakSeesawGate(*domainListB)
    setMinimumSuccess(25)
    rates.append(calcMetrics(gateA, gateB))
    setMinimumSuccess(2)
    rates.append(calcLeakMetrics(gateA, gateB))
    return rates



def outputRates(rates, time_taken):
    # this method assumes that the rates object is a 2D array
    try:
        rates_measured = len(rates[0])
        param_measured = len(rates)
        file_name = "data{}.txt".format(str(datetime.now()))
        with open(file_name, "w") as output:
            output.write("Dataset Created On: {}\n Took {} s\n".format(
                str(datetime.now()), time_taken))
            for i in range(0, rates_measured):
                output_str = ""
                #output.write("writing rate " + str(i) + "\n")
                for j in range(0, param_measured):
                    # i.e. measure the SAME RATE accross different parameters
                    output_str += ("{},{},{},".format(
                        rates[j][i].k1, rates[j][i].k1_bounds[0], rates[j][i].k1_bounds[1]))
                    # If an alternative rate exists for that reaction do output that!
                    if rates[j][i].k1Alt != None:
                        output_str += ("[{},{},{},]".format(
                            rates[j][i].k1Alt, rates[j][i].k1_alt_bounds[0], rates[j][i].k1_alt_bounds[1]))
                output_str += "\n"
                output.write(output_str)
    except TypeError:
        print "Rates Input Object Invalid"
        print "No output has been printed"


def printTimeElapsed():
    curr_time = time.time()
    elap_time = curr_time - start_time
    if(elap_time < 1000):
        print "Time since start: {} s\n".format(elap_time)
    elif(elap_time < 6000):
        print "Time since start: {} min\n".format(elap_time / 60)
    else:
        print "Time since start: {} hr\n".format(elap_time / 3600)


def runAndLogClamped():
    start_time = time.time()
    data = runClampedSimulations(CL_LONG_GATE_A_SEQ, CL_LONG_GATE_B_SEQ)
    time_taken = time.time() - start_time
    outputRates(data, time_taken)


def runAndLogAntiLeak():
    CL_LONG_GATE_A_SEQ.extend(['GA', 'GT'])
    CL_LONG_GATE_B_SEQ.extend(['GA', 'GT'])
    start_time = time.time()
    data = runAntiLeakSimulations(CL_LONG_GATE_A_SEQ, CL_LONG_GATE_B_SEQ)
    time_taken = time.time() - start_time
    outputRates(data, time_taken)


if __name__ == '__main__':
    gateA = ClampedSeesawGate(*CL_LONG_GATE_A_SEQ)
    print gateA.gate_input_complex
    print gateA.gate_fuel_complex
    print gateA.gate_output_complex


def calcLeakMetrics(gateA, gateB):
    rates = []
    print "\n **** Fuel Leak Rate s ****\n"
    rates.append(calculateGateFuelLeak(gateA))
    printTimeElapsed()
    rates.append(calculateGateFuelLeak(gateB))
    printTimeElapsed()

    print "\n **** Gate Leak Rates **** \n"
    rates.append(calculateGateGateLeak(gateA, gateB))
    printTimeElapsed()
    print "Finished Rates \n"
    return rates
