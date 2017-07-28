import sys
from os.path import expanduser
from datetime import datetime

# Tad dodgy, but using in non-root environments. Should work
dirs = ["~/workspace/multistrand", "~/workspace/multistrandPy",
        "~/multistrand", "~/multistrandPy"]
for x in dirs:
    i = expanduser(x)
    sys.path.append(i)

from LeakToolkit import calculateBaseOutputRate, calculateGateGateLeak, calculateBaseFuelRate, calculateGateFuelLeak
from SeesawGate import NormalSeesawGate, MismatchedSeesawGate

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
LONG_SEQ7 = "TCATTCCAACATTCA"
LONG_SEQ1 = "CCATTCAACTTAATC"
LONG_SEQT = SHORT_SEQT


SHORT_GATE_A_SEQ = [SHORT_SEQ1, SHORT_SEQ2, SHORT_SEQ5, SHORT_SEQ7, SHORT_SEQT]
SHORT_GATE_B_SEQ = [SHORT_SEQ2, SHORT_SEQ5, SHORT_SEQ6, SHORT_SEQ7, SHORT_SEQT]
LONG_GATE_A_SEQ = [LONG_SEQ1, LONG_SEQ2, LONG_SEQ5, LONG_SEQ7, LONG_SEQT]
LONG_GATE_B_SEQ = [LONG_SEQ2, LONG_SEQ5, LONG_SEQ6, LONG_SEQ7, LONG_SEQT]

def calcMetrics(gateA, gateB):
    rates = []
    print "\n **** Base Output Rates **** \n"
    rates.append(calculateBaseOutputRate(gateA))
    rates.append(calculateBaseOutputRate(gateB))

    print "\n **** Base Fuel Rates **** \n"
    rates.append(calculateBaseFuelRate(gateA))
    rates.append(calculateBaseFuelRate(gateB))

    print "\n **** Gate Leak Rates **** \n"
    rates.append(calculateGateGateLeak(gateA, gateB))

    print "\n **** Fuel Leak Rate s ****\n"
    rates.append(calculateGateFuelLeak(gateA))
    rates.append(calculateGateFuelLeak(gateB))
    print "Finished Rates \n"
    return rates

def setupNormalSimulations(trials, domainListA, domainListB):

    gateA = NormalSeesawGate(*domainListA)
    gateB = NormalSeesawGate(*domainListB)

    calculateGateGateLeak(gateA, gateB)


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

    outputRates(mismatched_rates)


def outputRates(rates):
    # this method assumes that the rates object is a 2D array
    try:
        rates_measured = len(rates[0])
        param_measured = len(rates)
        with open("data.txt", "w") as output:
            output.write("Dataset Created On: {}\n".format(
                str(datetime.now())))
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


# The actual main method
if __name__ == '__main__':
    print "Initializing Example 3"
    
    if USE_SHORT_DOMAINS:
        runMismatchSimulations(SHORT_GATE_A_SEQ, SHORT_GATE_B_SEQ)
    else:
        gateA = MismatchedSeesawGate(*LONG_GATE_A_SEQ)
        gateB = MismatchedSeesawGate(*LONG_GATE_A_SEQ)
