
import sys
from os.path import expanduser

# Tad dodgy, but using in non-root environments. Should work
dirs = ["~/workspace/multistrand", "~/workspace/multistrandPy", "~/multistrand", "~/multistrandPy"]
for x in dirs:
    i = expanduser(x)
    sys.path.append(i)

from multistrand.objects import Domain
from LeakToolkit import calculateBaseOutputRate, calculateGateGateLeak, calculateBaseFuelRate, calculateGateFuelLeak
from SeesawGate import NormalSeesawGate, MismatchedSeesawGate

DEFAULT_NUM_TRIALS = 200000

SEQ1 = "ACCTCT"
SEQ2 = "TCTTTA"
SEQ7 = "ACATCC"
SEQ5 = "TACTAC"
SEQ6 = "ACCATT"
SEQT = "CTCT"

S1 = Domain(name="S1", sequence=SEQ1)
S2 = Domain(name="S2", sequence=SEQ2)
S5 = Domain(name="S5", sequence=SEQ5)
S6 = Domain(name="S6", sequence=SEQ6)
S7 = Domain(name="S7", sequence=SEQ7)
T = Domain(name="T", sequence=SEQT)


def setupNormalSimulations(trials, domainListA, domainListB):

    gateA = NormalSeesawGate(*domainListA)
    gateB = NormalSeesawGate(*domainListB)

    # The base rates are far, far more likely that the leak rate
    # - why simulate them for so long?
    print trials
    print("Base Output")
    calculateBaseOutputRate(gateA, trials / 100)
    calculateBaseOutputRate(gateB, trials / 100)

    print("Base Fuel")
    calculateBaseFuelRate(gateA, trials / 100)
    calculateBaseFuelRate(gateB, trials / 100)

    print("Gate-Gate Leak")
    calculateGateGateLeak(gateA, gateB, trials)
    calculateGateGateLeak(gateB, gateA, trials)

    print("Gate-Fuel Leak")
    calculateGateFuelLeak(gateA, trials)
    calculateGateFuelLeak(gateB, trials)


def setupMismatchSimulations(trials, domainListA=None, domainListB=None):
    # Defaults to placing a C (unless the base their was a C before)
    # Will need to modify a second gate sequence slightly in order to
    # Ensure that the type of mismatch is always a C-C mismatch
    # this is a small change
    gateA = MismatchedSeesawGate(*domainListA)
    gateB = MismatchedSeesawGate(*domainListB)

    # Assume that it is required for a mismatch in the output of the first
    # gate and the input of the second gate

    for i in [4]:
        print "Mismatch in position %d" % i
        gateA.placeOutputMismatchInOutput(
            i, domainListA[2], domainListA[1], domainListA[4])
        gateB.placeInputMismatchInBase(
            i, domainListA[0], domainListA[1], domainListA[4])
        # NB: If an upstream gate has a mismatch in its output, this implies
        # that the relevant downstream gate should have a mismatch in its output
        # The numbering remains the same since we are always 5' to 3' and the
        # orientation is not changing
        print "Base Output"
        calculateBaseOutputRate(gateA, trials / 100)
        calculateBaseOutputRate(gateB, trials / 100)

        print("Base Fuel")
        calculateBaseFuelRate(gateA, trials / 100)
        calculateBaseFuelRate(gateB, trials / 100)

        print("Gate-Gate Leak")
        calculateGateGateLeak(gateA, gateB, trials)
        calculateGateGateLeak(gateB, gateA, trials)

        print("Gate-Fuel Leak")
        calculateGateFuelLeak(gateA, trials)
        calculateGateFuelLeak(gateB, trials)



# The actual main method
if __name__ == '__main__':
    print "Initializing Example 3"
    if len(sys.argv) == 2:
        # MS: Integer number of trials, I hope ...
        num_trials = int(sys.argv[1])
        print "Trials: %d" % num_trials
        # sequences from qian, winfree, 2011

    else:
        print "(Default) Trials: %d" % DEFAULT_NUM_TRIALS
        domainList_A = [S1, S2, S5, S7, T]
        domainList_B = [S2, S5, S6, S7, T]
        setupNormalSimulations(DEFAULT_NUM_TRIALS, domainList_A, domainList_B)
        setupMismatchSimulations(DEFAULT_NUM_TRIALS, domainList_A,  domainList_B)
