import sys
from os.path import expanduser

home = expanduser("~/multistrand")
model = expanduser("~/multistrandPy")
sys.path.append(model)
sys.path.append(home)


from multistrand.objects import Domain
from LeakToolkit import calculateBaseOutputRate, calculateGateGateLeak
from SeesawGate import NormalSeesawGate, MismatchedSeesawGate

default_num_trials = 1000

seq1 = "ACCTCT"
seq2 = "TCTTTA"
seq7 = "ACATCC"
seq5 = "TACTAC"
seq6 = "ACCATT"
seqT = "CTCT"

S1 = Domain(name="S1", sequence=seq1)
S2 = Domain(name="S2", sequence=seq2)
S5 = Domain(name="S5", sequence=seq5)
S6 = Domain(name="S6", sequence=seq6)
S7 = Domain(name="S7", sequence=seq7)
T = Domain(name="T", sequence=seqT)


def setupNormalSimulations(trials, s1_sequence, s2_sequence, s5_sequence,
                           s6_sequence, s7_sequence, toehold_sequence):

    gateA = NormalSeesawGate(S1, S2, S5, S7, T)
    gateB = NormalSeesawGate(S2, S5, S6, S7, T)

    calculateGateGateLeak(gateA, gateB, trials)
    calculateGateGateLeak(gateB, gateA, trials)


def setupMismatchSimulations(trials, s1_sequence, s2_sequence, s5_sequence,
                             s6_sequence, s7_sequence, toehold_sequence):

    # Defaults to placing a C (unless the base their was a C before)
    # Will need to modify a second gate sequence slightly in order to
    # Ensure that the type of mismatch is always a C-C mismatch
    gateA = MismatchedSeesawGate(S1, S2, S5, S7, T,
                                 MismatchedSeesawGate.base_mismatch, 2)
    print(gateA.input_strand)


# The actual main method
if __name__ == '__main__':
    print("Initializing Example 3")
    if len(sys.argv) == 2:
        # MS: Integer number of trials, I hope ...
        num_trials = int(sys.argv[1])
        print("Trials: %d" % num_trials)
        # sequences from qian, winfree, 2011

        setupMismatchSimulations(num_trials, seq1, seq2, seq5, seq6, seq7,
                                 seqT)
    else:
        print("(Default) Trials: %d" % default_num_trials)
