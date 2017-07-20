import sys
from os.path import expanduser

home = expanduser("~/multistrand")
model = expanduser("~/multistrandPy")
sys.path.append(model)
sys.path.append(home)


from multistrand.objects import Domain
from leakToolkit import calculateBaseOutputRate, calculateGateGateLeak
from SeesawGate import NormalSeesawGate

default_num_trials = 1000


def setupNormalSimulations(s1_sequence, s2_sequence, s5_sequence,
                           s6_sequence, s7_sequence, toehold_sequence,
                           trials):
    # domains
    S1 = Domain(name="S1", sequence=s1_sequence)
    S2 = Domain(name="S2", sequence=s2_sequence)
    S5 = Domain(name="S5", sequence=s5_sequence)
    S6 = Domain(name="S6", sequence=s6_sequence)
    S7 = Domain(name="S7", sequence=s7_sequence)
    T = Domain(name="T", sequence=toehold_sequence)

    gateA = NormalSeesawGate(S1, S2, S5, S7, T)
    gateB = NormalSeesawGate(S2, S5, S6, S7, T)
    # to do - modify so it just takes the gate please
    calculateGateGateLeak(gateA, gateB, trials)
    calculateGateGateLeak(gateB, gateA, trials)


# The actual main method
if __name__ == '__main__':
    print("Initializing Example 3")
    if len(sys.argv) == 2:
        # MS: Integer number of trials, I hope ...
        num_trials = int(sys.argv[1])
        print("Trials: %d" % num_trials)
        # sequences from qian, winfree, 2011

        seq1 = "ACCTCT"
        seq2 = "TCTTTA"
        seq7 = "ACATCC"
        seq5 = "TACTAC"
        seq6 = "ACCATT"
        seqT = "CTCT"

        setupNormalSimulations(seq1, seq2, seq5, seq6, seq7, seqT, num_trials)
    else:
        print("(Default) Trials: %d" % default_num_trials)
