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

from LeakToolkit import *
from SeesawGate import NormalSeesawGate, MismatchedSeesawGate, ClampedSeesawGate
from HiddenSeesawGate import *

# Non-clamped sequences
LONG_SEQ2 = "CCAAACAAAACCTAT"
LONG_SEQ5 = "AACCACCAAACTTAT"
LONG_SEQ6 = "CCTAACACAATCACT"
# some ive made up, but these shouldn't make much difference
LONG_SEQ7 = "CCACAAAACAAAACT"
LONG_SEQ1 = "CATCCATTCAACTAT"
SEQT = "CTCT"

# Clamped Sequences
CL_LONG_S18 = "TCTTCTAACAT"
CL_LONG_S5 = "CCACCAAACTT"
CL_LONG_S6 = "TAACACAATCA"
CL_LONG_S29 = "CCAATACTCCT"
CL_LONG_S53 = "TATCTAATCTC"
CL_LONG_S44 = "AAACTCTCTCT"
CL_SEQT = "TCT"
CLAMP_SEQ = "CA"

# Longer Toehold
CL_LONG_SEQT = "TCTCT"

def setupSims():
    setMaxTrials(25000000)
    setMinimumSuccess(50)
    # perhaps too high....
    setMinimumSuccess(30)

def runAndLogNormal():
    setupSims()
    domainListA = [LONG_SEQ1, LONG_SEQ2, LONG_SEQ5, LONG_SEQ7, SEQT]
    domainListB = [LONG_SEQ2, LONG_SEQ5, LONG_SEQ6, LONG_SEQ7, SEQT]
    gateA = NormalSeesawGate(*domainListA)
    gateB = NormalSeesawGate(*domainListB)
    calcGateMetrics(gateA)
    calculateGateGateLeak(gateA, gateB)


def runAndLogClamped():
    setupSims()
    domainListA = [CL_LONG_S44, CL_LONG_S18,
                  CL_LONG_S5, CL_LONG_S29, CL_SEQT, CLAMP_SEQ]
    domainListB = [CL_LONG_S53, CL_LONG_S5,
                      CL_LONG_S6, CL_LONG_S29, CL_SEQT, CLAMP_SEQ]
    gateA = ClampedSeesawGate(*domainListA)
    gateB = ClampedSeesawGate(*domainListB)
    calcGateMetrics(gateA)
    calculateGateGateLeak(gateA, gateB)

def runAndLogDoubleMismatchAntiLeak():
    setupSims()
    domainListA = [CL_LONG_S44, CL_LONG_S18,
                  CL_LONG_S5, CL_LONG_S29, CL_SEQT, CLAMP_SEQ]
    domainListB = [CL_LONG_S53, CL_LONG_S5,
                      CL_LONG_S6, CL_LONG_S29, CL_SEQT, CLAMP_SEQ]
    gateA = DoubleCentralMismatchSeesawGate(*domainListA)
    gateB = DoubleCentralMismatchSeesawGate(*domainListB)
    calcGateMetrics(gateA)
    calculateGateGateLeak(gateA, gateB)

def runAndLogSingleMismatchAntiLeak():
    setupSims()
    # note the use of long sequences here...
    domainListA = [CL_LONG_S44, CL_LONG_S18,
                  CL_LONG_S5, CL_LONG_S29, CL_LONG_SEQT, CLAMP_SEQ]
    domainListB = [CL_LONG_S53, CL_LONG_S5,
                      CL_LONG_S6, CL_LONG_S29, CL_LONG_SEQT, CLAMP_SEQ]
    gateA = CentralMismatchSeesawGate(*domainListA)
    gateB = CentralMismatchSeesawGate(*domainListB)
    calcGateMetrics(gateA)
    calculateGateGateLeak(gateA, gateB)

def calcGateMetrics(gateA):
    calculateBaseOutputRate(gateA)
    calculateReverseOutputRate(gateA)
    calculateBaseFuelRate(gateA)
    calculateReverseFuelRate(gateA)
    calculateBaseThresholdRate(gateA)
    calculateInputGateOcclusion(gateA)
    calculateInputGateOcclusionUnbind(gateA)
    calculateOutputThresholdOcclusion(gateA)
    calculateOutputThresholdOcclusionUnbind(gateA)
    calculateOutputGateOcclusion(gateA)
    calculateOutputGateOcclusionUnbind(gateA)
    calculateGateFuelLeak(gateA)
