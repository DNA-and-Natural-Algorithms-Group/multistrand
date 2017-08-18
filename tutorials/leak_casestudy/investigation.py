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

from multistrand.concurrent import myMultistrand, MergeSim, FirstStepRate, Bootstrap
from multistrand.objects import StopCondition, Domain, Strand, Complex
from multistrand.options import Options
from msArrhenius import setArrheniusConstantsDNA23
from LeakToolkit import getOptions, getRates
from SeesawGate import SeesawRates
import numpy as np


ATIME_OUT = 10.0
A_CONCENTRATION = 50e-9
TRIALS = 5000
DNA = "DNA"


myMultistrand.setNumOfThreads(8)



LONG_SEQ2 = "CCAAACAAAACCTAT"
LONG_SEQ5 = "AACCACCAAACTTAT"
LONG_SEQ6 = "CCTAACACAATCACT"
# some ive made up, but these shouldn't make much difference
LONG_SEQ7 = "CCACAAAACAAAACT"
LONG_SEQ1 = "CATCCATTCAACTAT"
LONG_SEQT = "CTCT"

# The actual main method
if __name__ == '__main__':
    toehold_dom = Domain(name="T", sequence=LONG_SEQT)
    base_dom = Domain(name="S2", sequence=LONG_SEQ2+LONG_SEQ2)
    input_dom = Domain(name="S5", sequence=LONG_SEQ5+LONG_SEQ5)
    output_dom = Domain(name="S6", sequence=LONG_SEQ6+LONG_SEQ6)

    input_strand = base_dom + toehold_dom + input_dom
    output_strand = output_dom + toehold_dom  + base_dom
    base_strand = toehold_dom.C + base_dom.C + toehold_dom.C

    gate_complex = Complex(strands=[base_strand, output_strand], structure = ".((+.))")
    input_complex = Complex(strands=[input_strand], structure="...")
    intermediate_complex = Complex(strands=[base_strand, output_strand, input_strand], structure = "(((+.).+.))")
    output_complex = Complex(strands=[output_strand], structure="...")

    success_stop_condition = StopCondition(
        Options.STR_SUCCESS, [(intermediate_complex, Options.looseMacrostate, 6)])
    failed_stop_condition = StopCondition(
        Options.STR_FAILURE, [(input_complex, Options.dissocMacrostate, 0)])

    start_time = time.time()
    myMultistrand.setOptionsFactory6(getOptions, TRIALS, DNA, gate_complex, input_complex, [
                                         success_stop_condition], [failed_stop_condition])
    print getRates()
    print "Took {} s".format(time.time() - start_time)

    success_stop_condition = StopCondition(
        Options.STR_SUCCESS, [(output_complex, Options.dissocMacrostate, 0)])
    start_time = time.time()
    myMultistrand.setOptionsFactory6(getOptions, TRIALS, DNA, gate_complex, input_complex, [
                                         success_stop_condition], [failed_stop_condition])
    print getRates()
    print "Took {} s".format(time.time() - start_time)
