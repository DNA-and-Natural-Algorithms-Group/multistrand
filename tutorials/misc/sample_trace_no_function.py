# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from multistrand.objects import Complex, Domain, Strand, StopCondition
from multistrand.options import Options
from multistrand.system import SimSystem


def doSims(strandSeq, numTraj=1):
    o = Options(temperature=25, dangles="Some")

    o.num_simulations = numTraj
    o.output_interval = 1
    o.simulation_time = 10.0
    o.verbosity = 3

    # Using domain representation makes it easier to write secondary structures.
    onedomain = Domain(name="itall", sequence=strandSeq)
    top = Strand(name="top", domains=[onedomain])
    bot = top.C

    # Note that the structure is specified to be single stranded, but this will be over-ridden when Boltzmann sampling is turned on.
    startTop = Complex(strands=[top], structure=".")
    startBot = Complex(strands=[bot], structure=".")

    o.start_state = [startTop, startBot]

    # Stop when the exact full duplex is achieved.
    success_complex = Complex(strands=[top, bot], structure="(+)")
    stopSuccess = StopCondition("success", [(success_complex, 0, 0)])

    # Declare the simulation unproductive if the strands become single-stranded again.
    failed_complex = Complex(strands=[top], structure=".")
    stopFailed = StopCondition("failure", [(failed_complex, 2, 0)])

    o.stop_conditions = [stopSuccess, stopFailed]
    o.initial_seed = 1783

    s = SimSystem(o)
    s.start()

    return o, s


def main():
    return doSims("GCGTTTCAC")


if __name__ == '__main__':
    main()
