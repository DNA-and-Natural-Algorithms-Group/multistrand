# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
Example of leak rates, AT ends vs CG ends.
Compatible with Multistrand 2.1 or higher
"""

# FD: This script is now set to use 4 threads and just 50,000 trajectories.
# FD: This is different from the results in case1.pdf
# FD: The results of this study heavily depend on the parameterization of the Metropolis model: JS or DNA23 (see below).

from multistrand.objects import StopCondition, Domain, Complex, Strand
from multistrand.options import Options, Literals
from multistrand.concurrent import MergeSim


ATIME_OUT = 10.0

myMultistrand = MergeSim()
myMultistrand.setNumOfThreads(8)
myMultistrand.setLeakMode()


def first_step_simulation(strand_seq, trials, T=25, material="DNA"):

    print(f"Running first step mode simulations for {strand_seq} (with Boltzmann sampling)...")

    # Using domain representation makes it easier to write secondary structures.
    onedomain = Domain(name="onedomain", sequence=strand_seq)
    gdomain = Domain(name="gdomain", sequence="TTTT")

    top = Strand(name="top", domains=[onedomain])
    bot = top.C
    dangle = Strand(name="Dangle", domains=[onedomain, gdomain])

    duplex_complex = Complex(strands=[top, bot], structure="(+)")
    invader_complex = Complex(strands=[dangle], structure="..")
    duplex_invaded = Complex(strands=[dangle, bot], structure="(.+)")

    # Declare the simulation complete if the strands become a perfect duplex.
    success_stop_condition = StopCondition(Literals.success, [(duplex_invaded, Literals.dissoc_macrostate, 0)])
    failed_stop_condition = StopCondition(Literals.failure, [(duplex_complex, Literals.dissoc_macrostate, 0)])

    for x in [duplex_complex, invader_complex]:
        x.boltzmann_count = trials
        x.boltzmann_sample = True

    # the first argument has to be trials.
    def getOptions(trials, duplex_complex, dangle, success_stop_condition, failed_stop_condition):

        o = Options(simulation_mode="First Step", substrate_type=material, rate_method="Metropolis",
                    num_simulations=trials, simulation_time=ATIME_OUT, temperature=T)
        o.DNA23Metropolis()
        o.start_state = [duplex_complex, dangle]
        o.stop_conditions = [success_stop_condition, failed_stop_condition]

        # FD: The result of this script depend significantly on JS or DNA23 parameterization.
#         o.JSMetropolis25()
        return o

    myMultistrand.setOptionsFactory5(getOptions, trials, duplex_complex,
                                     invader_complex, success_stop_condition, failed_stop_condition)
    myMultistrand.run()
    myFSR = myMultistrand.results
    

def doFirstStepMode(seq, T=25, material="DNA", numOfRuns=500):
    first_step_simulation(seq, numOfRuns, T=T, material=material)


def makePlots():
    seqs = list()
    seqs.append('GTCGATGC')   
    seqs.append('TCGAGTGA')
#    seqs.append('CCTACGTCTCACTAACG')
#     seqs.append('ACTACGTCTCACTAACG')

    for seq in seqs:
        doFirstStepMode(seq, numOfRuns=50000)


# The actual main method
if __name__ == '__main__':
    makePlots()
