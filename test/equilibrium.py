import sys, os, os.path
import pickle
import math

from multistrand.objects import Strand, Complex, Domain
from multistrand.options import Options
from multistrand.concurrent import MergeSim

import math
from nupack import *

# FD: This test attempts to compare the long-run partition function obtained from Multistrand and Nupack 


myMultistrand = MergeSim()

numOfPaths = 400.0
kBoltzmann = .00198717  # units of kcal/(mol*K)
RT = kBoltzmann * (273.15 + 37.0)

def setup_options(trials, seq, concentration):
    """
    setup_options( seq )

    creates an Options object using the sequence passed as a single
    domain with initially unpaired structure. 
    """
    d = Domain(name="initial", sequence=seq, length=len(seq))
    s = Strand(domains=[d])
    c = Complex(strands=[s], structure=".")
    
    o = Options(simulation_mode="Normal", parameter_type="Nupack", substrate_type="DNA",
                num_sims=trials, sim_time=0.008, start_state=[c])
    o.DNA23Metropolis()
    o.temperature = 310.15
    o.join_concentration = concentration
    return o


def nupack_pfunc(sequence, temperature):
    return pfunc([sequence], material="dna", T=temperature)

def run_distribution(seq):
    myMultistrand.setNumOfThreads(8)
    myMultistrand.setOptionsFactory4(setup_options, numOfPaths, seq, 2 * math.exp(-5.06 / RT), None)
    myMultistrand.run()    
    
    eq_dict = {}
    for end_state in myMultistrand.endStates:
        for cmplx in end_state:
            if cmplx[4] in eq_dict:
                count = eq_dict[cmplx[4]][1]
                eq_dict[cmplx[4]][1] = count + 1
            else:
                eq_dict[cmplx[4]] = [cmplx[5], 1]
    return eq_dict


NOW_TESTING = "GGGGAAACCCC"

eqd = run_distribution(NOW_TESTING)
print(eqd)
pf = nupack_pfunc(NOW_TESTING, 37.0)
print("NUPACK partition function is ", pf)

mySum = 0.0
for i in eqd.items():
    val = math.exp(-i[1][0] / RT)
    mySum += val
    
print("mySum is ", mySum)
print("myPart is ", math.log(mySum) * -RT)
