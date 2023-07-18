"""
This test compares Multistrand methods to compute dissociation rates for
duplexes.

We can either compute dissociation rate k- directly or compute k+ and use the
Boltzmann equation k+ / k- = exp( - dG / RT ) to infer k-. dG is obtained from
NUPACK. As it turns out, using either method is equivalent, both for Arrhenius
and Metropolis models.
"""

GAS_CONSTANT_R = 0.0019872036

import math
import numpy as np

# FIXME: referring to tutorials/compute/anneal.py
from anneal import compute as computeAnneal
# FIXME: referring to tutorials/compute/dissociation.py
from dissociation import compute as computeDissociation
from nupack import pfunc


def comparison(f):
    seq = "ATCCCAATCAACACCTTTCCTA"
    seqC = "TAGGAAAGGTGTTGATTGGGAT"

    temp = 25.0 + 273.15 + 40
    
    f.write(str(seq) + "   ")
    f.write(str("%0.3g" % (temp)) + "   ")
    
    # FD: just leaving this linking for now, meaning the simulation
    # settings will default to whatever is default for these functions
    
    predictedA = computeAnneal(seq, temp, 1.0)          
    predictedD = computeDissociation(seq, temp)
    
    # dotparen = "("*len(seq) + "+" + ")"*len(seq)
    dG = pfunc([seq, seqC], [1, 2], T=(temp - 273.15), material="dna")
    print(str(dG))
    
    kMinus = predictedA.k1() * math.exp(dG / (GAS_CONSTANT_R * temp))
    
    f.write(str("%0.3g" % np.log10(kMinus)) + "    ")
    f.write(str("%0.3g" % np.log10(predictedD.k1())) + "    ")
    
    diff = np.abs(np.log10(kMinus) - np.log10(predictedD.k1())) 
    
    f.write(str("%0.3g" % diff) + "    ")
    f.write("\n")
    f.flush()


if __name__ == '__main__':

    # For each column, compute the rate and print it to a file.
    resultFileName = "dissociation_comparison.txt"
    with open(resultFileName, 'w+') as f:
        f.write("Seq    temp     predicted-D    predicted-A   difference \n\n")
        comparison(f)
