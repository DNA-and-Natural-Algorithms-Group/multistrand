## Comparison of Multistrand methods to compute dissociation rates for duplexes.

## We can either compute dissociation rate k- directly or compute k+ and
## use the Boltzmann equation k+ / k- = exp(  - dG / RT  ) to infer k-.
## dG is obtained from NUPACK. 

## As it turns out, using either method is equivalent, both for Arrhenius and Metropolis models. 

GAS_CONSTANT_R = 0.0019872036

import sys, os
import math
import numpy as np

from anneal import compute as computeAnneal
from dissociation import compute as computeDissociation
from nupack import pfunc


## For each column, compute the rate and print it to a file. 
resultFileName = "dissociation_comparison.txt"
file = open(resultFileName, 'w+')

file.write("Seq    temp     predicted-D    predicted-A   difference \n\n") 


def comparison():

    for seq, seqC in zip( ["ACTGAAGATGACGA"], ["TCGTCATCTTCAGT"] ):
        
        temp = 25.0 + 273.15  + 40
        
        file.write(str(seq) + "   ")
        file.write(str( "%0.3g" %  (temp) ) + "   ")
        
        # FD: just leaving this linking for now, meaning the simulation
        # settings will default to whatever is default for these functions
        
        predictedA = computeAnneal(seq, temp, 1.0)          
        predictedD = computeDissociation(seq, temp)
        
        dotparen = "("*len(seq) + "+" + ")"*len(seq)
        
        dG = pfunc([seq, seqC], [1, 2], T=(temp-273.15), material="dna")
        print (str(dG)) 
        
        kMinus = predictedA.k1() * math.exp( dG / ( GAS_CONSTANT_R * temp) ) 
        
        file.write(str( "%0.3g" % np.log10(kMinus)) +     "    "  )
        file.write(str( "%0.3g" % np.log10(predictedD.k1())) +     "    "  )
        
        diff = np.abs(np.log10(kMinus) - np.log10(predictedD.k1())) 
        
        file.write(str( "%0.3g" % diff ) +     "    "  )
    
        file.write("\n")
    
        file.flush()


comparison() 


file.close()
