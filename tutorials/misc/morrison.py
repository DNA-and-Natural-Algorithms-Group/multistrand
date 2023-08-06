# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
Morrison and Stols, 2003. Comparison of Multistrand vs. reported values.
this computes the rates via k+ / k- = exp -dG / RT
"""

import os
import math
import xlrd         #excel sheets
import numpy as np

from multistrand.utils import GAS_CONSTANT, C2K

from anneal import compute
from nupack import pfunc


print("Morrison and Stols 1993")

def openDocument(document): 
    reader = xlrd.open_workbook(document, 'rb')
    row = reader.sheet_by_index(0)
    return row

def excelFind(row, col):
    dir = os.path.dirname(__file__)
    document = os.path.join(dir, 'data/helix/Fig6_0.xlsx')
    myDoc = openDocument(document)
    return myDoc.cell(row, col).value


## For each column, compute the rate and print it to a file. 
resultFileName = "morrison-results-ALT-METROPOLIS.txt"
file = open(resultFileName, 'w+')

file.write("Seq    temp    measured    predicted    pred-95-low    pred-95-high     alternative computation\n\n")

## FD: again, this will use the default settings of compute files

def doMorrison(myRange):

    diffSum = np.float(0.0)
    for i in myRange:

        seq = excelFind(i+1, 1)
        seqC = excelFind(i+1,2)
        temp = 1000/ float(excelFind(i+1, 3))
        measured = float(excelFind(i+1, 5))

        file.write(str(seq) + "   ")
        file.write(str( "%0.3g" %  (temp - C2K) ) + "   ")
        file.write(str(  "%0.3g" % measured) + "   ")

        predicted = compute(seq, temp)
        low, high = predicted.doBootstrap()

        dotparen = "("*len(seq) + "+" + ")"*len(seq)

        dG = pfunc([seq, seqC], [1,2], T=(temp - C2K), material="dna")
        print(str(dG))

        kMinus = predicted.k1() * math.exp( dG / ( GAS_CONSTANT * temp) )
        kMinusLow = low * math.exp( dG / ( GAS_CONSTANT * temp) )
        kMinusHigh = high * math.exp( dG / ( GAS_CONSTANT * temp) )
        
        file.write(str( "%0.3g" % np.log10(kMinus)) +     "    "  )
        file.write(str( "%0.3g" % np.log10(kMinusLow)) +     "    "  )
        file.write(str( "%0.3g" % np.log10(kMinusHigh)) +     "    "  )    
        
        diff = np.abs(np.log10(kMinus) - np.float(measured)) 
        diffSum = diffSum +  diff
        
        file.write(str( "%0.3g" % diff ) +     "    "  )
        file.write(str( "%0.3g" % diffSum ) +     "    "  )
        file.write("\n")
        file.flush()


doMorrison(list(range(6,12))) #do short seq first
doMorrison(list(range(6)))

file.close()
