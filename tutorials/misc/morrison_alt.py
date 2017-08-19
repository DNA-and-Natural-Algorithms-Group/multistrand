## FD: Morrison and Stols, 2003. Comparison of Multistrand vs. reported values.

## FD: this is the EXACT same as morrison but computes the rates via k+ / k- = exp -dG / RT

GAS_CONSTANT_R = 0.0019872036

import sys, os
import math
import xlrd         #excel sheets
import numpy as np

from anneal import compute
from nupack import energy


print "Morrison and Stolls 1993"

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
resultFileName = "morrison-results-ALT.txt"
file = open(resultFileName, 'w+')

file.write("Seq    temp    measured    predicted    pred-95-low    pred-95-high     alternative computation\n\n") 

def doMorrison(myRange):

    for i in myRange:
        
        seq = excelFind(i+1, 1)
        seqC = excelFind(i+1,2)
        temp = 1000/ float(excelFind(i+1, 3))
        measured = excelFind(i+1, 5) 
        
        
        file.write(str(seq) + "   ")
        file.write(str(temp) + "   ")
        file.write(str(measured) + "   ")
        
        predicted = compute(seq, temp)
        low, high = predicted.doBootstrap()
        
        dotparen = "("*len(seq) + "+" + ")"*len(seq)
        
        dG = energy([seq, seqC], dotparen, T=(temp-273.15), material="dna") + 2.44  # NUPACK energy is too negative.
        print (str(dG)) 
        
        kMinus = predicted.k1() * math.exp( dG / ( GAS_CONSTANT_R * temp) ) 
        kMinusLow = low * math.exp( dG / ( GAS_CONSTANT_R * temp) ) 
        kMinusHigh = high * math.exp( dG / ( GAS_CONSTANT_R * temp) ) 
        
        file.write(str(np.log10(kMinus)) +     "    "  )
        file.write(str(np.log10(kMinusLow)) +     "    "  )
        file.write(str(np.log10(kMinusHigh)) +     "    "  )    
    
        file.write("\n")
    
        file.flush()


doMorrison(range(6,12)) #do short seq first
doMorrison(range(6))

file.close()
