## FD: Morrison and Stols, 2003. Comparison of Multistrand vs. reported values.

import sys, os
import xlrd         #excel sheets
import numpy as np

from dissociation import compute

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
resultFileName = "morrison-results.txt"
file = open(resultFileName, 'w+')

file.write("Seq    temp    measured    predicted    pred-95-low    pred-95-high\n\n") 

def doMorrison(myRange):

    for i in myRange:
        
        seq = excelFind(i+1, 1)
        temp = 1000/ float(excelFind(i+1, 3))
        measured = excelFind(i+1, 5) 
        
        
        file.write(str(seq) + "   ")
        file.write(str(temp) + "   ")
        file.write(str(measured) + "   ")
        
        predicted = compute(seq, temp)
        file.write(str(np.log10(predicted.k1())) +     "    "  )
    
        low, high = predicted.doBootstrap()
        file.write(str(np.log10(low)) +     "    "  )
        file.write(str(np.log10(high)) +     "    "  )    
    
        file.write("\n")
    
        file.flush()


doMorrison(range(6,12)) #do short seq first
doMorrison(range(6))

file.close()
rate = compute("AGAT")
