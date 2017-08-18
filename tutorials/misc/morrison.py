## FD: Morrison and Stols, 2003. Comparison of Multistrand vs. reported values.

import sys, os
import xlrd         #excel sheets

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

for i in range(12):
    
    name = excelFind(i+1, 1)
    seq = excelFind(i+1, 2)
    
    file.write(str(name) + "    ")
    file.write(str(seq) + "   ")
    file.write(str(99.0) + "   ")

    file.write("\n")


file.close()
rate = compute("AGAT")
