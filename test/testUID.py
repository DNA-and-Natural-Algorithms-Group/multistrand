# FD, May 17th, 2017. 
# This demonstrates the pair-type functionality
# For a given complex, Pairtype returns a unique representation. 
# This is important for hashing functions 

import sys 
import os.path


from multistrand.objects import Complex, Domain, Strand, StopCondition
from multistrand.options import Options
from multistrand.system import SimSystem
from multistrand.utils import uniqueStateID, pairType
from multistrand.experiment import standardOptions, hybridization

import unittest
import warnings
    

class idTest(unittest.TestCase):    

    def testOne(self):
        names1 = "strand1,strand2"
        names2 = "strand2,strand1"
    
        struct1 = "..((....((..+))....))"
        struct2 = "((....((+..))....)).."
    
        uID1 = pairType(names1, struct1)
        uID2 = pairType(names2, struct2)
        
        print(str(uID1))
        print(str(uID2))
        
        self.assertTrue(uID1 == uID2)
        
    def testTwo(self):
        names1 = '11:invader,5:invader*,10:top'
        names2 = '11:invader,5:invader*,10:top'

        struct1 = ".(.+.)..(((+)))."
        struct2 = ".(.+.)..((.+.))."
        
        uID1 = pairType(names1, struct1)
        uID2 = pairType(names2, struct2)
        
        print(str(uID1))
        print(str(uID2))
        
        self.assertTrue(not uID1 == uID2)

    
# # The actual main method
if __name__ == '__main__':
    test = idTest()
    test.doTest()
