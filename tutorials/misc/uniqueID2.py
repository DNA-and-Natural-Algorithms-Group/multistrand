# FD, May 17th, 2017. 
# This demonstrates the pair-type functionality
# For a given complex, Pairtype returns a unique representation. 
# This is important for hashing functions 

import sys

from multistrand.objects import Complex, Domain, Strand, StopCondition
from multistrand.options import Options
from multistrand.system import SimSystem
from multistrand.utils import uniqueStateID, pairType
from multistrand.experiment import standardOptions, hybridization

    
ATIME_OUT = 0.0010

def doMain():    

    
    
    names1 = "strand1,strand2"
    names2 = "strand2,strand1"

    struct1 = "..((....((..+))....))"
    struct2 = "((....((+..))....)).."

    uID1 = pairType(names1, struct1)
    uID2 = pairType(names2, struct2)
    
    print uID1 
    print uID2
    
    if uID1 == uID2:
        print "Passed test 1"
    else:
        print "Fail test 1"            
    
    
    print "\n\n\n\n\n"
       
        
    names1 = '11:invader,5:invader*,10:top'
    names2 = '11:invader,5:invader*,10:top'

    
    struct1 = ".(.+.)..(((+)))."
    struct2 = ".(.+.)..((.+.))."
    
    uID1 = pairType(names1, struct1)
    print "\n\n\n\n"
    uID2 = pairType(names2, struct2)
#     
#     print names1 + struct1
    print uID1
#     
#     
#     
#     print names2 + struct2
    print uID2
#     
    
    if not uID1 == uID2:
        print "pass test2"
    else:
        print "Fail test 2"  

#     print "Equality: "  + str( uID1 == uID2 )
    
    
    
    
# # The actual main method
if __name__ == '__main__':
    
    print sys.argv
    doMain()
         
        
        
































