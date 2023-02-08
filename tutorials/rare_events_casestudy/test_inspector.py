from __future__ import print_function

from multistrand.objects import Strand, Complex, Domain, StopCondition
from multistrand.options import Literals
from multistrand.system import SimSystem
from multistrand.experiment import makeComplex, hybridization, dissociation, standardOptions
from multistrand.builder import Builder


def test0(unusedargs):

    o1 = standardOptions();
    dissociation(o1, "CCCCCATTAAC");
    o1.simulation_mode = Literals.first_passage_time
    
    return o1


def test1(unusedargs):
    
    o1 = standardOptions();
    hybridization(o1, "CCC");
    o1.simulation_mode = Literals.first_passage_time
    
    return o1


def main():
        
    b = Builder(test0, [])
    b.genAndSavePathsFile(inspecting=True)
    
    print(b)
        

    b = Builder(test1, [])
    b.genAndSavePathsFile(inspecting=True)
    
    print(b)


main()

