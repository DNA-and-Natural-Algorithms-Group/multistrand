from __future__ import print_function
import unittest

from multistrand.options import Options
from multistrand.system import SimSystem, initialize_energy_model, energy
from multistrand.objects import Complex, Strand

import sys
import nupack

#magic constantm, see multistrand.system.options
Complex_Energy = 2

class EnergyComparison(unittest.TestCase):
    
    def setUp(self):
        self.o = Options()
        initialize_energy_model(self.o)
        self.sections = {}
        self.sectionname = None
    
    
    def test_single_stranded(self):
        self.process_file('testSetSS.txt')
        
    def process_file(self, filename):
        lines = []
        f = open(filename)
        count = 0
        for l in f:
             if l.startswith('>'):
                if self.sectionname != None:
                    print("Section {0} Complete [{1}]".format(self.sectionname, len(self.sections[self.sectionname])))
                self.sectionname = l[1:].rstrip('\n')
             else:
                 lines.append(l.rstrip(' \n'))
             if len(lines) >= 2:
                self.process(lines[0], lines[1], count)
                lines[:] = []
                count += 1

        f.close()
     
    def process(self, sequence, structure, count):
        
        if (count % 1) == 0:
            print("{0} ".format(count) , end='\r')
            sys.stdout.flush()
            print(" sequence: {0}\nstructure: {1}".format(sequence, structure))  

        multistrand_energy = self.multistrand_ene(sequence, structure)
        nupack_energy = self.nupack_ene(sequence, structure)
        
        if self.sectionname not in self.sections:
            self.sections[self.sectionname] = []

        print("Test: {0} == {1} \n".format(nupack_energy, multistrand_energy))                 
        self.sections[self.sectionname].append((nupack_energy == multistrand_energy, nupack_energy, multistrand_energy, sequence, structure))
        
    def nupack_ene(self, sequence, structure):
        
        result = nupack.energy([sequence], structure, material='dna')        
        return round(result, 2)

    def multistrand_ene(self, SEQUENCE, STRUCTURE):
              
#         c = multistrand.objects.Complex(strands=[Strand(name="hairpin", sequence=SEQUENCE)], structure=STRUCTURE)
        c = Complex(strands=[Strand(name="hairpin", sequence=SEQUENCE)], structure=STRUCTURE)

        myE = energy([c], self.o, Complex_Energy)  
        return round(myE[0], 2)
        

