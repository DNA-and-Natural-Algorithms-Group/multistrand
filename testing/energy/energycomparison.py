from __future__ import print_function
import subprocess
import unittest
import time

import sys
sys.path.append('../../')

import multistrand.options
import multistrand.objects
import multistrand.system



class EnergyComparison( unittest.TestCase ):
    def __init__( self ):
        pass

    def setUp(self):
        o = multistrand.options.Options()
        multistrand.system.initialize_energy_model( o )
        self.sections = {}
        self.sectionname = None

    def tearDown(self):
        self.sections = {}
        self.sectionname = None
    
    def test_single_stranded( self ):
        self.process_file( 'testSetSS.txt' )
        for i in self.sections.keys():
            print("Section {0}:\n".format( i ))
            correct = len( [j for j in self.sections[i] if j[0]] )
            incorrect = [j for j in self.sections[i] if not j[0]]
            print("Correct: [{0}/{1}]\n".format( correct, len( self.sections[i] )))
            if len( incorrect ) > 0:
                print("Failed: \n")
                for k in incorrect:
                    print("{0[3]} {0[4]}: {0[1]} {0[2]}".format( k ))
            

    def process_file( self, filename ):
        lines = []
        f = open(filename)
        count = 0
        for l in f:
            if l.startswith('>'):
                self.sectionname = l[1:].rstrip('\n')
            else:
                lines.append(l.rstrip(' \n'))

            if len(lines) >= 2:
                self.process( lines[0], lines[1], count )
                lines[:] = []
                count += 1
            if count > 500:
                break
        f.close()
        
    def process( self, sequence, structure, count ):
        if (count % 10) == 0:
            print("{0} ".format(count) , end='\r')
            sys.stdout.flush()
        #print(" sequence: {0}\nstructure: {1}".format( sequence, structure ) )
        #return
        nupack_energy = self.nupack_ene( sequence, structure )
        multistrand_energy = self.multistrand_ene( sequence, structure )
        if self.sectionname not in self.sections:
            self.sections[self.sectionname] = []
        #print "Test: {0} == {1} \n".format( nupack_energy, multistrand_energy )
        self.sections[self.sectionname].append( (nupack_energy == multistrand_energy, nupack_energy, multistrand_energy, sequence, structure ))

    def nupack_ene( self, sequence, structure ):
        time.sleep(.05)
        p = subprocess.Popen(["energy", "-material", "dna"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        input_str = "{seq}\n{struc}\n".format( seq=sequence, struc=structure )
        
        # input_str = "{0}\n{1}\n{2}\n".format( len(self.strand_list),
        #                                 "\n".join( [i.sequence for i in self.strand_list] ),
        #                                 " ".join( [str(i+1) for i in range( len( self.strand_list ))])
        #                                   )
        stdout_chars = p.communicate(input_str)[0]
        flag = False
        energy = None
        # value we want is on the line after "% Energy"
        result = stdout_chars.split("\n")
        for l in result:
            if flag:
                energy = float( l )
                flag = False
            if l.startswith("% Energy"):
                flag = True
        return energy

    def multistrand_ene( self, sequence, structure ):
        c = multistrand.objects.Complex( sequence=sequence, structure=structure )
        energy = multistrand.system.energy([c])
        return round(energy[0],2)
    
