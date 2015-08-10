####################################################################
#                                                                  #
#  Copyright (c) 2010-2015 California Institute of Technology.     #
#  Distributed under the MIT License.                              #
#  (See accompanying file LICENSE or copy at                       #
#  http://opensource.org/licenses/MIT)                             #
#                                                                  #
####################################################################
#                                                                  #
#   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)         #
#                                                                  #
####################################################################

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
        self.display_results( self.sections )
        self.write_results( self.sections, 'testSS_nupack_multistrand.results')
    def display_results( self, sections_data ):
        for i in sections_data.keys():
            print("Section {0}:\n".format( i ))
            correct = len( [j for j in sections_data[i] if j[0]] )
            incorrect = [j for j in sections_data[i] if not j[0]]
            print("Correct: [{0}/{1}]\n".format( correct, len( sections_data[i] )))
            incorrect_rounding = 0
            rounding_mismatch = 0
            incorrect_other = 0
            if len( incorrect ) > 0:
                print("Failed: \n")
                for k in incorrect:
                    if round(k[1],2) == round(k[2],2):
                        incorrect_rounding += 1
                    elif -.01 <= round(k[1]-k[2],2) <= .01:
                        rounding_mismatch += 1
                    elif len(incorrect) < 50:
                        print("{0[3]} {0[4]}: {0[1]} {0[2]}".format( k ))
                    else:
                        incorrect_other += 1
            if incorrect_rounding > 0:
                print("Incorrect (Rounding): [{0}/{1}]".format( incorrect_rounding, len( sections_data[i] )))
            if rounding_mismatch > 0:
                print("Incorrect (Mismatch Rounding): [{0}/{1}]".format( rounding_mismatch, len( sections_data[i] )))
            if incorrect_other > 0:
                print("Incorrect (Other): [{0}/{1}]".format( incorrect_other, len( sections_data[i] )))

    def write_results( self, sections_data, fname=None):
        f = open(fname+'.txt','wt')
        f.write( repr( sections_data ))
        f.close()
        import cPickle
        f = open(fname + ".dat", 'wb')
        cPickle.dump( sections_data, f, -1 )
        f.close()

    def process_file( self, filename ):
        lines = []
        f = open(filename)
        count = 0
        for l in f:
            if l.startswith('>'):
                if self.sectionname != None:
                    print("Section {0} Complete [{1}]".format( self.sectionname, len(self.sections[self.sectionname] ) ))
                self.sectionname = l[1:].rstrip('\n')
            else:
                lines.append(l.rstrip(' \n'))

            if len(lines) >= 2:
                self.process( lines[0], lines[1], count )
                lines[:] = []
                count += 1
#            if count > 500:
#                break
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
        time.sleep(.02)
        p = subprocess.Popen(["energy", "-material", "dna"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        input_str = "{seq}\n{struc}\n".format( seq=sequence, struc=structure )
        
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
        return round(energy,2)

    def multistrand_ene( self, sequence, structure ):
        c = multistrand.objects.Complex( sequence=sequence, structure=structure )
        energy = multistrand.system.energy([c])
        return round(energy[0],2)
    
class Vienna_Comparison( EnergyComparison ):
    def __init__( self ):
        pass

    def setUp(self):
        self.vienna_proc = subprocess.Popen(["RNAeval", "-d0","-P","dna.par",'-logML'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        o = multistrand.options.Options()
        o.dangles = 0
        o.parameter_type = multistrand.options.Constants.ENERGYMODEL_TYPE['Vienna']
        
        multistrand.system.initialize_energy_model( o )
        self.sections = {}
        self.sectionname = None

    def tearDown(self):
        self.sections = {}
        self.sectionname = None
        v.vienna_proc.kill()
    
    def test_single_stranded( self ):
        self.process_file( 'testSetSS.txt' )
        self.display_results( self.sections )
        self.write_results( self.sections, 'testSS_vienna_multistrand.results')
        
    def process( self, sequence, structure, count ):
        if (count % 10) == 0:
            print("{0} ".format(count) , end='\r')
            sys.stdout.flush()

        vienna_energy = self.vienna_ene( sequence, structure )
        multistrand_energy = self.multistrand_ene( sequence, structure )
        if self.sectionname not in self.sections:
            self.sections[self.sectionname] = []

        self.sections[self.sectionname].append( (vienna_energy == multistrand_energy, vienna_energy, multistrand_energy, sequence, structure ))

    def vienna_ene( self, sequence, structure ):
        if self.vienna_proc == None:
            raise ValueError

        input_str = "{seq}\n{struc}\n".format( seq=sequence, struc=structure )
        self.vienna_proc.stdin.write(input_str)
        self.vienna_proc.stdin.flush()
        count = 0
        l = self.vienna_proc.stdout.readline()
        while not (l.startswith(" energy") or (l.startswith( structure ) and len(l)>len(structure))):
            l = self.vienna_proc.stdout.readline()
            count += 1

        if l.startswith(" energy"):
            energy = float(l[l.rfind('=')+1:])
        else:
            energy = float(l[l.rfind('(')+1:l.rfind(')')])
        return round(energy,2)


    # def multistrand_ene( self, sequence, structure ):
    #     c = multistrand.objects.Complex( sequence=sequence, structure=structure )
    #     energy = multistrand.system.energy([c])
    #     return round(energy[0],2)
    
