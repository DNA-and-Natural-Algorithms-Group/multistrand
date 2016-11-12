from __future__ import print_function
import unittest
import warnings
import subprocess
import timeit
import os, os.path
import cPickle
#import numpy
import random

# for IPython, some of the IPython libs used by unittest have a
# deprecated usage of BaseException, so we turn that specific warning
# off.
warnings.filterwarnings("ignore", r"BaseException[.]message has been deprecated as of Python 2[.]6", DeprecationWarning)

import multiprocessing
from multiprocessing import Pool

import sys
multihome = None
if 'MULTISTRANDHOME' in os.environ:
    if not os.path.isfile( os.path.join( os.environ['MULTISTRANDHOME'], 'setup.py') ):
        warnings.warn( ImportWarning("Could not find the file 'setup.py' in your MULTISTRANDHOME [{0}]; this environment variable is possibly out of date or not referring to the new Mulistrand distribution."))
        multihome=None
    else:
        if os.environ['MULTISTRANDHOME'] not in sys.path:
            multihome= os.environ['MULTISTRANDHOME']     

idx = os.getcwd().find( os.path.join('testing','speed') )
if idx == -1:
    idx = os.getcwd().find( os.path.join('test','speed') )
if idx > -1:
    rootpath = os.path.abspath(os.getcwd()[:idx])
    #abspath cleans up the tail of the pathname as needed.
    if rootpath not in sys.path:
        sys.path.append( rootpath )
elif multihome != None:
    sys.path.append(multihome)

try:
    from multistrand.objects import Strand, Complex, Domain
    from multistrand.options import Options, Constants
    from multistrand.system import SimSystem
    import multistrand.utils
except ImportError:
    # we want to tell the user how to fix this, but then reraise so it still fails. 
    print("Could not import Multistrand, please add it to your sys.path, or make sure that MULTISTRANDHOME is set correctly. This sub-program can also be run from the native test/speed/ directory.")
    raise



#Useful decorator for timing purposes. 
def timer( f ):
    def res( *args, **kargs ):
        start_time = os.times()
        retval = f(*args, **kargs)
        end_time = os.times()
        time_delta = [(j-i) for i,j in zip( start_time, end_time )]
        return (retval,time_delta)
    return res


class LengthResult( dict ):
    def __init__(self, in_dictionary):
        dict.__init__( self, in_dictionary )
        
        if 'maxtime' in self:
            self.maxtime = self['maxtime']
        else:
            self.maxtime = 1000.0
        if 'length' in self:
            self.length = self['length']
        else:
            self.length = None
        if 'load' in self:
            self.load = self['load']
        else:
            self.load = None

        keyvals = [i for i in self.iteritems()]
        for k,v in keyvals:
            try:
                if len(v) == 5:
                    if k.startswith('Kinfold'):
                        userplussys = v[2] + v[3]
                        self[k+'_user_sys']  = userplussys
                    if k.startswith('Multistrand'):
                        userplussys = v[0] + v[1]
                        self[k+'_user_sys']  = userplussys
                    if numpy.allclose( userplussys, v[4], .01 ):
                        self[k+'_close'] = True
                    else:
                        self[k+'_close'] = False
                    self[k+'_computed'] = v[4]

            except TypeError:
                pass

    def __str__(self):
        res = "{fields[0]:>11} : {0[Kinfold]:<20}\n\
{fields[1]:>11} : {0[Kinfold_init]:<}\n\
{fields[2]:>11} : {0[Multistrand]:<20}\n\
{fields[1]:>11} : {0[Multistrand_init]:<}\n\
{fields[3]:>11} : {0.maxtime}".format(self, fields=['Kinfold','(init)','Multistrand','Max Time'])
        if self.length != None:
            res += "\n{0:>11} : {1.length}".format( 'Length', self )
        if self.load != None:
            res += "\n{0:>11} : {1.load}".format( 'Load', self )
        return res

    def __repr__(self):
        res = "LengthResult({"
        
        keyvals = [i for i in self.iteritems()]
        for k,v in keyvals:
            if '_' in k and not k.endswith('_init'):
                continue
            if not res.endswith('{'):
                res += ','
            res += "{0!r} : {1!r}".format(k,v)

        res += "})"
        return res
        

class Speedtest_FromFile( unittest.TestCase ):
    """ This test case class handles the comparison between Kinfold and Multistrand on
    a particular input file of random sequences. """

    def setUp( self ):
        """ Does nothing at the moment """
        pass

        
    def __getattr__(self, seq):
        parts = seq.split(':')
        structure = None
        if len(parts) == 2:
            idx = int( parts[0] )
            seq = parts[1].upper()
            file_prefix = ''
            time = 1000.0
        elif len(parts) == 3:
            idx = int( parts[0] )
            time = float(parts[1])
            file_prefix = ''
            seq = parts[2].upper()
        elif len(parts) == 4:
            idx = int( parts[0] )
            time = float(parts[1])
            file_prefix = parts[2]
            seq = parts[3].upper()
        elif len(parts) == 5:
            idx = int( parts[0] )
            time = float(parts[1])
            file_prefix = parts[2]
            seq = parts[3].upper()
            structure = parts[4]
        else:
            return super(Speedtest_FromFile, self).__getattr__(name)

        class stub_test_runner(object):
            def __init__(self,runner, file_prefix):
                self.my_test_runner = runner
                self.prefix = file_prefix
            def __call__(self,*args):
                print("{0} ...".format(self.__doc__))
                self.my_test_runner( seq, idx, time, self.prefix, structure )
                
        if len(seq) > 40:
            shortname = str(len(seq)) + ': ' + seq[:40] + '...'
        else:
            shortname = str(len(seq)) + ': ' + seq

        mystub = stub_test_runner( self.my_test_runner, file_prefix )
        mystub.__doc__ = "Random Sequence (#{1},t={2}) [{0}]".format(shortname,idx,time)
        if structure != None and len(structure)<80:
            mystub.__doc__  += "\nStructure [{0}]".format(structure)
        return mystub
            
        #raise AttributeError("{0}: Invalid name for a test case file.".format(name))

    def my_test_runner( self, seq, idx, time, prefix, structure ):
        filename = prefix + 'len_{0}_sequence_{1}.out'.format( len(seq), idx )
        if os.path.isfile( filename ):
            print("File [{filename}] already exists, skipping test.".format(filename=filename))
            return
        times_kin = self.setup_kinfold( seq, time, 100, structure)
        times_ms = self.setup_multistrand( seq, time, 100, structure )

        print("Sequence {2} [{3}]: {0:>35} | {1}".format(  times_kin[0][2] + times_kin[0][3], times_ms[0][0] + times_ms[0][1], idx, len(seq)),sep="")
        f = open(filename,'wt')
        f.write( repr( LengthResult(
            {'Kinfold':times_kin[0],
             'Multistrand':times_ms[0],
             'Kinfold_init':tuple([j-i for i,j in zip(times_kin[0], times_kin[1])]),
             'Multistrand_init':tuple([j-i for i,j in zip(times_ms[0], times_ms[1])]),
             'maxtime':time,
             'length':len(seq),
             'load':os.getloadavg()})))
        f.close()

    @timer
    def setup_kinfold( self, sequence, time, count ):
        kinfoldproc = subprocess.Popen(["Ksim","--noShift","--logML","--start","--fpt","--time","{0:f}".format(time),"--num","{0:d}".format(count),"--silent","--dangle","0","--Par","dna.par"], stdin=subprocess.PIPE, stdout=subprocess.PIPE )

        sequence = sequence.replace("T","U")
        if structure == None:
            input_str = "{0}\n{1}\n".format( sequence, "."*len(sequence) )
        else:
            input_str = "{0}\n{1}\n".format( sequence,structure)                                         

        @timer
        def runOnce():
            self.output_kinfold, _ = kinfoldproc.communicate( input_str )

        res = runOnce()[1]
        return res

    @timer
    def setup_old_multistrand( self, sequence, time, count, structure ):
        multistrandproc = subprocess.Popen(["Multistrand"],stdin=subprocess.PIPE, stdout=subprocess.PIPE )
        input_str = self.helper_create_multistrand_infile( sequence, time, count)

        @timer
        def runOnce():
            self.output_multistrand, _ = multistrandproc.communicate( input_str )

        return runOnce()[1]

    @timer
    def setup_multistrand( self, sequence, time, count, structure ):
        input_o = self.helper_create_Multistrand_options( sequence, time, count, structure)
        s = SimSystem( input_o )

        @timer
        def runOnce():
            s.start()
            self.output_multistrand = str(input_o.interface.results)
            print(self.output_multistrand)

        res = runOnce()[1]
        return res
    
    def helper_create_multistrand_infile(self,sequence,time,count):
        """ reuses some code from strandedit to create Multistrand input files. """
        infile = '#Strands\n' \
                 'Strand1, {0}\n'.format( sequence )
        infile += '#StartStructure\n'\
                  'Strand1\n'\
                  '{0}\n'.format( "." * len(sequence))

        #Other options.
        infile += '#StopOptions=0\n' \
                  '#Energymodel=VIENNA_DNA\n##\n' \
                  '#OutputInterval=-1\n' \
                  '#Temperature=37\n' \
                  '#NumSims={0:d}\n'.format( count )
        infile += '#SimTime={0:f}\n'.format( time )
        infile += '#Ratemethod=2\n'  \
                  '#dangles=0\n' \
                  '#Inter=1.0\n' \
                  '#Intra=1.0\n' \
                  '##'
        return infile

    def helper_create_Multistrand_options( self, sequence, time, count, struc):
        """ helper """
        
        return Options( num_sims= count,
                        sim_time= time,
                        start_state= [Complex(sequence=sequence,
                                              structure = struc or "." * len( sequence ) )],
                        dangles = 'None',
                        biscale = 1.0,
                        uniscale = 1.0,
                        substrate_type = 'DNA',
                        parameter_type = 'Vienna' )

class Multistrand_Suite_Base( object ):
    """ Base class for test suites - defines async run, etc. """
    def runTests_Async(self, shuffle_tasks = True):
        """ This runner runs the tests via an asynchronous pool of
        workers.

        The number of workers is equal to the CPU core count, and this
        function only returns once every task has been completed. It's
        possible we should up the chunksize so that tasks are
        partitioned in bigger groups, but since all the cores being
        used are local, it seems like overkill."""
        starttime = os.times()
        k = multiprocessing.cpu_count()
        p = Pool( processes = k )

        if shuffle_tasks:
            random.shuffle( self._suite._tests )
        
        p.map( MyRunner , iter(self._suite), chunksize = 1 )
        p.close()
        p.join()
        endtime = os.times()
        print("Async run complete! Processing took [{0[4]}] of real time before completion. [u/s/cu/cs]:[{0[0]}/{0[1]}/{0[2]}/{0[3]}]".format( [j-i for i,j in zip( starttime, endtime )] ))
    

class MyRunner( object ):
    def __init__( self, testcase ):
        class DevNull(object):
            def write(self, _): pass
        #unittest.TextTestRunner( descriptions=0,verbosity=0,stream=DevNull()).run( testcase )
        unittest.TextTestRunner( descriptions=1,verbosity=1).run( testcase )
    

class Length_Tests( Multistrand_Suite_Base ):
    """ Uses the Speedtest_FromFile testcase class to run a bunch of test cases. """
    def __init__(self, lengths, times, file_prefix = ""):
        self._suite = unittest.TestSuite()
        self._lengths = {}
        self._times = {}
        
        for l in lengths:
            self._lengths[str(l)] = []
        for l,t in zip(lengths, times):
            self._times[str(l)] = t

        for n in self._lengths.keys():
            if not os.path.isfile(file_prefix + 'length_{0}_sequences_random.dat'.format( n )):
                f = open(file_prefix + 'length_{0}_sequences_random.txt'.format(n),'wt')
                for i in range(100):
                    self._lengths[n].append( multistrand.utils.generate_sequence( int(n)))
                    f.write( "{0}\n".format(self._lengths[n][-1] ))
                f.close()
                f = open(file_prefix + 'length_{0}_sequences_random.dat'.format(n),'wb')
                cPickle.dump( self._lengths[n], f, protocol=-1)
                f.close()
            else:
                f = open(file_prefix + 'length_{0}_sequences_random.dat'.format(n),'rb')
                self._lengths[n] = cPickle.load( f )

        for k,v in self._lengths.iteritems():
            for i in range( len(v) ):
                if k in self._times:
                    time_to_sim = self._times[k]
                else:
                    time_to_sim = (len(self._lengths[k][i]) <= 40 and 5000.0) or 1000.0
                if i > 2:
                    break
                self._suite.addTest( Speedtest_FromFile('{idx}:{time}:{prefix}:{seq}'.format( idx=i, time=time_to_sim, prefix=file_prefix, seq=self._lengths[k][i])))

        
if __name__ == '__main__':
    pass
    #short_lengths = Length_Tests( range(20,100,2), [5000.0]*11 + [1000.0]*39, 'length_short/')
    #long_lengths = Length_Tests( range(100,205,5), [1000.0] * 21, 'length_longs/')
    #very_long_lengths = Length_Tests( range(210,310,10), [100.0] * 10, 'length_very_longs/')
    #single_short = Length_Tests( [30], [5000.0], 'newblock/')
    #del single_short._suite._tests[3:]
    #single_long = Length_Tests( [105], [1000.0], 'length_longs/')
    #single_long.runTests_Async(shuffle_tasks=False)
    #very_long_lengths.runTests_Async()
    #long_lengths.runTests_Async()
    #single_short.runTests_Async()
class Fourway_BM_Tests( Multistrand_Suite_Base ):
    """ Uses the Speedtest_FromFile testcase class to run tests on
    4-way branch migration systems."""

    def __init__(self, lengths, times, file_prefix = ""):
        self._suite = unittest.TestSuite()
        self._lengths = {}
        self._times = {}

        for l in lengths:
            self._lengths[str(l)] = []
        for l,t in zip(lengths, times):
            self._times[str(l)] = t

        a = Domain(name='a', length=4, sequence='GTTC')
        b = Domain(name='b', length=4, sequence='GCCC')
        c = Domain(name='c', length=4, sequence='TCAC')
        d = Domain(name='d', length=4, sequence='AAGG')
        T = Domain(name='T', length=3, sequence='AAA')
        L = Domain(name='L')
        
        primary_strand = a + L + b + T + b.C + L.C + c + T + c.C + L + d + T + d.C + L.C + a.C
        primary_strand.name = "Primary"

        for n in self._lengths.keys():
            if n > 0 and not os.path.isfile(file_prefix + 'length_{0}_sequences_random.dat'.format( n )):
                f = open(file_prefix + 'length_{0}_sequences_random.txt'.format(n),'wt')

                for i in range(100):
                    self._lengths[n].append( multistrand.utils.generate_sequence( int(n)))

                    f.write( "{0}\n".format(self._lengths[n][-1] ))
                f.close()
                f = open(file_prefix + 'length_{0}_sequences_random.dat'.format(n),'wb')
                cPickle.dump( self._lengths[n], f, protocol=-1)
                f.close()
            elif n == 0:
                for i in range(100):
                    self._lengths[n].append("")
            else:
                f = open(file_prefix + 'length_{0}_sequences_random.dat'.format(n),'rb')
                self._lengths[n] = cPickle.load( f )

        bm_complex = Complex( strands = [primary_strand], structure = "(((.))(.)((.)))")
        for k,v in self._lengths.iteritems():
            L.length = k
            for i in range( len(v) ):
                L.sequence = self._lengths[k][i]
                if i == 0:
                    bm_complex._init_parse_structure("(((.))(.)((.)))")
                if k in self._times:
                    time_to_sim = self._times[k]
                else:
                    time_to_sim = (len(self._lengths[k][i]) <= 40 and 5000.0) or 1000.0
                if len(bm_complex.sequence) != len(bm_complex.structure):
                    raise ValueError("Sequence mismatch.")

                self._suite.addTest( Speedtest_FromFile('{idx}:{time}:{prefix}:{seq}:{struc}'.format( idx=i, time=time_to_sim, prefix=file_prefix, seq=bm_complex.sequence, struc=bm_complex.structure)))

class Threeway_BM_Tests( Multistrand_Suite_Base ):
    """ Uses the Speedtest_FromFile testcase class to run tests on
    3-way branch migration systems."""

    def __init__(self, lengths, times, file_prefix = ""):
        self._suite = unittest.TestSuite()
        self._lengths = {}
        self._times = {}
        
        for l in lengths:
            self._lengths[str(l)] = []
        for l,t in zip(lengths, times):
            self._times[str(l)] = t

        a = Domain(name='a', length=4, sequence='GTTC')
        b = Domain(name='b', length=4, sequence='GCCC')

        T = Domain(name='T', length=3, sequence='AAA')
        L = Domain(name='L')

        primary_strand = L + a + T + a.C + L.C + b + T + b.C + L
        primary_strand.name = "Primary"
        
        for n in self._lengths.keys():
            if n > 0 and not os.path.isfile(file_prefix + 'length_{0}_sequences_random.dat'.format( n )):
                f = open(file_prefix + 'length_{0}_sequences_random.txt'.format(n),'wt')
                
                for i in range(100):
                    self._lengths[n].append( multistrand.utils.generate_sequence( int(n)))

                    f.write( "{0}\n".format(self._lengths[n][-1] ))
                f.close()
                f = open(file_prefix + 'length_{0}_sequences_random.dat'.format(n),'wb')
                cPickle.dump( self._lengths[n], f, protocol=-1)
                f.close()
            elif n == 0:
                for i in range(100):
                    self._lengths[n].append("")
            else:
                f = open(file_prefix + 'length_{0}_sequences_random.dat'.format(n),'rb')
                self._lengths[n] = cPickle.load( f )
        bm_complex = Complex( strands = [primary_strand], structure = "((.))(.).")
        for k,v in self._lengths.iteritems():
            L.length = k
            for i in range( len(v) ):
                L.sequence = self._lengths[k][i]
                if i == 0:
                    bm_complex._init_parse_structure("((.))(.).")
                if k in self._times:
                    time_to_sim = self._times[k]
                else:
                    time_to_sim = (len(self._lengths[k][i]) <= 15 and 5000.0) or 1000.0
                if len(bm_complex.sequence) != len(bm_complex.structure):
                    raise ValueError("Sequence mismatch.")
                
                self._suite.addTest( Speedtest_FromFile('{idx}:{time}:{prefix}:{seq}:{struc}'.format( idx=i, time=time_to_sim, prefix=file_prefix, seq=bm_complex.sequence, struc=bm_complex.structure)))

        
if __name__ == '__main__':
    #pass
    short_lengths = Length_Tests( range(20,100,2), [5000.0]*11 + [1000.0]*39, 'length_short/')
    short_lengths.runTests_Async()
    #long_lengths = Length_Tests( range(100,205,5), [1000.0] * 21, 'length_longs/')
    # very_long_lengths = Length_Tests( range(210,310,10), [100.0] * 10, 'length_very_longs/')
    #single_short = Length_Tests( [30], [5000.0], 'length_short/')
    #bm_4way_tests = Fourway_BM_Tests( range(0,42,2), [5000.0]*10+[1000.0]*11, 'short_4way/')
    #bm_3way_tests = Threeway_BM_Tests( range(0,60), [5000.0]*10+[1000.0]*50, 'short_3way/' )
    #very_long_lengths.runTests_Async()
    #short_lengths.runTests_Async()
    #long_lengths.runTests_Async()
    #single_short.runTests_Async(False)
    #bm_3way_tests.runTests_Async()


