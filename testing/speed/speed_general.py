from __future__ import print_function
import unittest
import warnings
import subprocess
import timeit
import os, os.path
import cPickle
import numpy
import random

# for IPython, some of the IPython libs used by unittest have a
# deprecated usage of BaseException, so we turn that specific warning
# off.
warnings.filterwarnings("ignore", r"BaseException[.]message has been deprecated as of Python 2[.]6", DeprecationWarning)

import multiprocessing
from multiprocessing import Pool

import sys
sys.path.append('../../')

try:
    from multistrand.objects import Strand, Complex
    from multistrand.options import Options, Constants
    from multistrand.system import SimSystem
    import multistrand.utils
except ImportError:
    # we want to tell the user how to fix this, but then reraise so it still fails. 
    print("Could not import Multistrand, please add it to your sys.path, or run this program from the native testing/speed/ directory.")
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


class Length_Result( dict ):
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
                    self[k] = v[4]

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
        
        

class Speedtest_FromFile( unittest.TestCase ):
    """ This test case class handles the comparison between Kinfold and Multistrand on
    a particular input file of random sequences. """

    def setUp( self ):
        print("setting up!")
        
    def __getattr__(self, seq):
        parts = seq.split(':')
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
        else:
            return super(Speedtest_FromFile, self).__getattr__(name)

        class stub_test_runner(object):
            def __init__(self,runner, file_prefix):
                self.my_test_runner = runner
                self.prefix = file_prefix
            def __call__(self,*args):
                print("{0} ...".format(self.__doc__))
                self.my_test_runner( seq, idx, time, self.prefix )
                
        if len(seq) > 40:
            shortname = str(len(seq)) + ': ' + seq[:40] + '...'
        else:
            shortname = str(len(seq)) + ': ' + seq

        mystub = stub_test_runner( self.my_test_runner, file_prefix )
        mystub.__doc__ = "Random Sequence (#{1},t={2}) [{0}]".format(shortname,idx,time)
        return mystub
            
        #raise AttributeError("{0}: Invalid name for a test case file.".format(name))

    def my_test_runner( self, seq, idx, time, prefix ):
        filename = prefix + 'len_{0}_sequence_{1}.out'.format( len(seq), idx )

        times_kin = self.setup_kinfold( seq, time, 100)
        times_ms = self.setup_multistrand( seq, time, 100 )

        print("Sequence {2} [{3}]: {0:>35} | {1}".format(  times_kin[0][2] + times_kin[0][3], times_ms[0][0] + times_ms[0][1], idx, len(seq)),sep="")
        f = open(filename,'wt')
        f.write( "Length_Result({0})".format(
            repr( {'Kinfold':times_kin[0],
                   'Multistrand':times_ms[0],
                   'Kinfold_init':tuple([j-i for i,j in zip(times_kin[0], times_kin[1])]),
                   'Multistrand_init':tuple([j-i for i,j in zip(times_ms[0], times_ms[1])]),
                   'maxtime':time,
                   'length':len(seq),
                   'load':os.getloadavg()} ) )
                 )
        f.close()

    def setUp(self):
        pass
    
    @timer
    def setup_kinfold( self, sequence, time, count ):
        kinfoldproc = subprocess.Popen(["Kinfold","--noShift","--logML","--start","--fpt","--time","{0:f}".format(time),"--num","{0:d}".format(count),"--silent","--dangle","0","--Par","dna.par"], stdin=subprocess.PIPE, stdout=subprocess.PIPE )

        input_str = "{0}\n{1}\n".format( sequence, "."*len(sequence) )

        @timer
        def runOnce():
            self.output_kinfold, _ = kinfoldproc.communicate( input_str )

        return runOnce()[1]

    @timer
    def setup_old_multistrand( self, sequence, time, count ):
        multistrandproc = subprocess.Popen(["Multistrand"],stdin=subprocess.PIPE, stdout=subprocess.PIPE )
        input_str = self.helper_create_multistrand_infile( sequence, time, count)

        @timer
        def runOnce():
            self.output_multistrand, _ = multistrandproc.communicate( input_str )

        return runOnce()[1]

    @timer
    def setup_multistrand( self, sequence, time, count ):
        input_o = self.helper_create_Multistrand_options( sequence, time, count)
        s = SimSystem( input_o )

        @timer
        def runOnce():
            s.start()
            self.output_multistrand = str(input_o.interface.results)

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

    def helper_create_Multistrand_options( self, sequence, time, count):
        """ helper """
        
        return Options( num_sims= count,
                        sim_time= time,
                        start_state= [Complex(sequence=sequence,
                                              structure = "." * len( sequence ) )],
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
        unittest.TextTestRunner( descriptions=0,verbosity=0,stream=DevNull()).run( testcase )
    


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
                
                self._suite.addTest( Speedtest_FromFile('{idx}:{time}:{prefix}:{seq}'.format( idx=i, time=time_to_sim, prefix=file_prefix, seq=self._lengths[k][i])))

        
if __name__ == '__main__':
    short_lengths = Length_Tests( range(20,100,2), [5000.0]*11 + [1000.0]*39, 'length_short/')
    long_lengths = Length_Tests( range(100,205,5), [1000.0] * 21, 'length_longs/')
    very_long_lengths = Length_Tests( range(210,310,10), [100.0] * 10, 'length_very_longs/')
    single_short = Length_Tests( [30], [5000.0], 'length_short/')
    very_long_lengths.runTests_Async()
    #long_lengths.runTests_Async()
    #single_short.runTests_Async()



# # testing area
# class alpha( unittest.TestCase ):
#     def setUp( self ):
#         if hasattr( self, 'name' ):
#             print("Name found")
#         else:
#             print("Noname")  # always this case
#         print("setting up!")
        
#     def __getattr__(self, name):
#         if name.endswith("_txt"):
#             def stub_test_runner():
#                 self.myrun( name )
#             stub_test_runner.__doc__ = "Random Sequences [{0}]".format( name.replace('_','.'))
#             return stub_test_runner
#         else:
#             raise AttributeError
#     def runfunct( self, name ):
#         return lambda: self.myrun( name )
#     def myrun( self, name):
#         self.name = name
#         print("Processing name: {0}".format( name ) )

# # # test code:       
# # import unittest

# a = alpha('blah_txt')

# suite = unittest.TestSuite()
# suite.addTest( a)
# unittest.TextTestRunner(verbosity=2).run( suite )
        
