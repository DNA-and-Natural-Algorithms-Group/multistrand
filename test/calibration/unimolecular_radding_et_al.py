from __future__ import print_function
import unittest
import warnings
#import subprocess
#import timeit
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

idx = os.getcwd().find( os.path.join('test','calibration') )
if idx > -1:
    rootpath = os.path.abspath(os.getcwd()[:idx])
    #abspath cleans up the tail of the pathname as needed.
    if rootpath not in sys.path:
        sys.path.append( rootpath )
elif multihome != None:
    sys.path.append(multihome)

try:
    from multistrand.objects import Strand, Complex, Domain, StopCondition
    from multistrand.options import Options, Constants
    from multistrand.system import SimSystem, initialize_energy_model
    import multistrand.utils
except ImportError:
    # we want to tell the user how to fix this, but then reraise so it still fails. 
    print("Could not import Multistrand, please add it to your sys.path, or make sure that MULTISTRANDHOME is set correctly. This sub-program can also be run from the native test/calibration/ directory.")
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

class ThreewayBMTest( unittest.TestCase ):
    """ This test case class handles the comparison between Kinfold and Multistrand on
    a particular input file of random sequences. """

    def setUp( self ):
        """ Does nothing at the moment """
        pass

        
    def __getattr__(self, seq):
        parts = seq.split(':')
        structure = None
        if len(parts) == 9:
            idx = int( parts[0] )
            time = float(parts[1])
            file_prefix = parts[2]
            seq = parts[3].upper()
            dangles = parts[4]
            ratemethod = parts[5]
            temperature = float(parts[6])
            energymodel = parts[7]
            substratetype = parts[8]
        else:
            return super(ThreewayBMTest, self).__getattr__(name)

        class stub_test_runner(object):
            def __init__(self,runner, file_prefix):
                self.my_test_runner = runner
                self.prefix = file_prefix
            def __call__(self,*args):
                print("{0} ...".format(self.__doc__))
                self.my_test_runner( seq, idx, time, self.prefix, dangles, ratemethod, temperature, energymodel, substratetype )
                
        if len(seq) > 40:
            shortname = str(len(seq)) + ': ' + seq[:40] + '...'
        else:
            shortname = str(len(seq)) + ': ' + seq

        mystub = stub_test_runner( self.my_test_runner, file_prefix )
        mystub.__doc__ = "Random Sequence (#{1},t={2}) [{0}]".format(shortname,idx,time)
        # if structure != None and len(structure)<80:
        #     mystub.__doc__  += "\nStructure [{0}]".format(structure)
        return mystub
     

    def my_test_runner( self, seq, idx, time, prefix, dangles, ratemethod, temperature, energymodel, substratetype ):
        filename = prefix + 'len_{0}_sequence_{1:03}.out'.format( len(seq), idx )
        if os.path.isfile( filename ):
            print("File [{filename}] already exists, skipping test.".format(filename=filename))
            return

        count = 1000
        times_ms = self.setup_multistrand( filename, seq, time, count, dangles, ratemethod, temperature, energymodel, substratetype )

        print("Sequence {1} [{2}]: {0}".format(  times_ms[0][0] + times_ms[0][1], idx, len(seq)),sep="")
        f = open(filename,'wt')
        f.write(
            repr( {'Multistrand':times_ms[0],
                   'Multistrand_init':tuple([j-i for i,j in zip(times_ms[0], times_ms[1])]),
                   'maxtime':time,
                   'length':len(seq),
                   'load':os.getloadavg()}))
        f.close()

    @timer
    def setup_multistrand( self, filename, sequence, time, count, dangles, ratemethod, temperature, energymodel, substratetype ):
        input_o = self.helper_create_Multistrand_options( sequence, time, count, dangles, ratemethod, temperature, energymodel, substratetype )
        result_filename_txt = filename.replace('.out','.result.txt')
        result_filename_dat = filename.replace('.out','.result.dat')
        #initialize_energy_model( input_o )
        
        s = SimSystem( input_o )

        @timer
        def runOnce():
            s.start()
            self.output_multistrand = str(input_o.interface.results)
            f = open(result_filename_txt,'wt')
            f.write( self.output_multistrand )
            f.close()

            f = open( result_filename_dat, 'wb')
            cPickle.dump( input_o.interface.results, f, protocol=-1)
            cPickle.dump( input_o , f, protocol=-1)
            f.close()
            if input_o.verbosity > 0:
                print(self.output_multistrand)

        res = runOnce()[1]
        return res

    def helper_create_Multistrand_options( self, sequence, time, count, dangles, ratemethod, temperature, energymodel, substratetype ):
        """ helper """

        toehold = Domain(name='toehold', length=8, sequence='GCTGTAGC')

        bm_region = Domain( name='bm', length = len(sequence), sequence=sequence )

        base_strand = toehold + bm_region + toehold
        incumbent_strand = bm_region.C + toehold.C
        invading_strand = toehold.C + bm_region.C

        bm_complex = Complex( strands = [invading_strand, incumbent_strand, base_strand], structure = "(.+((+)))" )
        end_complex = Complex( strands = [invading_strand, incumbent_strand, base_strand], structure = "((+.(+)))" )
        disassoc_complex = Complex( strands = [invading_strand], structure="..")
        
        stop_conditions = []
        stop_conditions.append( StopCondition("forward_bm", [(end_complex,0,0)] ))
        stop_conditions.append( StopCondition("disassoc", [(disassoc_complex,2,0)] ))
        

        o = Options( simulation_mode="Normal",
                     substrate_type = substratetype,
                     parameter_type = energymodel,
                     num_sims = count,
                     sim_time = time,
                     biscale = 1.0,
                     uniscale = 1.0,
                     start_state = [bm_complex],
                     temperature = temperature,
                     verbosity = 0)

        # o.parameter_type = Constants.ENERGYMODEL_TYPE[energymodel]
        # o.substrate_type = Constants.SUBSTRATE_TYPE[substratetype]
        o.dangles = Constants.DANGLES[dangles]
        o.rate_method = Constants.RATEMETHOD[ratemethod]
        o.stop_conditions = stop_conditions

        return o
    
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

class Threeway_BM_Calibration( Multistrand_Suite_Base ):
    """ insert stuff here.

    Options that need to be set for different runs.

    Dangles
    Rate Method
    Temperature?
    Energy Model?

    """

    def __init__(self, lengths, times, dangles, ratemethod, temperature, energymodel, substrate_type, file_prefix="" ):
        if not os.path.isdir( file_prefix ):
            os.mkdir( file_prefix )

        self._suite = unittest.TestSuite()
        self._lengths = {}
        self._times = {}
        
        for l in lengths:
            self._lengths[str(l)] = []
        for l,t in zip(lengths, times):
            self._times[str(l)] = t


        self.generate_or_load_sequences(file_prefix)
        self.add_test_cases(file_prefix, dangles, ratemethod, temperature, energymodel, substrate_type)

    def add_test_cases( self, file_prefix, dangles, ratemethod, temperature, energymodel, substrate_type ):
        for k,v in self._lengths.iteritems():
            for i in range( len(v) ):
                sequence = self._lengths[k][i]
                if len( self._lengths[k][i]) != int(k):
                    raise ValueError("Actual sequence length doesn't match sequence length value.")

                self._suite.addTest( ThreewayBMTest("{idx}:{time}:{prefix}:{seq}:{dangle}:{rate}:{temp}:{em}:{st}".format(idx=i, time=self._times[k], prefix = file_prefix, seq= sequence, dangle=dangles, rate=ratemethod, temp=float(temperature), em=energymodel, st=substrate_type)))
                

    def generate_or_load_sequences(self,file_prefix):
        """ For every length we should be running, either create (if it doesn't exist) or load the sequence data from file. """
        
        for n in self._lengths.keys():
            if not os.path.isfile(file_prefix + 'length_{0}_sequences_random.dat'.format( n )):
                f = open(file_prefix + 'length_{0}_sequences_random.txt'.format(n),'wt')
                
                for i in range(1000):  # 1000 sequences for each length.
                    self._lengths[n].append( multistrand.utils.generate_sequence( int(n)))

                    f.write( "{0}\n".format(self._lengths[n][-1] ))
                f.close()
                f = open(file_prefix + 'length_{0}_sequences_random.dat'.format(n),'wb')
                cPickle.dump( self._lengths[n], f, protocol=-1)
                f.close()
            else:
                f = open(file_prefix + 'length_{0}_sequences_random.dat'.format(n),'rb')
                self._lengths[n] = cPickle.load( f )


    

# class Threeway_BM_Tests( Multistrand_Suite_Base ):
#     """ Uses the Speedtest_FromFile testcase class to run tests on
#     3-way branch migration systems."""

#     def __init__(self, lengths, times, file_prefix = ""):



#         self._suite.addTest( Speedtest_FromFile('{idx}:{time}:{prefix}:{seq}:{struc}'.format( idx=i, time=time_to_sim, prefix=file_prefix, seq=bm_complex.sequence, struc=bm_complex.structure)))

# WARNING: This is copied into analysis_unimolecular.py, any changes should go to both files
def list_of_parameters():

    dangles_types = ["None", "Some", "All"]
    rate_methods = ["Kawasaki", "Metropolis"]
    temperatures = [310.15, 273.15+25.0]
    models = ["Nupack","Vienna"]
    substrates = ["DNA","RNA"]
    # With these 5 variables, there are 48 total parameter sets we wish to calibrate for.
    # These are divided into groups of 6 which have the same temperature, model and substrate.

    parameters = []
    for e in substrates:
        for d in models:
            for c in temperatures:
                for b in rate_methods:
                    for a in dangles_types:
                        parameters.append( (a,b,c,d,e) )
    return parameters

# WARNING: This is copied into analysis_unimolecular.py, any changes should go to both files
def args_from_param_number( num ):
    """ Returns the proper argument list for Threeway_BM_Calibration, using the index num into the list of all parameter combinations """
    
    parameter_list = list_of_parameters()
    pm = parameter_list[num]
    args = ([10,20,30], [1.0e5, 1.0e6, 1.0e6]) + pm + ("calibration_{0[3]}_{0[4]}_{0[1]}_{0[0]}_{1}C/".format( pm, int(pm[2] - 273.15)),)
    return args

if __name__ == '__main__':
    if len(sys.argv) < 2:
        idx = 0
    else:
        idx = int(sys.argv[1])

    calibration = Threeway_BM_Calibration( *args_from_param_number( idx ))
    calibration.runTests_Async()

    #    calibration = Threeway_BM_Calibration( [10,20,30], [1.0e5,1.0e6,1.0e6],  "None", "Kawasaki", 310.15, "Nupack", "DNA", 'calibration_test/' )


