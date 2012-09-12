from __future__ import print_function
import unittest
import warnings

import os, os.path
import cPickle

import random

# for IPython, some of the IPython libs used by unittest have a
# deprecated usage of BaseException, so we turn that specific warning
# off.
warnings.filterwarnings("ignore", r"BaseException[.]message has been deprecated as of Python 2[.]6", DeprecationWarning)

import multiprocessing
from multiprocessing import Pool


## Set Multistrand paths to be able to import multistrand module.

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
  

# class LengthResult( dict ):
#     def __init__(self, in_dictionary):
#         dict.__init__( self, in_dictionary )
        
#         if 'maxtime' in self:
#             self.maxtime = self['maxtime']
#         else:
#             self.maxtime = 1000.0
#         if 'length' in self:
#             self.length = self['length']
#         else:
#             self.length = None
#         if 'load' in self:
#             self.load = self['load']
#         else:
#             self.load = None

#         keyvals = [i for i in self.iteritems()]
#         for k,v in keyvals:
#             try:
#                 if len(v) == 5:
#                     if k.startswith('Kinfold'):
#                         userplussys = v[2] + v[3]
#                         self[k+'_user_sys']  = userplussys
#                     if k.startswith('Multistrand'):
#                         userplussys = v[0] + v[1]
#                         self[k+'_user_sys']  = userplussys
#                     if numpy.allclose( userplussys, v[4], .01 ):
#                         self[k+'_close'] = True
#                     else:
#                         self[k+'_close'] = False
#                     self[k+'_computed'] = v[4]

#             except TypeError:
#                 pass

#     def __str__(self):
#         res = "{fields[0]:>11} : {0[Kinfold]:<20}\n\
# {fields[1]:>11} : {0[Kinfold_init]:<}\n\
# {fields[2]:>11} : {0[Multistrand]:<20}\n\
# {fields[1]:>11} : {0[Multistrand_init]:<}\n\
# {fields[3]:>11} : {0.maxtime}".format(self, fields=['Kinfold','(init)','Multistrand','Max Time'])
#         if self.length != None:
#             res += "\n{0:>11} : {1.length}".format( 'Length', self )
#         if self.load != None:
#             res += "\n{0:>11} : {1.load}".format( 'Load', self )
#         return res

#     def __repr__(self):
#         res = "LengthResult({"
        
#         keyvals = [i for i in self.iteritems()]
#         for k,v in keyvals:
#             if '_' in k and not k.endswith('_init'):
#                 continue
#             if not res.endswith('{'):
#                 res += ','
#             res += "{0!r} : {1!r}".format(k,v)

#         res += "})"
#         return res

class Bimolecular_Calibration_TestCase( unittest.TestCase ):
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
            length = int(parts[3])
            dangles = parts[4]
            ratemethod = parts[5]
            temperature = float(parts[6])
            energymodel = parts[7]
            substratetype = parts[8]
        else:
            return super(Bimolecular_Calibration_TestCase, self).__getattr__(name)

        class stub_test_runner(object):
            def __init__(self,runner, file_prefix):
                self.my_test_runner = runner
                self.prefix = file_prefix
            def __call__(self,*args):
                print("{0} ...".format(self.__doc__))
                self.my_test_runner( length, idx, time, self.prefix, dangles, ratemethod, temperature, energymodel, substratetype )

        mystub = stub_test_runner( self.my_test_runner, file_prefix )
        mystub.__doc__ = "Random Trial (#{1},t={2}) [Length {0}]".format(length,idx,time)
        return mystub
     

    def my_test_runner( self, length, idx, time, prefix, dangles, ratemethod, temperature, energymodel, substratetype ):
        filename = prefix + 'bimolecular_len_{0}_run_{1:03}.out'.format( length, idx )
        if os.path.isfile( filename ):
            print("File [{filename}] already exists, skipping test.".format(filename=filename))
            return

        count = 100
        times_ms = self.setup_multistrand( filename, length, time, count, dangles, ratemethod, temperature, energymodel, substratetype )

        print("Trial {1} [{2}]: {0}".format(  times_ms[0][0] + times_ms[0][1], idx, length),sep="")
        f = open(filename,'wt')
        f.write(
            repr( {'Multistrand':times_ms[0],
                   'Multistrand_init':tuple([j-i for i,j in zip(times_ms[0], times_ms[1])]),
                   'maxtime':time,
                   'length':length,
                   'load':os.getloadavg()}))
        f.close()

    @timer
    def setup_multistrand( self, filename, length, time, count, dangles, ratemethod, temperature, energymodel, substratetype ):
        input_o = self.helper_create_Multistrand_options( length, time, count, dangles, ratemethod, temperature, energymodel, substratetype )
        result_filename_txt = filename.replace('.out','.result.txt')
        result_filename_dat = filename.replace('.out','.result.dat')

        # need to run init here as our test suite is very likely to be
        # across multiple energy model parameters. EDIT: Actually,
        # with async running, this is unclear - depends on where the
        # model gets forked / whether each new process has its own
        # energy model. Should be able to check quickly by analysis of
        # the memory footprint of 2 async runs with differing energy
        # models (preferably temperature difference) vs two separate
        # runs.
        
        initialize_energy_model( input_o )
        
        s = SimSystem( input_o )

        @timer
        def runOnce():
            s.start()
            self.output_multistrand = str(input_o.interface.results)
            f = open(result_filename_txt,'wt')
            f.write( input_o.calibration_string + "\n\n")
            f.write( self.output_multistrand )
            f.close()

            f = open( result_filename_dat, 'wb')
            cPickle.dump( input_o.interface.results, f, protocol=-1)
            cPickle.dump( input_o.calibration_string, f, protocol = -1 )
            #cPickle.dump( input_o , f, protocol=-1)
            f.close()
            if input_o.verbosity > 0:
                print(self.output_multistrand)

        res = runOnce()[1]
        return res

    def helper_create_Multistrand_options( self, length, time, count, dangles, ratemethod, temperature, energymodel, substratetype ):
        """ helper """


        ### Sequences taken from Morrison & Stols 1993:
        ### "The complementary 10-mers,d(TTGGTGATCC) and
        ### d(GGATCACCAA), and the complementary 20-mers,
        ### d(AGATTAGCAGGTTTCCCACC) and d(GGTGGGAAACCTGCTAATCT),
        ### were..."

        sequence_10 = Domain(name="hybrid_10",length=10, sequence='TTGGTGATCC')
        sequence_20 = Domain(name="hybrid_20",length=20, sequence='AGATTAGCAGGTTTCCCACC')
        if length == 10:
            seq = sequence_10
        else:
            seq = sequence_20
        
        strand_A = Strand(name="hybridization_A", domains=[seq])
        strand_B = strand_A.C  # the complementary strand to hybridize to.

        start_complex_A = Complex(strands=[strand_A],structure=".")
        start_complex_B = Complex(strands=[strand_B],structure=".")
        ## to enable initial state boltzmann sampling, uncomment the following lines:
        # start_complex_A.boltzmann_sample = True
        # start_complex_A.boltzmann_count = 1000  # we expect to run 1000 trajectories, so might as well generate the samples all at the start as it's faster.

        # start_complex_B.boltzmann_sample = True
        # start_complex_B.boltzmann_count = 1000  # we expect to run 1000 trajectories, so might as well generate the samples all at the start as it's faster.

        ## End boltzmann sampling setup.
        
        complete_complex = Complex(strands=[strand_A,strand_B], structure = "(+)")
         ## fully paired

        complete_stopcondition = StopCondition("complete",[(complete_complex,0,0)])
         ## exact structure stop condition, e.g. fully hybridized

        touching_stopcondition = StopCondition("touching",[(complete_complex,2,0)])
        
        o = Options( simulation_mode="First Passage Time",
                     substrate_type = substratetype,
                     parameter_type = energymodel,
                     dangles = dangles,
                     rate_method = ratemethod,
                     num_sims = count,
                     sim_time = time,
                     start_state = [start_complex_A,start_complex_B],
                     temperature = temperature,
                     concentration = 1.0,
                     verbosity = 0)
        
        o.rate_scaling = "Bimolecular_Calibrate"
        o.stop_conditions = [complete_stopcondition] # must be a list of StopCondition items.
        # for testing purposes, see how long it takes to "touch":
        #o.stop_conditions = [touching_stopcondition]

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

class Bimolecular_Calibration( Multistrand_Suite_Base ):
    """ Test Suite for calibrating the bimolecular parameters.
    """

    def __init__(self, param_numbers, length, num_trials, times=None ):
        
        self._suite = unittest.TestSuite()
        if times == None:
            times = [2.0e2] * len(param_numbers) # 'seconds' max time.
            
        for p,time in zip(param_numbers,times):
            pm, file_prefix = args_from_param_number(p)
            if not os.path.isdir( file_prefix ):
                os.mkdir( file_prefix )

            self.add_test_cases(file_prefix, length, time, pm[0], pm[1], pm[2], pm[3], pm[4], num_trials)

    def add_test_cases( self, file_prefix, length_chosen, time_chosen, dangles, ratemethod, temperature, energymodel, substrate_type, num_trials ):
        
        for i in range(num_trials):
            self._suite.addTest( Bimolecular_Calibration_TestCase("{idx}:{time}:{prefix}:{length}:{dangle}:{rate}:{temp}:{em}:{st}".format(idx=i, time=time_chosen, prefix = file_prefix, length=length_chosen, dangle=dangles, rate=ratemethod, temp=float(temperature), em=energymodel, st=substrate_type)))
                


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


def args_from_param_number( num ):
    """ Returns values needed for Bimolecular_Calibration for a given parameter index so that it can generate the test cases."""
    
    parameter_list = list_of_parameters()
    pm = parameter_list[num]
    args = (pm, "calibration_{0[3]}_{0[4]}_{0[1]}_{0[0]}_{1}C/".format( pm, int(pm[2] - 273.15)))
    return args


### Following gets executed if we're running this program directly, e.g. it sets up the test cases and runs them.


if __name__ == '__main__':
    
    calibration = Bimolecular_Calibration([6,7,8,9,10,11], 10, 10) #param type 0, length 10, 1 trial run (1000 trajectories)
    calibration.runTests_Async()
 



### NOTES
###
### Running at conc = .01M, 1 trial at 1e1 max seconds, using param 0:
##    Calibration data set for determining bimolecular parameters.
##    Unimolecular: 514000.0  Bimolecular: 1.0. 
##     Model: [Nupack] Substrate: [DNA] Dangles: [None] Rate Method: [Kawasaki]
##
## it took 125s running time on laptop to get a single trajectory timeout.
## I think perhaps we need a much higher concentration if we want the sim times to
## be reasonable.

## Same result, ~125s running time at conc 1M, no complete trajectories.
