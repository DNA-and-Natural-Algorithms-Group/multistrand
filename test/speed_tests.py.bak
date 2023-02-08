import unittest
import subprocess
import timeit

import multiprocessing
from multiprocessing import Pool
from multistrand.objects import Strand, Complex
from multistrand.options import Options
from multistrand.utils import generate_sequence

class Results_Store( object ):
    def __init__(self):
        self.results = []
        
    def add_results( self, item ):
        self.results.append(item)

Results = Results_Store()

class Speedtest_Random_Sequences( unittest.TestCase ):
    """ This test case handles the comparison between Kinfold and Multistrand on
    random sequences. """
    def setUp(self):
        f = open('test_sequences.txt')
        self.sequences = [i.strip('\n') for i in f.readlines()]
        self.indexed_sequences = {}
        for i in self.sequences:
            print i
            if hasattr( self.indexed_sequences, str(len(i)) ):
                self.indexed_sequences[len(i)].append(i)
            else:
                self.indexed_sequences[len(i)] = [i]
        f.close()

    def tearDown(self):
        self.sequences = []

    def helper_create_Multistrand_infile(self,sequence,time,count):
        """ reuses some code from strandedit to create Multistrand input files. """
        infile = ""
        
        infile += '#Strands\n'
        infile += 'Strand1, {0}\n'.format( sequence )

        #Starting Structure(s)
        infile += '#StartStructure\n'
        infile += 'Strand1\n'
        infile += '{0}\n'.format( "." * len(sequence))
        #Other options.

        infile += '#StopOptions=0\n' # no stopping conditions.
        infile += '#Energymodel=NUPACK_DNA_2_3\n##\n'
        infile += '#OutputInterval=-1\n'
        infile += '#Temperature=37\n'
        infile += '#NumSims={0:d}\n'.format( count )
        infile += '#SimTime={0:f}\n'.format( time )
        infile += '#Ratemethod=2\n'  # Kawasaki, like Kinfold default.
        infile += '#dangles=2\n' # Kinfold default

        infile += '#Inter=1.0\n'
        infile += '#Intra=1.0\n'
        infile += '#Logfile=sample.log\n'
        infile += '##'

        return infile
    
    def helper_create_Multistrand_options( self, sequence, time, count):
        """ helper """
        
        s1 = Strand("test_strand1",  str(sequence), None )
        c1 = Complex( 1, "start", [s1], "." * len( sequence ) )
        
        o = Options()
        o.simulation_mode = Constants.SIMULATION_MODE['Python Module']
        o.use_stop_states = False
        o.num_simulations = count
        o.simulation_time = time
        o.start_state = [c1]
        o.dangles = Constants.DANGLES['All']
        o.rate_method = Constants.RATEMETHOD['Kawasaki']
        o.bimolecular_scaling = 1.0
        o.unimolecular_scaling = 1.0
        o.temperature = 37.0
        o.boltzmann_sample = False

        return o
    
    def helperRunSingle_Kinfold( self, sequence, time, count=1 ):
        kinfoldproc = subprocess.Popen(["Kinfold","--noShift","--logML","--start","--fpt","--time","{0:f}".format(time),"--num","{0:d}".format(count),"--silent"], stdin=subprocess.PIPE, stdout=subprocess.PIPE )

        input_str = "{0}\n{1}\n".format( sequence, "."*len(sequence) )
        def runOnce():
            self.output_kinfold, _ = kinfoldproc.communicate( input_str )
        t=timeit.Timer( runOnce )

        result_time = t.timeit(1) # only run it once.
        return result_time


    def helper_RunSingle_Multistrand( self, sequence, time, count=1 ):
        input_str = self.helper_create_Multistrand_infile( sequence, time, count)
        multistrandproc = subprocess.Popen(["Multistrand"],stdin=subprocess.PIPE, stdout=subprocess.PIPE )
        def runOnce():
            self.output_multistrand, _ = multistrandproc.communicate( input_str )
        t=timeit.Timer(runOnce)
        result_time = t.timeit(1)
        return result_time

    def helper_RunSingle_Multistrand_I( self, sequence, time, count=1 ):
        # import pdb
        # pdb.set_trace()
        input_o = self.helper_create_Multistrand_options( sequence, time, count)
        from multistrand.system import SimSystem

        def runOnce():
            s = SimSystem( input_o )
            s.start()
            self.output_multistrand = str(input_o.interface.results)
            del s
            
        t=timeit.Timer(runOnce)
        result_time = t.timeit(1)
        #print(self.output_multistrand)
        return result_time
    
    def helper_run_single( self, sequence, time, count = 1, print_flag = False ):
        if print_flag:
            print("Executing Kinfold on {0} ...".format( sequence ))
        time_kinfold = self.helperRunSingle_Kinfold( sequence, time, count )
        if print_flag:
            print("{0:>70}".format( time_kinfold ))
            
        if print_flag:
            print("Executing Multistrand [I] on {0} ...".format( sequence ))
        time_multistrand_I = self.helper_RunSingle_Multistrand_I( sequence, time, count )
        if print_flag:
            print("{0:>70}".format( time_multistrand_I ))

        if print_flag:
            print("Executing Multistrand on {0} ...".format( sequence ))
        time_multistrand = self.helper_RunSingle_Multistrand( sequence, time, count )
        if print_flag:
            print("{0:>70}".format( time_multistrand ))

        return time_kinfold, time_multistrand_I, time_multistrand

    def helper_run_k_trials_single_seq( self, seq, n, k, time, verbose=True ):
        """ Helper to run k started sims. """
        timeresults = []
        count = k
        
        if verbose:
            print("Time: [{0} (simulated time units)]\nTrajectories per simulation (n): {1}\nTotal simulation calls (k): {2}\n".format( time, n, k) )
            if k > 1:
                print("First Trial:")
            times = self.helper_run_single( seq, time, n,  print_flag=True )
            timeresults.append( (time, times[0],times[0] / time,  times[1], times[1]/time ))
            if k > 1:
                print("Running Remaining...")
            count -= 1
            
        for i in range(count):
            times = self.helper_run_single( seq, time, n, print_flag=verbose )
            timeresults.append( (time, times[0],times[0] / time,  times[1], times[1]/time ))
            if verbose and (i+1) % 25 == 0:   # modulus binds weakly
                print("{0} complete...".format(i+1))
        if verbose:
            print("\nComplete!")
        return timeresults

    def helper_print_results( self,time, n,k, timeresults, verbose=True, sequence=None ):
        #print("Real Time per Simulated Time Unit:")
        print("-"*70)
        print("{0:<7}|{1:<12}{2:>3}{3:>14} |{4:<12}{5:>3}{6:>14} |".format("Time","Kinfold","[R]","[R]/n*k","Multistrand","[R]","[R]/n*k"))
        tot_kin = 0.0
        tot_mul = 0.0
        for i in timeresults:
            tot_kin += i[1]
            tot_mul += i[3]
            if verbose:
                print("{0[0]:<7}|{0[1]:>15e}{0[2]:>14e} |{0[3]:>15e}{0[4]:>14e} |".format( i ))
        if verbose:
            print("-"*70)
            
        print("{1:<7}|{0[0]:>15e}{0[1]:>14e} |{0[2]:>15e}{0[3]:>14e} |".format( (tot_kin,tot_kin/(time*n*k),tot_mul,tot_mul/(time*n*k)) , "Tot:"))
        Results.add_results( (sequence, n, k, time, (tot_kin,tot_kin/(time*n*k),tot_mul,tot_mul/(time*n*k)), timeresults ))
        print("-"*70)

    def test_single_sequence_100_trials_t100( self ):
        """ Test [time=100.0, n=1, k=100]: Pure run, single seq, 1x time. """
        return
        verbose = False
        print("")
        res = self.helper_run_k_trials_single_seq(self.sequences[0] , 1, 100, 100.0, verbose )
        self.helper_print_results( 100.0, 1, 100, res, verbose )

    def test_single_sequence_100_trials_t1000( self ):
        """ Test [time=1000.0, n=1, k=100]: Pure run, single seq, 10x time. """
        return
        verbose = False
        print("")
        res = self.helper_run_k_trials_single_seq(self.sequences[0] , 1, 100, 1000.0, verbose )
        self.helper_print_results( 1000.0, 1, 100, res, verbose )

    def test_single_again( self ):
        return
        self.test_single_sequence_100_trials_t100()

    def test_num_100(self):
        """ Test [time=1000.0, n=100, k=1]: Single run, single seq, 100 trajectories. """
        return
        verbose = False
        print("")
        res = self.helper_run_k_trials_single_seq(self.sequences[0] , 100, 1, 1000.0, verbose )
        self.helper_print_results( 1000.0, 100, 1, res, verbose )

    def test_num_1000(self):
        """ Test [time=1000.0, n=10000, k=1]: Single run, single seq, 1000 trajectories. """
        return
        verbose = False
        print("")
        res = self.helper_run_k_trials_single_seq(self.sequences[0] , 1000, 1, 1000.0, verbose )
        self.helper_print_results( 1000.0, 1000, 1, res, verbose )

    def test_begin_random( self ):
        """ Test [time=1000.0, n=100, k=10]: 10 runs of 100 trajectories of a single sequence. Composition test [20]."""
        return
        verbose = True
        n = 100
        k = 1
        time = 1000.0
        print("")
        for seq in self.sequences[5:105]:
            print("Sequence [{0}]:\n{1}".format( len(seq), seq ))
            res = self.helper_run_k_trials_single_seq(seq , n, k, time, verbose )
            self.helper_print_results( time, n, k, res, verbose, sequence=seq )

    def test_begin_random_2( self ):
        """ Test [time=1000.0, n=100, k=10]: 10 runs of 100 trajectories of a single sequence. Composition test [40]."""
        verbose = True
        n = 100
        k = 1
        time = 1000.0
        print("")
        #        for seq in self.sequences[105:205]:
        for seq in self.sequences[105:110]:
            print("Sequence [{0}]:\n{1}".format( len(seq), seq ))
            res = self.helper_run_k_trials_single_seq(seq , n, k, time, verbose )
            self.helper_print_results( time, n, k, res, verbose, sequence=seq )

        
    def test_length_variance(self):
        """ Test [time=1000.0, n=100, k=10]: 10 runs of 100 trajectories of a single sequence. Length test."""
        return
        verbose = False
        n = 100
        k = 10
        time = 1000.0
        print("")
        for seq in self.sequences[0:5]:
            print("Sequence [{0}]:\n{1}".format( len(seq), seq ))
            res = self.helper_run_k_trials_single_seq(seq , n, k, time, verbose )
            self.helper_print_results( time, n, k, res, verbose, sequence=seq )

    def test_python_interface(self):
        """ Test [time=1000.0, n=100, k=1]: 1 run, 100 trajectories, single sequence. Python interface vs Multistrand executable vs Kinfold. """
        return
        n = 100
        k = 1
        time = 1000.0
        print("")
        seq = self.sequences[0]
        print("Sequence [{0}]:\n{1}".format( len(seq), seq ))
        res = self.helper_run_k_trials_single_seq(seq , n, k, time )
        self.helper_print_results( time, n, k, res, sequence=seq )
    


class SetupSuite( object ):
    """ Container for basic testing. """

    def __init__(self):
        self._suite = unittest.TestSuite()
        self._suite.addTests(
            unittest.TestLoader().loadTestsFromTestCase(
                Speedtest_Random_Sequences ))
        
    def runTests(self,print_results=False):
        import os.path
        if not os.path.isfile('test_sequences.txt'):
            f = open('test_sequences.txt','wt')
            for i in range(5):
                f.write( "{0}\n".format(utils.generate_sequence( 20 * (i+1) ) ))
            for i in range(100):
                f.write( "{0}\n".format(utils.generate_sequence( 20 ) ))
            for i in range(100):
                f.write( "{0}\n".format(utils.generate_sequence( 40 ) ))
            f.close()
            
        if hasattr(self, "_suite") and self._suite is not None:
            unittest.TextTestRunner(verbosity=2).run( self._suite )
        # if print_results:
        #     self.printResults()

    def runTests_Async(self):
        import os.path
        if not os.path.isfile('test_sequences.txt'):
            f = open('test_sequences.txt','wt')
            for i in range(5):
                f.write( "{0}\n".format(utils.generate_sequence( 20 * (i+1) ) ))
            for i in range(100):
                f.write( "{0}\n".format(utils.generate_sequence( 20 ) ))
            for i in range(100):
                f.write( "{0}\n".format(utils.generate_sequence( 40 ) ))
            f.close()
        k = multiprocessing.cpu_count()
        p = Pool( processes = k )
        p.map( MyRunner , iter(self._suite), chunksize = 1 )
        p.close()
        p.join()

class MyRunner( object ):
    def __init__( self, testcase ):
        unittest.TextTestRunner( verbosity=2).run( testcase )
        
if __name__ == '__main__':
    suite = SetupSuite()
    #suite.runTests()
    suite.runTests_Async()

