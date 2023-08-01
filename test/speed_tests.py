# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
This test generates random sequences and runs a number of trajectories.
"""

import unittest
import subprocess
import os.path
import timeit

from multiprocessing import Pool, cpu_count

from multistrand.system import SimSystem
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
    """
    This test case handles the comparison between Kinfold and Multistrand on
    random sequences.
    """
    def setUp(self):
        f = open('test_sequences.txt')
        self.sequences = [i.strip('\n') for i in f.readlines()]
        self.indexed_sequences = {}
        for i in self.sequences:
            print(i)
            if hasattr( self.indexed_sequences, str(len(i)) ):
                self.indexed_sequences[len(i)].append(i)
            else:
                self.indexed_sequences[len(i)] = [i]
        f.close()

    def tearDown(self):
        self.sequences = []

    def helper_create_Multistrand_infile(self,sequence,time,count):
        """
        Reuses some code from strandedit to create Multistrand input files.
        """
        infile = ""
        
        infile += '#Strands\n'
        infile += f'Strand1, {sequence}\n'

        #Starting Structure(s)
        infile += '#StartStructure\n'
        infile += 'Strand1\n'
        infile += f"{'.' * len(sequence)}\n"
        #Other options.

        infile += '#StopOptions=0\n' # no stopping conditions.
        infile += '#Energymodel=NUPACK_DNA_2_3\n##\n'
        infile += '#OutputInterval=-1\n'
        infile += '#Temperature=37\n'
        infile += f'#NumSims={count:d}\n'
        infile += f'#SimTime={time:f}\n'
        infile += '#Ratemethod=2\n'  # Kawasaki, like Kinfold default.
        infile += '#dangles=2\n' # Kinfold default

        infile += '#Inter=1.0\n'
        infile += '#Intra=1.0\n'
        infile += '#Logfile=sample.log\n'
        infile += '##'

        return infile
    
    def helper_create_Multistrand_options( self, sequence, time, count):
        """
        Helper
        """
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
        kinfoldproc = subprocess.Popen(
            ["Kinfold","--noShift","--logML","--start","--fpt","--time",
             f"{time:f}","--num","{count:d}","--silent"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE)

        input_str = f"{sequence}\n{'.'*len(sequence)}\n"
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
        input_o = self.helper_create_Multistrand_options( sequence, time, count)

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
            print(f"Executing Kinfold on {sequence} ...")
        time_kinfold = self.helperRunSingle_Kinfold( sequence, time, count )
        if print_flag:
            print(f"{time_kinfold:>70}")
            
        if print_flag:
            print(f"Executing Multistrand [I] on {sequence} ...")
        time_multistrand_I = self.helper_RunSingle_Multistrand_I( sequence, time, count )
        if print_flag:
            print(f"{time_multistrand_I:>70}")

        if print_flag:
            print(f"Executing Multistrand on {sequence} ...")
        time_multistrand = self.helper_RunSingle_Multistrand( sequence, time, count )
        if print_flag:
            print(f"{time_multistrand:>70}")

        return time_kinfold, time_multistrand_I, time_multistrand

    def helper_run_k_trials_single_seq( self, seq, n, k, time, verbose=True ):
        """
        Helper to run k started sims.
        """
        timeresults = []
        count = k
        
        if verbose:
            print(f"Time: [{time} (simulated time units)]\n"
                  f"Trajectories per simulation (n): {n}\n"
                  f"Total simulation calls (k): {k}\n")
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
                print(f"{i+1} complete...")
        if verbose:
            print("\nComplete!")
        return timeresults

    def helper_print_results( self,time, n,k, timeresults, verbose=True, sequence=None ):
        hline = "-" * 70
        fmt_th = "{0:<7}|{1:<12}{2:>3}{3:>14} |{4:<12}{5:>3}{6:>14} |"
        fmt_tr = "{0[0]:<7}|{0[1]:>15e}{0[2]:>14e} |{0[3]:>15e}{0[4]:>14e} |"
        print(hline)
        print(fmt_th.format("Time","Kinfold","[R]","[R]/n*k","Multistrand","[R]","[R]/n*k"))

        tot_kin = 0.0
        tot_mul = 0.0
        for i in timeresults:
            tot_kin += i[1]
            tot_mul += i[3]
            if verbose:
                print(fmt_tr.format(i))
        if verbose:
            print(hline)
            
        print(fmt_tr.format(("Tot:", tot_kin,tot_kin/(time*n*k),tot_mul,tot_mul/(time*n*k))))
        Results.add_results( (sequence, n, k, time, (tot_kin,tot_kin/(time*n*k),tot_mul,tot_mul/(time*n*k)), timeresults ))
        print(hline)

    def test_single_sequence_100_trials_t100( self ):
        """
        Test [time=100.0, n=1, k=100]: Pure run, single seq, 1x time.
        """
        verbose = False
        print("")
        res = self.helper_run_k_trials_single_seq(self.sequences[0] , 1, 100, 100.0, verbose )
        self.helper_print_results( 100.0, 1, 100, res, verbose )

    def test_single_sequence_100_trials_t1000( self ):
        """
        Test [time=1000.0, n=1, k=100]: Pure run, single seq, 10x time.
        """
        verbose = False
        print("")
        res = self.helper_run_k_trials_single_seq(self.sequences[0] , 1, 100, 1000.0, verbose )
        self.helper_print_results( 1000.0, 1, 100, res, verbose )

    def test_single_again( self ):
        self.test_single_sequence_100_trials_t100()

    def test_num_100(self):
        """
        Test [time=1000.0, n=100, k=1]: Single run, single seq, 100
        trajectories.
        """
        verbose = False
        print("")
        res = self.helper_run_k_trials_single_seq(self.sequences[0] , 100, 1, 1000.0, verbose )
        self.helper_print_results( 1000.0, 100, 1, res, verbose )

    def test_num_1000(self):
        """
        Test [time=1000.0, n=10000, k=1]: Single run, single seq, 1000
        trajectories.
        """
        verbose = False
        print("")
        res = self.helper_run_k_trials_single_seq(self.sequences[0] , 1000, 1, 1000.0, verbose )
        self.helper_print_results( 1000.0, 1000, 1, res, verbose )

    def test_begin_random( self ):
        """
        Test [time=1000.0, n=100, k=10]: 10 runs of 100 trajectories of a single
        sequence. Composition test [20].
        """
        verbose = True
        n = 100
        k = 1
        time = 1000.0
        print("")
        for seq in self.sequences[5:105]:
            print(f"Sequence [{len(seq)}]:\n{seq}")
            res = self.helper_run_k_trials_single_seq(seq , n, k, time, verbose )
            self.helper_print_results( time, n, k, res, verbose, sequence=seq )

    def test_begin_random_2( self ):
        """
        Test [time=1000.0, n=100, k=10]: 10 runs of 100 trajectories of a single
        sequence. Composition test [40].
        """
        verbose = True
        n = 100
        k = 1
        time = 1000.0
        print("")
        #        for seq in self.sequences[105:205]:
        for seq in self.sequences[105:110]:
            print(f"Sequence [{len(seq)}]:\n{seq}")
            res = self.helper_run_k_trials_single_seq(seq , n, k, time, verbose )
            self.helper_print_results( time, n, k, res, verbose, sequence=seq )

    def test_length_variance(self):
        """
        Test [time=1000.0, n=100, k=10]: 10 runs of 100 trajectories of a single
        sequence. Length test.
        """
        verbose = False
        n = 100
        k = 10
        time = 1000.0
        print("")
        for seq in self.sequences[0:5]:
            print(f"Sequence [{len(seq)}]:\n{seq}")
            res = self.helper_run_k_trials_single_seq(seq , n, k, time, verbose )
            self.helper_print_results( time, n, k, res, verbose, sequence=seq )

    def test_python_interface(self):
        """
        Test [time=1000.0, n=100, k=1]: 1 run, 100 trajectories, single
        sequence. Python interface vs Multistrand executable vs Kinfold.
        """
        n = 100
        k = 1
        time = 1000.0
        print("")
        seq = self.sequences[0]
        print(f"Sequence [{len(seq)}]:\n{seq}")
        res = self.helper_run_k_trials_single_seq(seq , n, k, time )
        self.helper_print_results( time, n, k, res, sequence=seq )


class SetupSuite( object ):
    """
    Container for basic testing.
    """
    def __init__(self):
        self._suite = unittest.TestSuite()
        self._suite.addTests(
            unittest.TestLoader().loadTestsFromTestCase(
                Speedtest_Random_Sequences ))
        
    def runTests(self):
        if not os.path.isfile('test_sequences.txt'):
            f = open('test_sequences.txt','wt')
            for i in range(5):
                f.write(f"{generate_sequence(20 * (i+1))}\n")
            for i in range(100):
                f.write(f"{generate_sequence(20)}\n")
            for i in range(100):
                f.write(f"{generate_sequence(40)}\n")
            f.close()
            
        if hasattr(self, "_suite") and self._suite is not None:
            unittest.TextTestRunner(verbosity=2).run( self._suite )

    def runTests_Async(self):
        if not os.path.isfile('test_sequences.txt'):
            f = open('test_sequences.txt','wt')
            for i in range(5):
                f.write(f"{generate_sequence(20 * (i+1))}\n")
            for i in range(100):
                f.write(f"{generate_sequence(20)}\n")
            for i in range(100):
                f.write(f"{generate_sequence(40)}\n")
            f.close()
        k = cpu_count()
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
