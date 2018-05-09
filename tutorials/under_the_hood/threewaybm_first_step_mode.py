# threewaybm_first_step_mode.py
#
# adapted from three_way_sims_example.py by Niranjan Srinivas and Joseph Schaeffer.
#
# This is the type of code used in the case studies... time to take a look at the innards...
# The main concept being shown here is the process for doing simulations in parallel on multicore machines.
# If you want a total number T of trajectories to be simulated, you can divide it into J jobs of S simulations each
# (with T = J*S) and then spawn off J instances of Multistrand, each of which creates and runs one SimSystem 
# that performs S simulation trials. The J jobs might not all be the same setup.  We demonstrate a simple batch scheduler
# that takes a list of jobs, counts the number of cores on your machine, and runs that many in parallel -- starting a new
# process whenever an old one completes. The results of each run are stored in data files in appropriately named
# directories "Data_toehold*" which can later be all loaded back and recombined into a single data set for each simulated
# system.  The analysis to extract rates from the simulations is similar to that in hybridization_first_step_mode.py,
# except that here we don't use Boltzmann sampling and we don't calculate error bars, which simplified things.
#
#
# Usage for a quick-running test using just toehold length 3 and a just a few trials, taking about 5 minutes:
# python -i threewaybm_first_step_mode.py test
#
# Usage for a full simulation with Metropolis kinetics and toehold lengths 1-4, about a 3 hours:
# python -i threewaybm_first_step_mode.py short
#
# Usage for a full simulation with Metropolis kinetics and toehold lengths 5-9, almost an hour:
# python -i threewaybm_first_step_mode.py medium
#
# Usage for a full simulation with Kawasaki kinetics and toehold 10-15, under half an hour:
# python -i threewaybm_first_step_mode.py long
#
# After the simulations are done, and you are dumped into interactive python, you can analyze the results:
# >>> final_results(3,"Metropolis")
#
# If you have already run the simulations, you can do the analysis for all toehold lengths:
# python -i threewaybm_first_step_mode.py analyze
#
# Times are for a 4-core Intel Core i7 MacBook Pro, 2.6 GHz, year 2014 model.

# Import things you need
import sys, os

import cPickle
import random
import numpy as np

if False:  # only needed if you're having trouble with your Multistrand installation
    import multistrand_setup

try:
    from multistrand.objects import *
    from multistrand.options import Options
    from multistrand.system import SimSystem

except ImportError:
    print("Could not import Multistrand.")
    raise

try:
    import multiprocessing
    from multiprocessing import Pool
    
except ImportError:
    print("Could not import multiprocessing.")
    raise

#########

# for StopCondition and Macrostate definitions:
Exact_Macrostate = 0  # match a secondary structure exactly (i.e. any system state that has a complex with this exact structure)
Bound_Macrostate = 1  # match any system state in which the given strand is bound to another strand
Dissoc_Macrostate = 2  # match any system state in which there exists a complex with exactly the given strands, in that order
Loose_Macrostate = 3  # match a secondary structure with "don't care"s, allowing a certain number of disagreements
Count_Macrostate = 4  # match a secondary structure, allowing a certain number of disagreements
# see Schaeffer's PhD thesis, chapter 7.2, for more information


def create_setup(toehold_length, num_traj, rate_method_k_or_m):
    # essentially, creates the options object and prepares to simulate
    
    toehold_seq = "TCTCCATGTCACTTC"  # toehold sequence

    bm_design = "CCCTCATTCAATACCCTACG"  # branch migration domain

    # creating Domain objects
    toehold = Domain(name="toehold", sequence=toehold_seq[0:toehold_length])
    branch_migration = Domain(name="bm", sequence=bm_design)
     
    incoming = branch_migration + toehold  # you can concatenate domains to make domains
    incoming.name = "incoming"
    # here the "invader" (according to Zhange & Winfree 2009) is called the "incoming" strand.
    
    toehold_extra = Domain(name="toehold_extra", sequence=toehold_seq[toehold_length:]) 
    substrate = toehold_extra.C + toehold.C + branch_migration.C
    
    incumbent = Strand(name="incumbent", domains=[branch_migration])

    # creates the substrate-incumbent start complex  
    start_complex_substrate_incumbent = Complex(strands=[substrate, incumbent], structure="..(+)")

    # creates incoming ("invader") complex. 
    start_complex_incoming = Complex(strands=[incoming], structure="..") 
    # by default, Boltzmann sampling is *NOT* used in first step mode.  See hybridization_first_step_mode.py for Boltzmann sampling.

    # creates a complex for a "succcessful displacement" stop condition. This is the incumbent strand forming a complex of its own which means it has been displaced.
    complete_complex_success = Complex(strands=[incumbent], structure=".")
    success_stop_condition = StopCondition("SUCCESS", [(complete_complex_success, Dissoc_Macrostate, 0)])
    # creates the successful displacement stop condition

    # complex to create failed displacement stop condition; incumbent falls off.   
    failed_complex = Complex(strands=[incoming], structure="..")  
    failed_stop_condition = StopCondition("FAILURE", [(failed_complex, Dissoc_Macrostate, 0)]) 
    # creates failed stop condition
    
    o = Options(simulation_mode="First Step", parameter_type="Nupack", substrate_type="DNA",
                rate_method=rate_method_k_or_m, num_simulations=num_traj, simulation_time=10.0,  # note the 10 second simulation time, to make sure simulations finish
                dangles="Some", temperature=25 + 273.15, rate_scaling="Calibrated", verbosity=0)

    o.start_state = [start_complex_incoming, start_complex_substrate_incumbent]
    o.stop_conditions = [success_stop_condition, failed_stop_condition]
    return o

# The following is somewhat general code for running a bunch of Multistrand simulations on as many cores as this computer has available.
# Directories "Data_toehold_*" are created and simulation results are stored there.
# When all the simulations are done, the data is all read back in, and processed together.
#
# Note the following important restriction:  all the simulations must use "the same energy model", that is, they must all request
# the same temperature, join_concentration, substrate_type (RNA vs DNA), dangles option, parameter_type (NUPACK vs Vienna),
# rate method (Kawasaki or Metropolis), and rate parameters (e.g. Calibrated, Unitary, etc).
# This is because all multithreaded parallel simulations will use the same energy/kinetics model.  The first one to run will initiallize
# the energy model, and all the rest will use it.  If the simulation script (e.g. actual_simunlation() below) were to call
# system.initialize_energy_model(o) -- either with the same conditions or different onese -- then the likely result would be
# a segmentation fault.  Don't do it. 


class Multistrand_Suite_Base(object):
    """ Base class for test suites - defines async run, etc. """

    def runTests_Async(self, shuffle_tasks=True):
        """ This runner runs the tests via an asynchronous pool of
        workers.

        The number of workers is equal to the CPU core count, and this
        function only returns once every task has been completed. It's
        possible we should up the chunksize so that tasks are
        partitioned in bigger groups, but since all the cores being
        used are local, it seems like overkill."""
        starttime = os.times()
        k = multiprocessing.cpu_count()
        p = Pool(processes=k)

        if shuffle_tasks:
            random.shuffle(self.jobs)
        
        p.map(MyRunner , iter(self.jobs), chunksize=1)
        p.close()
        p.join()
        endtime = os.times()
        print("Async run complete! Processing took [{0[4]}] seconds of real time before completion. [u/s/cu/cs]:[{0[0]}/{0[1]}/{0[2]}/{0[3]}]".format([j - i for i, j in zip(starttime, endtime)]))
    

class MyRunner(object):

    def __init__(self, parameter_set):
        actual_simulation(*parameter_set)


def actual_simulation(toehold_length, num_traj, rate_method_k_or_m, index):
    print "Starting %d simulations for toehold length %d and %s kinetics." % (num_traj, toehold_length, rate_method_k_or_m)
    o = create_setup(toehold_length, num_traj, rate_method_k_or_m)
    s = SimSystem(o)
    s.start()
    prefix = "Data_toehold_{0}".format(toehold_length)
    filename = "DATA_{0}_{1}_{2:04}.dat".format(rate_method_k_or_m, toehold_length, index)
    full_filename = os.path.join(prefix, filename)
    f = open(full_filename, 'wb')
    cPickle.dump(o.interface.results, f, protocol=-1)
    f.close()


class ThreeWaySimRunner(Multistrand_Suite_Base):

    def __init__(self, list_toeholds, list_counts, rate_method_k_or_m, jobcount=32, offset_for_unique_name=0):
        self.jobs = []

        def jobsfromtoelength(toelength, toecount):
            return [(toelength, toecount, rate_method_k_or_m, i + offset_for_unique_name) for i in range(jobcount)]

        for i, j in zip(list_toeholds, list_counts):
            self.jobs = self.jobs + (jobsfromtoelength(i, j))
            if not os.path.isdir("Data_toehold_{0}".format(i)):
                os.mkdir("Data_toehold_{0}".format(i))

########## Stuff above deals with running the simulations.
########## Stuff below deals with analyzing the simulations.


def load_file(dirname, filename):
    fullname = os.path.join(dirname, filename)

    f = open(fullname, 'rb')
    res = cPickle.load(f)
    f.close()

    return res


def process_dataset(dataset):
    # essentially identifies trajectories of each kind and takes the first passage times
    
    forward = [i for i in dataset if i.tag == "SUCCESS"]
    forward_array = np.zeros(len(forward))
    forward_array[:] = [i.time for i in forward]

    # if after 10 simulated seconds, the simulation didn't reach either stop state, its tag is "None", and it is lumped with failures.
    reverse = [i for i in dataset if i.tag == "FAILURE" or i.tag == Options.STR_NOINITIAL or i.tag == Options.STR_TIMEOUT]  
    reverse_array = np.zeros(len(reverse))
    reverse_array[:] = [i.time for i in reverse]

    collision_array = np.zeros(len(dataset))
    collision_array[:] = [i.collision_rate for i in dataset]
    
    return (forward_array, reverse_array, collision_array)


def load_dataset(toehold_length, rate_method):
    # Looks in the standard directory for all matching filenames,
    # loads each file and condenses the results.
    # If there are N files, the result_list will have N non-None items.
    # Each item is a tuple of arrays: (completion times for all successful trials, completion times for all failed trials, initial collition rates for all trials)
    # Thus these three arrays have lengths (S, F, S+F).
    
    dirname = "Data_toehold_{0}".format(toehold_length)
    if not os.path.isdir(dirname):
        print("Error: directory {0} does not exist.\n".format(dirname))
        return None

    files = os.listdir(dirname)
    sample_files = [i for i in files if i.startswith('DATA_' + rate_method)]
    sample_files.sort()
    
    def number_from_fname(filename):
        return filename.rstrip('.dat').split('_')[-1]

    result_list = [None] * (int(number_from_fname(sample_files[-1])) + 1)
    for f in sample_files:
        print("Processing file: {0}".format(f))
        result_list[ int(number_from_fname(f)) ] = process_dataset(load_file(dirname, f))

    return result_list


def combine_dataset(sampleset):
    if sampleset == None :
        print "No data, no rates computed."
        return None

    f, r, c = zip(*[i for i in sampleset if i is not None])
    forward_times = np.concatenate(f, axis=0)
    reverse_times = np.concatenate(r, axis=0)
    # these are essentially "forward" (successful displacement) and "reverse" (failed displacement) first passage times.
    
    collision_rates = np.concatenate(c, axis=0)
    
    print "In all: %d collisions, %d successes (forward), %d failures (reverse), %d unfinished" % \
        (len(collision_rates), len(forward_times), len(reverse_times), len(collision_rates) - len(forward_times) - len(reverse_times))

    # calculations from Joseph's PhD thesis.
    # Note: we are not using Boltzmann sampling, so the formulas below are fine.  But with Boltzmann sampling, we would need to partition the kcollision values.
     
    dTsuccess_uni = np.mean(forward_times)
    k2 = 1.0 / dTsuccess_uni
    dTfail_uni = np.mean(reverse_times)
    k2prime = 1.0 / dTfail_uni
    kcollision = np.mean(collision_rates)
    N_forward = len(forward_times)
    N_fail = len(reverse_times)
    k1 = (float(N_forward) / (N_forward + N_fail)) * kcollision
    k1prime = (float(N_fail) / (N_forward + N_fail)) * kcollision

    z = 50e-9  # concentration 50 nM used for k_eff calculation
    
    dTcoll = 1 / ((k1 + k1prime) * z)
    dTfail = dTcoll + dTfail_uni
    dTforward = dTcoll + dTsuccess_uni
    dTcorrect = dTfail * k1prime / k1 + dTforward
    
    keff = (1 / dTcorrect) * (1 / z)

    print "   k1 = %g /M/s, k2 = %g /s, k1' = %g /M/s, k2' = %g /s, and k_eff = %g /M/s at 50 nM" % \
        (k1, k2, k1prime, k2prime, keff)
    print "   incoming + substrate -->{k1}    succesful_intermediate -->{k2}  waste + incumbent"
    print "   incoming + substrate -->{k1'} unsuccesful_intermediate -->{k2'} incoming + substrate"

    return [N_forward, N_fail, keff, k1, k1prime, k2, k2prime, kcollision]    


# How you do the analysis
def final_results(toehold_length, rate_method):
    # given a toehold length and a rate method, this function returns [N_forward, N_fail, keff, k1, k1prime, k2, k2prime, kcollision]    
    print
    print "Results for toehold length %d (with %s kinetics)" % (toehold_length, rate_method)
    return combine_dataset(load_dataset(toehold_length, rate_method))

###### If run from the command line, the simulations are run but not analyzed.


if __name__ == '__main__':
    if len(sys.argv) < 2:
        idx = 0
        print "Must give argument 'test', 'short', 'medium', 'long', or 'analyze'."
    else:
        idx = str(sys.argv[1])
     
    if idx == 'test':
        testrun = ThreeWaySimRunner([3], [50], "Metropolis", jobcount=3, offset_for_unique_name=399)
        # testrun = ThreeWaySimRunner([list of toehold lengths],
        #                            [list of #sims per job, for each toehold length], 
        #                            Rate method ("Kawasaki" or "Metropolis"), 
        #                            jobcount=#jobs,
        #                            offset_for_unique_name=index you should use for naming data file so you don't over-write previous saved data files)
        testrun.runTests_Async()
        
    if idx == 'short':  # 20K total trajectories each
        testrun = ThreeWaySimRunner([1, 2, 3, 4], [1000, 1000, 1000, 1000], "Metropolis", jobcount=20, offset_for_unique_name=0)
        testrun.runTests_Async()

    if idx == 'medium':  # 4000 total trajectories each
        testrun = ThreeWaySimRunner([5, 6, 7, 8, 9], [1000, 1000, 1000, 1000, 1000], "Metropolis", jobcount=4, offset_for_unique_name=0)
        testrun.runTests_Async()
        
    if idx == 'long':  # 400 total trajectories each
        testrun = ThreeWaySimRunner([10, 11, 12, 13, 14, 15], [100, 100, 100, 100, 100, 100], "Metropolis", jobcount=4, offset_for_unique_name=0)
        testrun.runTests_Async()    

    if idx == 'analyze':
        toelengths = range(1, 16)
        k1s = range(1, 16)
        for n in toelengths:
            rates = final_results(n, "Metropolis")
            if rates == None:
                k1s[n - 1] = None
            else:
                k1s[n - 1] = rates[3]
        
        import matplotlib
        import matplotlib.pylab as plt
        
        plt.figure(1)
        plt.semilogy(toelengths, k1s, 'ko')
        plt.hold(True)
        plt.semilogy(toelengths, k1s, 'k-')
        plt.hold(False)
        plt.title("Toehold-mediated three-way strand displacement")
        plt.xlabel("Toehold Length (nt)", fontsize='larger')
        plt.ylabel("Bimolecular rate constant k1 (/M/s)", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.xticks(fontsize='larger')
        plt.show()

        print "As noted in Srinivas et al NAR 2013, Multistrand 2.0 can qualitatively but not quantitatively reproduce Zhang & Winfree JACS 2009."
        print "Please wait for Multistrand 3.0 for better quantitative results."
