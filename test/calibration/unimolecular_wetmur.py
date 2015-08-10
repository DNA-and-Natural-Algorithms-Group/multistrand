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
import unittest
import warnings

import os, os.path
import cPickle

import random

#from scipy.constants import R
import numpy as np
R_gas_constant = 1.9858775e-3  #units of kcal/K/mol


# for IPython, some of the IPython libs used by unittest have a
# deprecated usage of BaseException, so we turn that specific warning
# off.
warnings.filterwarnings("ignore", r"BaseException[.]message has been deprecated as of Python 2[.]6", DeprecationWarning)

# import multiprocessing
# from multiprocessing import Pool

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
    from multistrand.system import initialize_energy_model, energy
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


def set_energy_options( dangles, temperature ):
    o = Options( substratetype = "DNA",
                 parameter_type = "Nupack",
                 dangles = dangles,
                 rate_method = "Metropolis",
                 num_sims = 1,
                 sim_time = 100.0,
                 temperature = temperature,
                 concentration = 1.0e-3,
                 verbosity = 0 )
    o.rate_scaling = "Bimolecular_Calibrate"

    initialize_energy_model( o )
    return o

def rate_Kawasaki( start_energy, end_energy, temperature_in_K ):
    dG = end_energy - start_energy
    rate = np.exp( -dG / (2 * R_gas_constant * temperature_in_K ) )
    return rate

def rate_Metropolis( start_energy, end_energy, temperature_in_K ):
    dG = end_energy - start_energy
    if dG <= 0:
        return 1.0
    rate = np.exp( -dG / (R_gas_constant * temperature_in_K ))
    return rate

def rate_zippering( temperature_in_K ):
    ## following Wetmur 1967, eqn 40,41, pg 367+368.
    # eqn 40:   kf = A exp( -E_act / RT )
    # eqn 41:   kf = 3e9 s^-1   [at L=N=1, T=353K, E_act=7.5 kcal/mol]

    A = 3e9 / np.exp( -7.5 / (R_gas_constant * 353.0))

    zippering_rate = A * np.exp( -7.5 / (R_gas_constant * temperature_in_K))

    return zippering_rate,A

def calculate_rates( sequence, options, secondary_sequence = None ):
    ### sequence should always be length 3.

    ### One might ask whether we should add a G-C pair to
    ### the sequences, to eliminate any weirdness from extra A-T
    ### penalties. However, we should expect this to not change the
    ### rates in any way, as the top strand 5' side A-T pair (if present)
    ### has an AT penalty in the /open loop/, thus it gets cancelled back out
    ### when we subtract energies.

    short_seq = Domain( name="test_region", sequence=sequence )

    strand_top = Strand(name="test_strand",domains=[short_seq])

    if secondary_sequence is not None:
        short_seq_2 = Domain( name="test_region2",sequence=secondary_sequence)
        strand_bottom = Strand(name="test_strand2",domains=[short_seq_2])
    else:
        strand_bottom = strand_top.C

    complex_start = Complex( strands=[strand_top,strand_bottom], structure="(..+..)")
    complex_end = Complex( strands=[strand_top,strand_bottom], structure="((.+.))")
    
    
    energy_start = energy( [complex_start] )
    energy_end   = energy( [complex_end]   )

    rate = rate_Kawasaki( energy_start[0], energy_end[0], options.temperature )
    rate_M = rate_Metropolis( energy_start[0], energy_end[0], options.temperature )
    return rate, rate_M


def build_sequences_non_complementary_dangles():
    seqs_nc_primary = []
    seqs_nc_secondary = []
    for a in ["G","C","A","T"]:
        for b in ["G","C","A","T"]:
            for c in ["G","C","A","T"]:
                for d in ["G","C","A","T"]:
                    seqs_nc_primary.append( a + b + c )
                    seqs_nc_secondary.append( d + comp(b) + comp(a) )
    return seqs_nc_primary, seqs_nc_secondary

def build_sequences():
    seqs = []
    for a in ["G","C","A","T"]:
        for b in ["G","C","A","T"]:
            for c in ["G","C","A","T"]:
                seqs.append( a + b + c )
    return seqs


def build_rate_list( temperature ):
    results = {}
    seqs = build_sequences()
    for dangles in ["None","Some","All"]:
        
        o = set_energy_options( dangles, temperature )
        rates_kawa, rates_metro = zip( * [calculate_rates(i, o) for i in seqs])
        results[dangles] = {"Kawasaki":rates_kawa,"Metropolis":rates_metro}

    return results

def final_results( temperature):
    zip_rate, _ = rate_zippering( temperature )
    print("Zippering rate at {0:3.2f}K: {1:.2e}".format( temperature, zip_rate ))
    
    final_rate_dict = build_rate_list( temperature )
    for dangles in ["None","Some","All"]:
        print("Dangles = {0}".format(dangles))
        kawa_rates = final_rate_dict[dangles]["Kawasaki"]
        metro_rates = final_rate_dict[dangles]["Metropolis"]

        kawa_mean = np.mean( kawa_rates )
        kawa_min, kawa_max = np.min(kawa_rates) , np.max( kawa_rates )

        metro_mean = np.mean( metro_rates )
        metro_min, metro_max = np.min(metro_rates) , np.max( metro_rates )

        print("  Kawasaki (mean/min/max)           : {0:.2e}/{1:.2e}/{2:.2e}".format( kawa_mean, kawa_min, kawa_max ))
        print("  Kawasaki (k_uni/maxrate/minrate)  : {0:.1e}/{1:.2e}/{2:.2e}".format( zip_rate / kawa_mean, (zip_rate / kawa_mean) * (kawa_max / kawa_mean), (zip_rate / kawa_mean) * (kawa_min / kawa_mean) ))

        print("  Metropolis (mean/min/max)         : {0:.2e}/{1:.2e}/{2:.2e}".format( metro_mean, metro_min, metro_max ))
        print("  Metropolis (k_uni/maxrate/minrate): {0:.1e}/{1:.2e}/{2:.2e}".format( zip_rate / metro_mean,(zip_rate / metro_mean) * (metro_max / metro_mean), (zip_rate / metro_mean) * (metro_min / metro_mean) ))


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

# if __name__ == '__main__':
#     if len(sys.argv) < 2:
#         idx = 0
#     else:
#         idx = int(sys.argv[1])

#     calibration = Threeway_BM_Calibration( *args_from_param_number( idx ))
#     calibration.runTests_Async()

#    calibration = Threeway_BM_Calibration( [10,20,30], [1.0e5,1.0e6,1.0e6],  "None", "Kawasaki", 310.15, "Nupack", "DNA", 'calibration_test/' )


