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

import sys, os, os.path
import cPickle

import numpy as np

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
    # if os.getcwd[0] != '/':
    #     rootpath = os.path.expanduser('~' + os.getcwd()[:idx])
    rootpath = os.path.abspath(os.getcwd()[:idx])
    #abspath cleans up the tail of the pathname as needed.
    if rootpath not in sys.path:
        sys.path.append( rootpath )

if multihome != None:
    sys.path.append(multihome)

try:
    from multistrand.options import Options, Constants


except:
    # we want to tell the user how to fix this, but then reraise so it still fails. 
    print("Could not import Multistrand, please add it to your sys.path, or make sure that MULTISTRANDHOME is set correctly. This sub-program can also be run from the native test/calibration/ directory.")
    raise

# WARNING: This is copied from unimolecular.py, any changes should go to both files
def list_of_parameters():

    dangles_types = ["None", "Some", "All"]
    rate_methods = ["Kawasaki", "Metropolis"]
    temperatures = [310.15, 273.15+25]
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


def filename_from_param_number( num ):
    """ Returns the proper argument list for Threeway_BM_Calibration, using the index num into the list of all parameter combinations """
    if num == -1:
        return 'calibration_test'
    
    parameter_list = list_of_parameters()
    pm = parameter_list[num]
    return "calibration_{0[3]}_{0[4]}_{0[1]}_{0[0]}_{1}C/".format( pm, int(pm[2] - 273.15))


def load_file( dirname, filename):
    fullname = os.path.join( dirname, filename )

    f = open( fullname, 'rb')
    res = cPickle.load(f)
    cal_string = cPickle.load(f)
    #note: could read off the options object here as well
    f.close()

    return res,cal_string

def load_dataset( num, length, prefix=None, verbose = True ):
    if prefix is None:
        dirname = filename_from_param_number( num )
    else:
        dirname = prefix
    if not os.path.isdir( dirname ):
        print("Error: directory {0} does not exist.\n".format( dirname ))
        return None

    files = os.listdir( dirname )
    sample_files = [i for i in files if i.endswith('.result.dat') and i.startswith('bimolecular_len_{0}'.format(length))]
    sample_files.sort()
    if len(sample_files) == 0:
        print("No sample files to load for parameter set [{0}]".format(num))
        return []
    
    def number_from_fname( filename ):
        return filename.rstrip('.result.dat').split('_')[-1]

    max_index = max([number_from_fname(fname) for fname in sample_files])
    result_list = [None] * (int(max_index)+1)
    for f in sample_files:
        if verbose:
            print("Processing file: {0}".format( f ))
        data, cal_string = load_file( dirname, f)
        result_list[ int( number_from_fname(f)) ] = process_dataset(data,cal_string)
        

    return result_list

def process_dataset( dataset, cal_string ):
    forward = [i for i in dataset if i.tag == "complete"]
    forward_array = np.zeros( len(forward))
    forward_array[:] = [i.time for i in forward]

    std = np.std(forward_array,ddof = 1)
    err = std / np.sqrt( len( forward) )
    return (np.mean(forward_array), std, err, len(forward), len(dataset), forward_array, cal_string)


def collate_result_list( dataset ):
    statistics = [dataset[i][:-2] for i in range(len(dataset)) if dataset[i] is not None]
    data = [dataset[i][5] for i in range(len(dataset)) if dataset[i] is not None]
    calibration = [parse_calibration_string(dataset[i]) for i in range(len(dataset)) if dataset[i] is not None]
    data_array = np.concatenate(data,axis=0)
    return statistics, data, data_array, calibration

def parse_calibration_string( dataset_item ):
    # why did I save the string as text? no clue. Easier to parse now than go back and also add some options data?
    cal = dataset_item[6]

    unimolecular_start = cal.find("Unimolecular: ")+14
    unimolecular_end   = cal.find(" Bimolecular: ")
    unimolecular = float(cal[unimolecular_start:unimolecular_end])

    bimolecular_start =  cal.find("Bimolecular: ")+13
    bimolecular_end   = cal.find(". Concentration: ")
    if bimolecular_end == -1:
        # old style without concentration
        bimolecular_end = cal.find(". \n Model:")
    bimolecular = float(cal[bimolecular_start:bimolecular_end])

    concentration_start = cal.find("Concentration: ")+15
    concentration_end   = cal.find("M\n Model:")
    if concentration_start == 14:
        concentration = 1.0
    else:
        concentration = float(cal[concentration_start:concentration_end])

    return unimolecular,bimolecular,concentration

def round_to_sig_digits( value, digits ):
    # -floor(log10(value)) is the position of the leftmost digit (single digit values are position 0, left is more negative, right is positive)
    # we then need to round starting at this position + (digits-1) to get the desired value

    return np.around( value, -int(np.floor(np.log10(value)))+(digits-1))

def compute_bimolecular_rate( temperature, length ):
    # Values for E_act and ln A from Morrison & Stols 1993
    # via equation 7:
    #
    # ln (k) = ln A - E_act / RT
    #
    # Values from Table IV for ln A and E_act:
    #
    # length 10: E_act: 9.94   ln A: 32.7
    # length 20: E_act: 16.4   ln A: 43.3

    R_gas_constant = 1.9858775e-3  #units of kcal/K/mol
    
    if length == 10:
        E_act = 9.94
        lnA = 32.7
    elif length == 20:
        E_act = 16.4
        lnA = 43.3

    if temperature < 100.0: #wtf, that's a low temperature for Kelvin.
        #which means, it's probably a mistake and they wanted Celsius.
        temperature = temperature + 273.15
        import warnings
        warnings.warn( "Temperature converted to {0} in Kelvin.".format( temperature))

    lnk = lnA - (E_act / (R_gas_constant * temperature))

    bi_full_rate = np.exp(lnk)

    # now round to 3 sig digits
    bi_rate = round_to_sig_digits( bi_full_rate, 3)
    
    return bi_rate
    

def compute_expected_time( temperature, length, concentration ):
    rate = compute_bimolecular_rate( temperature, length)

    stochastic_rate = rate * concentration
    
    return 1.0/stochastic_rate

def compute_k_bi_value( dataset, temperature, length ):
    stats, _, data, calibration = collate_result_list( dataset )

    data_expected_time = np.mean( data )
    means = [i[0] for i in stats]
    bimolecular = calibration[0][1]
    concentration = calibration[0][2]

    experimental_expected_time = compute_expected_time( temperature, length, concentration )
    k_bi_adjusted = bimolecular * (data_expected_time / experimental_expected_time)
    k_bi_values = [bimolecular * (i / experimental_expected_time) for i in means]
    k_bi_std = np.std(k_bi_values)

    return data_expected_time, experimental_expected_time, round_to_sig_digits( k_bi_adjusted, 3), k_bi_std, k_bi_std / np.sqrt( len( k_bi_values))



def final_values( params, length ):
    print("Bimolecular rate values:")
    lp = list_of_parameters()
    for i in params:
        param_values = lp[i]
        dangles = param_values[0]
        method = param_values[1]
        temperature = param_values[2]
        model = param_values[3]
        substrate = param_values[4]

        dataset = load_dataset( i, length,verbose=False )
        if len(dataset) > 0:
            
            data_et, exp_et, k_bi, k_bi_std, k_bi_err = compute_k_bi_value( dataset, temperature, length )

            uni, bi, conc = parse_calibration_string( [j for j in dataset if j is not None][0])
            print("Param set [{0}]:\n  Model: [{1}] Method: [{2}] Substrate: [{3}] Temperature: [{4}]\n  Simulation Parameters: k_uni: {5:.2e}  k_bi: {6:.2e}  Concentration: {7:.2e} (M)\n  Expected time (sim): {8:.3e}\n  Expected time (experiment): {9:.3e}\n  Calculated k_bi: {10:.2e} [std:{11:.3e}] [err:{12:.3e}]".format(i, model, method, substrate, temperature, uni, bi, conc, data_et, exp_et, k_bi, k_bi_std, k_bi_err))

def ratios( params, length ):
    lp = list_of_parameters()
    ratio_list = []
    for i in params:
        param_values = lp[i]
        dangles = param_values[0]
        method = param_values[1]
        temperature = param_values[2]
        model = param_values[3]
        substrate = param_values[4]

        dataset = load_dataset( i, length,verbose=False )
        if len(dataset) > 0:
            
            data_et, exp_et, _,_,_ = compute_k_bi_value( dataset, temperature, length )
            ratio_list.append( data_et / exp_et )

    return ratio_list



