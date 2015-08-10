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
    rootpath = os.path.abspath(os.getcwd()[:idx])
    #abspath cleans up the tail of the pathname as needed.
    if rootpath not in sys.path:
        sys.path.append( rootpath )
elif multihome != None:
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
    #note: could read off the options object here as well
    f.close()

    return res

def load_dataset( num, length ):
    dirname = filename_from_param_number( num )
    if not os.path.isdir( dirname ):
        print("Error: directory {0} does not exist.\n".format( dirname ))
        return None

    files = os.listdir( dirname )
    sample_files = [i for i in files if i.endswith('.result.dat') and i.startswith('len_{0}'.format(length))]
    sample_files.sort()
    
    def number_from_fname( filename ):
        return filename.rstrip('.result.dat').split('_')[-1]

    result_list = [None] * 1000 #len(sample_files)
    for f in sample_files:
        print("Processing file: {0}\n".format( f ))
        result_list[ int( number_from_fname(f)) ] = process_dataset(load_file( dirname, f ))

    return result_list

def process_dataset( dataset ):
    forward = [i for i in dataset if i.tag == "forward_bm"]
    forward_array = np.zeros( len(forward))
    forward_array[:] = [i.time for i in forward]

    std = np.std(forward_array,ddof = 1)
    err = std / np.sqrt( len( forward) )
    return (np.mean(forward_array), std, err, len(forward), len(dataset), forward_array)


def process_sample_set( samples ):
    return [process_dataset(i) for i in samples if i]

def aggregate_data( processed_samples ):
    no_nones = [i for i in processed_samples if i is not None]
    means = np.zeros( len( no_nones ))
    stds = np.zeros( len(no_nones ))
    means[:] = [i[0] for i in no_nones]
    stds[:] = [i[1] for i in no_nones]

    return (np.mean(means), np.std(means, ddof=1), np.mean( stds ), np.std( stds, ddof=1))


def compute_unimolecular( raw, length ):
    target = 11.6e-6 #seconds, e.g. 11.6 microseconds
    steps = length * length

    step_time = raw / float(steps)
    ratio = step_time / target
    return ratio
    

def compile_result_data( params, length ):
    results = []
    datasets = []
    for p_number in params:
        dataset = load_dataset( p_number, length)
        results.append( (p_number, aggregate_data(dataset) ))
        datasets.append( (p_number, dataset ))

    return results, datasets
            


# calibration on means


###Sample usage:
# __IP.Completer.omit__names = 1
# _ip.magic("run analysis_unimolecular.py")
# d0 = load_dataset( 0, 10)
# len(d0)
# d0_results = aggregate_data(d0)
# d0_results
# filename_from_param_number(0)
# d1 = load_dataset( 1, 10)
# d1_results = aggregate_data(d1)
# d1_results
# filename_from_param_number(1)
# d2 = load_dataset( 2, 10)
# d2_results = aggregate_data(d2)
# d2_results
# _ip.magic("run analysis_unimolecular.py")
# res, ds = compile_result_data( [0,1,2,3,4,5], 10 )

# import matplotlib
# matplotlib.use('macosx')
# import matplotlib.pylab as pyp

# ds0_means = [i[0] for i in ds[0][1] ]
# ds0_means
# pyp.hist(ds0_means)
# help(pyp.hist)
# pyp.hist(ds0_means, bins=50 )
# ds1_means = [i[0] for i in ds[1][1] ]
# pyp.hist(ds1_means, bins=50 )
# ds2_means = [i[0] for i in ds[2][1] ]
# pyp.hist(ds1_means, bins=50 )
# pyp.hist(ds2_means, bins=50 )
# ds3_means = [i[0] for i in ds[3][1] ]
# pyp.hist(ds3_means, bins=50 )
# ds4_means = [i[0] for i in ds[4][1] ]
# pyp.hist(ds4_means, bins=50 )
# ds5_means = [i[0] for i in ds[5][1] ]
# ds5_means = [i[0] for i in ds[5][1] if i is not None]
# pyp.hist(ds5_means, bins=50 )
# res
# compute_unimolecular( res[0][1][0], 10 )
# compute_unimolecular( res[1][1][0], 10 )
# compute_unimolecular( res[2][1][0], 10 )
# compute_unimolecular( res[3][1][0], 10 )
# compute_unimolecular( res[4][1][0], 10 )
# compute_unimolecular( res[5][1][0], 10 )
# _ip.system("ls -F ")
# res_20, ds_20 = compile_result_data( [0,1,2,3,4,5], 20 )
# res_20
# for i in res:
#     print "{0}:  {1} (time) {2} (uni_scale)".format( filename_from_param_number( i[0]), i[1][0], compute_unimolecular( i[1][0], 10 ) )
# print In[:]
