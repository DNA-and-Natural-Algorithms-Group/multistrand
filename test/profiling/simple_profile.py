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

from multistrand.options import Options, Constants
from multistrand.system import SimSystem
from multistrand.objects import Strand, Complex

def helper_create_Multistrand_options( sequence, time, count):
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

def helper_RunSingle_Multistrand_I( sequence, time, count=1 ):
    input_o = helper_create_Multistrand_options( sequence, time, count)

    s = SimSystem( input_o )
    s.start()
    output_multistrand = input_o.interface.results
    del s

    return output_multistrand

