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

from multistrand.system import *
from multistrand.options import *
from multistrand.objects import *

o = Options(num_simulations=1, simulation_time=10.0, start_state=[Complex(sequence='ACGATGTGGAGGCCAAGACATCGATAGCCGAAGCCAAACTCGGCACAGTC', structure='.'*50)], bimolecular_scaling=1.0, unimolecular_scaling=1.0)

def run(o):
	s = SimSystem(o)
	s.start()

