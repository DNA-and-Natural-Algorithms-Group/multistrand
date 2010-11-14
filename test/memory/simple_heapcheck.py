from multistrand.system import *
from multistrand.options import *
from multistrand.objects import *

o = Options(num_simulations=1, simulation_time=10.0, start_state=[Complex(sequence='ACGATGTGGAGGCCAAGACATCGATAGCCGAAGCCAAACTCGGCACAGTC', structure='.'*50)], bimolecular_scaling=1.0, unimolecular_scaling=1.0)

def run(o):
	s = SimSystem(o)
	s.start()

