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

from multistrand.objects import *
from multistrand.options import Options
from multistrand.system import run_system

o = Options()

o.simulation_mode = 128
o.num_simulations = 1
o.simulation_time = .001
o.output_interval = 1

d = Domain(name="d1",sequence="AGCACATTAG",length=10)
d2 = Domain(name="hp",sequence="AAAAA", length=5)
s1 = d + d2 + d.C
s1.name = "Hairpin"

start_c = Complex(strands=[s1],structure="...")

o.start_state = [start_c]


SA_t = Domain(name="SA_toe",sequence="ACCCAC",length=6)
SB_t = Domain(name="SB_toe",sequence="GCTAAC",length=6)
HA_t = Domain(name="HA_toe",sequence="GTGGGT",length=6)
HB_t = Domain(name="HB_toe",sequence="GTTAGC",length=6)
SA_t2 = Domain(name="SA_toe2",sequence = "ACC",length=3)
SA_d = Domain(name="SA",sequence="GCACGTCCACGGTGTCGC",length=18)
SB_d = Domain(name="SB",sequence="GCGACACCGTGGACGTGC",length=21)
HA_d = Domain(name="HA",sequence="GCGACACCGTGGACGTGC",length=18)
HB_d = Domain(name="HB",sequence="GCACGTCCACGGTGTCGC",length=18)
SA = SA_t2 + SA_d + SA_t
SB = SB_t + SA_d.C + SA_t2.C
HA = SA_t.C + HA_d
HB = HA_d.C + SB_t.C

Start_S = Complex(strands=[SB,SA],structure=".((+)).")
Start_H = Complex(strands=[HA,HB],structure=".(+).")
Final = Complex(strands=[SB,SA,HA,HB],structure="(((+)((+))+))")
stopc = StopCondition("Final",[ (Final,0,0)])
o2 = Options()
o2.simulation_mode = 128 #16
o2.num_simulations = 1
o2.simulation_time = .1
o2.initial_seed = 0xf67db97dedf3070
o2.output_interval = 100
#o2.output_interval = 1
o2.start_state = [Start_H, Start_S]
o2.stop_conditions = [stopc]
