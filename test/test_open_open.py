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
from multistrand.system import *
dd = Domain(name="short",sequence="TGCCGGTAT")
st = Strand(name="base", domains=[dd])
c = Complex( strands=[st,st.C], structure="....(....+....)....")
c2 = Complex( strands=[st], structure=".")
stt = StopCondition( "break", [(c2,2,0)])
o = Options( dangles = "All")
#o.start_state = [c]
#o.stop_conditions = [stt]

o.simulation_time = .01
#run_system(o)

polyA = Domain(name="polyA",sequence="AAAA")
singleG = Domain(name="G", sequence="G")

st1 = polyA + singleG + polyA
st2 = polyA + singleG.C + polyA + singleG.C + polyA
st3 = polyA + singleG + polyA

c_start = Complex( strands = [st1,st2,st3], structure=".(.+.).(.+.).")
c_end1 = Complex(strands=[st1], structure = "...")
c_end2 = Complex(strands=[st3], structure = "...")
sc_end1 = StopCondition("break1", [(c_end1,2,0)])
sc_end2 = StopCondition("break2", [(c_end2,2,0)])

o.start_state = [c_start]
o.stop_conditions = [sc_end1, sc_end2]


