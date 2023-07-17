
from multistrand.objects import Complex, Domain, Strand, StopCondition
from multistrand.options import Options
from multistrand.system import SimSystem

o1 = Options(temperature=25, dangles="Some")
    
o1.num_simulations = 1
o1.output_interval = 1 
o1.simulation_time = 10.0
o1.verbosity = 3

mySeq = "GCGTTTCAC"

# Using domain representation makes it easier to write secondary structures.
onedomain = Domain(name="itall", sequence=mySeq)
top = Strand(name="top", domains=[onedomain])
bot = top.C

# Note that the structure is specified to be single stranded, but this will be over-ridden when Boltzmann sampling is turned on.
startTop = Complex(strands=[top], structure=".")
startBot = Complex(strands=[bot], structure=".")

o1.start_state = [startTop, startBot]

# Stop when the exact full duplex is achieved.
success_complex = Complex(strands=[top, bot], structure="(+)")
stopSuccess = StopCondition("success", [(success_complex, 0, 0)])

# Declare the simulation unproductive if the strands become single-stranded again.
failed_complex = Complex(strands=[top], structure=".")
stopFailed = StopCondition("failure", [(failed_complex, 2, 0)])

o1.stop_conditions = [stopSuccess, stopFailed]
o1.initial_seed = 1783

s = SimSystem(o1)
s.start()
