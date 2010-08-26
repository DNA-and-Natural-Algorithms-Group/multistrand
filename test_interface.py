from python_objects import Domain, Strand, Complex, StopCondition, RestingState
import python_options
import multistrand


print "creating options..."
d1  = Domain("d1", "d", 5, False)
d1p = Domain("d1p", "d", 5, True)
s1  = Strand("s1", "s1",  "ACTTG", [d1])
s2  = Strand("s2", "s2",  "CAAGT", [d1p])
c1 = Complex("c1", "c1", [s1], ".....")
c2 = Complex("c2", "c2", [s2], ".....")
c3 = Complex("c3", "c3", [s1, s2], "(((((+)))))")
rev_tag = "REVERSE"
sc_rev = StopCondition(rev_tag, [(c1, 2, 0), (c2, 2, 0)])
for_tag = "END"
sc_for = StopCondition(for_tag, [(c3, 4, 1)])
o = python_options.MultistrandOptions()
o.simulation_mode = 3
o.use_stop_states = True
o.parameter_type  = 1
o.substrate_type  = 2
o.num_simulations = 3
o.simulation_time = 0.5
#o.initial_seed = 2
o.start_state = [c1, c2]
#o.start_state = [RestingState("1", [c1]), RestingState("2", [c2])]
o.stop_conditions = [sc_rev, sc_for]
print "finished options. creating simsystem..."
s = multistrand.SimSystem(o)
print "finished simsystem. now starting..."
s.start()
print "after call to start."
print o.interface


