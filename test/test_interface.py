###
### Note the changed stuff below to add the 'multistrand' package.
### If you run this file from a path other than the one it's located in, it may not
### be able to find the package correctly.
###

# fix relative paths so we can import multistrand package.
import os.path
import sys

if sys.path[0] == '':
  sys.path.append(os.path.realpath('../../'))
else:
  sys.path.append(os.path.realpath( os.path.join(sys.path[0], '../../')))

from multistrand.objects import Domain, Strand, Complex, StopCondition, RestingState
from multistrand.options import Options
from multistrand.system  import SimSystem

#from multiprocessing import Pool
#import cPickle
import random


def create_options():
  print("creating options...")
  d1 = Domain(name="d1", sequence="GTTGGTTTGTGTTTGGTGGG")
  s1 = Strand(name="s1", domains=[d1])
  c1 = Complex(name="c1", strands=[s1], structure=".")
  c2 = Complex(name="c2", strands=[s1.C], structure=".")
  c3 = Complex(name="c3", strands=[s1, s1.C], structure="(+)")
  
  sc_rev = StopCondition("REVERSE", [(c1, 2, 0), (c2, 2, 0)])
  sc_for = StopCondition("END", [(c3, 4, 6)])
  o = Options(simulation_mode = 'First Step', num_simulations = 100,
              simulation_time = 0.5, start_state = [RestingState("1", [c1]), RestingState("2", [c2])],
              stop_conditions = [sc_rev, sc_for])
  #o.initial_seed = random.SystemRandom().randrange(-2147483648, 2147483647)
  print("finished options. creating simsystem...")
  #print o.interface
  return o

def run_trajectory(o):
  s = SimSystem(o)
  print("finished simsystem. now starting...")
  s.start()
  print("after call to start.")
  
#  tag, time, rate = o.interface.current_tag, o.interface.current_time, o.interface.collision_rate
#  del s
#  return tag, time, rate


if __name__ == '__main__':
  num_sims = 1
  
  #pool = Pool()
  #results = [pool.apply_async(run_trajectory, ()) for i in range(num_sims)]
  #pool.close()
  #pool.join()
  #results = [r.get() for r in results]
  
  for i in range(num_sims):
    o = create_options()
    run_trajectory(o)
