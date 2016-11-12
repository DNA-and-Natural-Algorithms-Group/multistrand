import sys,os,os.path
import cPickle

multihome = None
if 'MULTISTRANDHOME' in os.environ:
    if not os.path.isfile( os.path.join( os.environ['MULTISTRANDHOME'], 'setup.py') ):
        warnings.warn( ImportWarning("Could not find the file 'setup.py' in your MULTISTRANDHOME [{0}]; this environment variable is possibly out of date or not referring to the new Mulistrand distribution."))
        multihome=None
    else:
        if os.environ['MULTISTRANDHOME'] not in sys.path:
            multihome= os.environ['MULTISTRANDHOME']     

if multihome != None:
    sys.path.append(multihome)

try:

    from multistrand.objects import *
    from multistrand.options import Options
    from multistrand.system import SimSystem

    # import multiprocessing
    # from multiprocessing import Pool
except ImportError:
    print("Could not import Multistrand.")
    raise




def create_setup():

    toehold_seq = "GTGGGT"

    bm_design_A = "ACCGCACGTCCACGGTGTCG"

    bm_design_B = "ACCGCACGTCACTCACCTCG"


    toehold = Domain(name="toehold",sequence=toehold_seq,length=6)
    branch_migration_A = Domain(name="bm_A", sequence=bm_design_A,
                                seq_length=20)

    branch_migration_B = Domain(name="bm_B", sequence=bm_design_B,
                                seq_length=20)
    
    substrate_A = toehold + branch_migration_A
    substrate_B = toehold + branch_migration_B
    incumbent_A = Strand(name="incumbent",domains=[branch_migration_A.C])
    incumbent_B = Strand(name="incumbent",domains=[branch_migration_B.C])

    incoming_A = substrate_A.C
    incoming_B = substrate_B.C

    start_complex_A = Complex(strands=[incoming_A, substrate_A, incumbent_A],
                              structure=".(+)(+)")
    start_complex_B = Complex(strands=[incoming_B, substrate_B, incumbent_B],
                              structure=".(+)(+)")
    

    o = Options()

    o.simulation_mode = 0x0080 # trajectory mode
    o.num_simulations = 1
    o.simulation_time = .05
    o.temperature = 37.0
    o.dangles = 1
    o.output_interval = 100
    o.start_state = [start_complex_A]


    o2 = Options()

    o2.simulation_mode = 0x0080 # trajectory mode
    o2.num_simulations = 1
    o2.simulation_time = .05
    o2.temperature = 37.0
    o2.dangles = 1
    o2.start_state = [start_complex_B]
    o2.output_interval = 100
    
    return o,o2



# Resulting trajectory pieces used in MS thesis.
# Note that in the thesis, we displayed them in 5'->3'
# counterclockwise. Oops. 

# In [156]: o.full_trajectory_times[-1]
# Out[156]: 0.04999384007252141

# In [157]: o.full_trajectory[-1]
# Out[157]: 
# [(-1641083587287114401,
#   10,
#   '14:Automatic_28*,28:Automatic_28,30:incumbent',
#   'CGACACCGTGGACGTGCGGTACCCAC+GTGGGTACCGCACGTCCACGGTGTCG+CGACACCGTGGACGTGCGGT',
#   '....((.(((....))).)).(((((+))))).((((((((((((((((((((+))))))))))))))))))))',
#   -34.180000000000007)]

# In [158]: o.full_trajectory[-4994]
# Out[158]: 
# [(-1641083587287114401,
#   0,
#   '14:Automatic_28*,28:Automatic_28,30:incumbent',
#   'CGACACCGTGGACGTGCGGTACCCAC+GTGGGTACCGCACGTCCACGGTGTCG+CGACACCGTGGACGTGCGGT',
#   '....(((((......)))))((((..+..)))).(((((((((((((((((((+))))))))))))))))))).',
#   -37.390000000000008)]

# In [159]: o.full_trajectory_times[-4994]
# Out[159]: 0.0099962550721091736

# In [215]: o2.full_trajectory_times[-4341]
# Out[215]: 0.0099989013661066287

# In [216]: o2.full_trajectory[-4341]
# Out[216]: 
# [(-3633056388795064589,
#   40,
#   '37:Automatic_37,19:Automatic_37*',
#   'GTGGGTACCGCACGTCACTCACCTCG+CGAGGTGAGTGACGTGCGGTACCCAC',
#   '((((((((((((((((((((((((((+))))))))))))))))))))))))))',
#   -37.070000000000007),
#  (-3633056388795064589,
#   39,
#   '39:incumbent',
#   'CGAGGTGAGTGACGTGCGGT',
#   '((..........))......',
#   0.77999999999999992)]

# In [219]: o2.full_trajectory_times[-2]
# Out[219]: 0.049990959322094936

# In [220]: o2.full_trajectory[-2]
# Out[220]: 
# [(-3633056388795064589,
#   434,
#   '37:Automatic_37,19:Automatic_37*',
#   'GTGGGTACCGCACGTCACTCACCTCG+CGAGGTGAGTGACGTGCGGTACCCAC',
#   '((((((((((((((((((((((((((+))))))))))))))))))))))))))',
#   -37.070000000000007),
#  (-3633056388795064589,
#   433,
#   '39:incumbent',
#   'CGAGGTGAGTGACGTGCGGT',
#   '....................',
#   0.0)]
