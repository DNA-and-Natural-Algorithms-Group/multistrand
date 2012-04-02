import sys,os,os.path
import cPickle
import subprocess # for Popen to call some nupack functions

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

except ImportError:
    print("Could not import Multistrand.")
    raise

kBoltzmann = .00198717 #units of kcal/(mol*K)
RT = kBoltzmann * (273.15+37.0)

def setup_options( seq, concentration ):
    """
    setup_options( seq )

    creates an Options object using the sequence passed as a single
    domain with initially unpaired structure. 
    """

    d = Domain(name="initial",sequence=seq,length=len(seq))
    s = Strand(domains=[d])
    c = Complex(strands=[s], structure=".")
    c2 = Complex(strands=[s.C], structure=".")
    

    o = Options(simulation_mode="Normal",parameter_type="Nupack",substrate_type="DNA",
                num_sims = 100, sim_time=10.0, start_state=[c,c2])

    o.temperature = 310.15
    o.join_concentration = concentration
    return o


def nupack_pfunc( sequence, temperature ):
    #time.sleep(.02)
    p = subprocess.Popen(["pfunc", "-material", "dna","-T","{0}".format(temperature)], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    input_str = "{seq}\n".format( seq=sequence )
    
    stdout_chars = p.communicate(input_str)[0]
    flag = False
    energy = None
    pfunc = None
    # value we want is on the line after "% Energy"
    result = stdout_chars.split("\n")
    for l in result:
        if flag and energy is None:
            energy = float( l )
        elif flag:
            pfunc = float( l )
            flag = False
        if l.startswith("% Free energy"):
            flag = True
    return (energy,pfunc)

def run_distribution( seq ):
    o = setup_options(seq, 2*math.exp(-5.06 / RT))
    s = SimSystem(o)
    s.start()
    eq_dict = {}
    
    for end_state in o.interface.end_states:
        for cmplx in end_state:
            if cmplx[4] in eq_dict:
                count = eq_dict[cmplx[4]][1]
                eq_dict[cmplx[4]][1] = count + 1
            else:
                eq_dict[cmplx[4]] = [cmplx[5], 1]
    return eq_dict

        
# NOTE: should also look at what % of total PF was covered by the
# states seen.

# some input lines:

# import os
# os.environ['MULTISTRANDHOME'] = "/Users/zifnab/Projects/MultistrandMainClone"

# eqd = run_distribution("GGGGGGAAAACCCCCC")
# eqd
# eqd.items()
# import math
# pf = nupack_pfunc("GGGGGGAAAACCCCCC",37.0)
# equiv = [(i[0],math.exp(-i[1][0]/RT)/pf, float(i[1][1])/1000.0) for i in eqd.items()]
# equiv
# eq_testing_1 = ("1000 sims, .01s",eqd,equiv)

# eqd2 = run_distribution("GGGGGGAAAACCCCCC")
# equiv2 = [(i[0],math.exp(-i[1][0]/RT)/pf, float(i[1][1])/1000.0) for i in eqd2.items()]
# equiv2
# total = 0.0
# for a,b,c in equiv2:
#     total += abs(b-c)
# total
# total = 0.0
# for a,b,c in equiv:
#     total += abs(b-c)
# total
# eq_testing_2 = ("1000 sims, .1s",eqd2,equiv2)
# eq_testing_1 = ("1000 sims, .01s",eqd,equiv)

# eqd3 = run_distribution("GGGGGGAAAACCCCCC")
# eqd3
# equiv3 = [(i[0],math.exp(-i[1][0]/RT)/pf, float(i[1][1])/10000.0) for i in eqd3.items()]
# equiv3
# total = 0.0
# for a,b,c in equiv3:
#     total += abs(b-c)
# total
# eq_testing_3 = ("10000 sims, .01s",eqd3,equiv3,total)
# f = open('three_datasets.dat','wb')
# import cPickle
# cPickle.dump((eq_testing_1, eq_testing_2, eq_testing_3),f,-1)
# f.close()

# eqd4 = run_distribution("GGGGGGAAAACCCCCC")
# eqd4
# equiv4 = [(i[0],math.exp(-i[1][0]/RT)/pf, float(i[1][1])/100000.0) for i in eqd4.items()]
# equiv4
# total = 0.0
# for a,b,c in equiv4:
#     total += abs(b-c)
# total
# probtot = 0.0
# for a,b,c in equiv4:
#     probtot += b
# probtot
# probtot2 = 0.0
# for a,b,c in equiv2:
#     probtot2 += b
# probtot
# probtot2
# probtot4 = probtot
# probtot = 0.0
# for a,b,c in equiv:
#     probtot += b
# probtot3 = 0.0
# for a,b,c in equiv3:
#     probtot3 += b
# probtot3
# probtot2
# probtot
# probtot4
# eq_testing_4 = ("100000 sims, .01s",eqd4,equiv4,total)
# f = open('four_datasets.dat',"wb")
# cPickle.dump((eq_testing_1, eq_testing_2, eq_testing_3,eq_testing_4),f,-1)
# f.close()
# f = open('bad_equilibrium.dat','rb')
# dd = cPickle.load( f )
# dd
# equivbad = [(i[0],math.exp(-i[1][0]/RT)/pf, float(i[1][1])/10000.0) for i in dd.items()]
# equivbad
# eq_testing_1
# eq_testing_2
# eq_testing_3
# eq_testing_4
