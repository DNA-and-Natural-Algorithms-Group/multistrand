
# Code contributors:  Erik Winfree, Chris Thachuk, Justin Bois.

import math
import subprocess as sub

# From Jason at Stackoverflow:
class Memoize:
    def __init__(self, f):
        self.f = f
        self.memo = {}
    def __call__(self, *args):
        if not args in self.memo:
            self.memo[args] = self.f(*args)
        return self.memo[args]
# Note that the memoization here is useful when algorithms make repeated calls to evaluate the same sequences.
# It comes at the expense of using an unbounded amount of memory.
# A more sophisticated memoization would limit the number of memoized items, and forget the older ones.

# from flebool on Stackoverflow
def index_max(values):
    """Gives the index of the maximum element."""
    return max(xrange(len(values)),key=values.__getitem__)

def index_min(values):
    """Gives the index of the minimum element."""
    return min(xrange(len(values)),key=values.__getitem__)

def cyclic_sort(mylist):
    """Rotates the list to put the "least" element first."""
    i=index_min(mylist)
    return mylist[i:]+mylist[:i]

def rev(dom) : 
    """Reverses a string."""
    return dom[::-1]

def WC(seq):
    """Computes the Watson-Crick complement for a DNA sequence."""
    return seq.translate(string.maketrans('ACTG','TGAC'))[::-1]

def dGadjust(T,N):
    """Adjust NUPACK's native free energy (with reference to mole fraction units) to be appropriate for molar units, assuming N strands in the complex."""
    R=0.0019872041 # Boltzmann's constant in kcal/mol/K 
    water=55.14    # molar concentration of water at 37 C, ignore temperature dependence, which is about 5%
    K=T+273.15     # Kelvin
    adjust = R*K*math.log(water) # converts from NUPACK mole fraction units to molar units, per association
    return adjust*(N-1)



def energy(seq,structure,material,T):
    """Calls NUPACK's 'energy' for a given sequence and secondary structure pair.  Returns the microstate dG.  T is in Celsius."""
    nupack_input = str(1) + '\n' + seq + '\n' + str(1) + '\n' + structure + '\n'
    p = sub.Popen(['energy', '-T', str(T), '-multi', '-material', material], stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.PIPE)
    output = p.communicate(nupack_input)[0]
    lines = output.split('\n')

    if lines[-3] != "% Energy (kcal/mol):" :
       raise ValueError('NUPACK output parsing problem')

    return float(lines[-2])


@Memoize
def pfunc_tuple(seqtuple,T):
    """Calls NUPACK's pfunc on a complex consisting of the unique strands in seqtuple, returns dG.  T is in Celsius."""
    user_input = str(len(seqtuple)) + '\n' + '\n'.join(seqtuple) + '\n' + str(range(1,len(seqtuple)+1))[1:-1].replace(',','')
    p=sub.Popen(['pfunc','-T',str(T),'-multi','-material','dna'],stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
    output,error = p.communicate(user_input)
    lines = output.split('\n')

    while len(lines) < 4 : # can't figure out why, but occasionally NUPACK returns empty-handed.  Subsequent tries seem to work...
        print 'Retrying in pfunc_tuple: NUPACK failed with output ' + `lines` + ' and error ' + `error` +" ."
        p=sub.Popen(['pfunc','-T',str(T),'-multi','-material','dna'],stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
        output,error = p.communicate(user_input)
        lines = output.split('\n')

    if lines[-4] != "% Free energy (kcal/mol) and partition function:" :
        raise NameError('NUPACK output parsing problem')

    if float(lines[-3])==float('inf') : return 0               # if these strands can't base-pair
    else: return float(lines[-3]) + dGadjust(T,len(seqtuple))

# We have pfunc_tuple because tuples but not lists can be memoized.  But lists are clearer notationally.
# An alternative syntax to the "@Memoize" line would be, here: pfunc_tuple=Memoize(pfunc_tuple)
# Note that since pfunc give identical results for any rotation of strand order, we use a canonical order that aids memoization.

def pfunc(seqlist,T) : return pfunc_tuple(tuple(cyclic_sort(seqlist)),T)

@Memoize
def duplex_energy(seq,T):
    """Calls NUPACK's 'energy' on a complex consisting of two complementary strands, returns the microstate dG.  T is in Celsius."""
    seqtuple=(seq,WC(seq))
    user_input = str(len(seqtuple)) + '\n' + '\n'.join(seqtuple) + '\n' + str(range(1,len(seqtuple)+1))[1:-1].replace(',','') + '\n' + '('*len(seq)+'+'+')'*len(seq)
    p=sub.Popen(['energy','-T',str(T),'-multi','-material','dna'],stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
    output = p.communicate(user_input)[0]
    lines = output.split('\n')

    while len(lines) < 3 : # can't figure out why, but occasionally NUPACK returns empty-handed.  Subsequent tries seem to work...
        print 'Retrying in duplex_energy: NUPACK failed with output ' + `lines` + ' and error ' + `error` +" ."
        p=sub.Popen(['energy','-T',str(T),'-multi','-material','dna'],stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
        output,error = p.communicate(user_input)
        lines = output.split('\n')

    if lines[-3] != "% Energy (kcal/mol):" :
       raise ValueError('NUPACK output parsing problem')

    return float(lines[-2]) + dGadjust(T,len(seqtuple))

# PyUtils has:
# getEnergy    energy
# getProb      prob
# getMFEStruct mfe
# getPairProbs pairs

# we want:
# sample
# defect
# energy --- need to allow structures other than duplex // multistranded things.  add salt.
# pfunc  --- need to make sure it works for any number of complexes.  add material and salt.
# prob
# mfe
# pairs

# other nupack functions that aren't as high priority right now, are:
# complexes
# count
# concentrations
# design
# distributions
# subopt
