
import os, random
from functools import reduce

import numpy as np

from ._objects.strand import Strand
from nupack import mfe

""" Returns the melting temperature in Kelvin for a duplex of the given sequence.
    Sequences should be at least 8 nt long for the SantaLucia model to reasonably apply.
    For shorter sequences, see notes on "Melting Temperature (Tm) Calculation" by biophp.org
    Specifically, look at basicTm vs Base-stacking Tm. FD mar 2018
"""


def meltingTemperature(seq, concentration=1.0e-9):
        GAS_CONSTANT = 0.001987
        
        strand = Strand(sequence=seq)

        energy20 = float(mfe([strand.sequence, strand.C.sequence ], material='dna', T=20.0)[0][1]) + GAS_CONSTANT * (273.15 + 20) * np.log(55.5)
        energy30 = float(mfe([strand.sequence, strand.C.sequence ], material='dna', T=30.0)[0][1]) + GAS_CONSTANT * (273.15 + 30) * np.log(55.5)
        
        dS = (energy20 - energy30) / 10.0  # kcal/ K mol
        dH = energy30 + (273.15 + 30.0) * dS  # kcal/mol
        return  (dH / (dS + GAS_CONSTANT * np.log(concentration / 4.0)))


def concentration_string(concentration):
    """ An easy print function to format concentration in M
    """
    if concentration < 1e-12: 
        return "{} fM".format(1e15 * concentration)
    if concentration < 1e-9: 
        return "{} pM".format(1e12 * concentration)
    if concentration < 1e-6: 
        return "{} nM".format(1e9 * concentration)
    if concentration < 1e-3: 
        return "{} uM".format(1e6 * concentration)
    if concentration < 1: 
        return "{} mM".format(1e3 * concentration)
    return "{} M".format(concentration)


def seqComplement(sequence):
    complement = {'G':'C',
                  'C':'G',
                  'A':'T',
                  'T':'A'}
    return "".join([complement[i] for i in reversed(sequence)])


def standardFileName(SCRIPT_DIR, mySeq=None, extraTitle=None, runs=None):
    fileName = str(SCRIPT_DIR)
    
    if not mySeq == None:
        fileName += str("/" + mySeq + '/' + mySeq)
    else:
        fileName += "/"
        
    if not runs == None:
        fileName += "-" + str(runs) 
    if not extraTitle == None:
        fileName += "-" + extraTitle
    if not os.path.exists(os.path.dirname(fileName)):
        try:
            os.makedirs(os.path.dirname(fileName))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    return fileName


# Takes a list of ids and structures, and computes the pairtype for each of the complexes. 
# Then returns a list of pairtypes that is alphabetically ordered.
def uniqueStateID(idsList, structsList):
    pairTypes = []
    
    for ids, struct in zip(idsList, structsList):
        myPairType = pairType(ids, struct)
        pairTypes.append(myPairType)

    # now sort the list of lists by the first element of each list (which is 
    mySortedList = sorted(pairTypes, key=lambda x: x[0])
    
    # to make this hashable, we make it into a tuple.
    return tuple(mySortedList)


# ## Pairtype util
def pairType(ids, structs):
    """Given identifiers and dot-parens for a complex, 
        pairType returns a unique identifier for that secondary structure
    """

    # utility function
    def generatePairing(dotParen, stack, offset, output):
        index = 1
        for c in dotParen:
            if c == '(':
                # pushing the first end of the basepair
                stack.append(offset + index)
            elif c == ')':
                # popping the stack, setting two locations
                currIndex = offset + index
                otherIndex = stack.pop()
                                
                output[currIndex - 1] = otherIndex
                output[otherIndex - 1] = currIndex
            elif not c == '.':
                raise Warning('There is an error in the dot paren structure -- PairType function') 
            index += 1
    
    idList = ids.split(',')
    nList = []
    for id in idList:
        id = id.split(":")[1]
        nList.append(id)
    idList = nList
    
    dotParens = structs.split('+') 
    N = len(dotParens)
    
    # the new ordering, for example: 3 0 1 2, so that idList[3] < idList[0] < idList[1] < idList[2]
    ordering = sorted(range(len(idList)), key=idList.__getitem__)
    
    idString = ''        
    
    newLengths = [len(dotParens[ordering[i]]) for i in range(N)]   
    newOffsets = [ sum(newLengths[0:i]) for i in range(N)]  # the offsets under the new ordering
    offsets = [0, ] * N  # the new offsets under the old ordering
    for i in range(N):
        newPosition = ordering[i]
        offsets[newPosition] = newOffsets[i]
        idString += idList[newPosition]

    myStack = []
    output = [0, ] * sum([len(dp) for dp in dotParens])
    myEnum = enumerate(idList)
    for index, val in myEnum:
        generatePairing(dotParens[index], myStack, offsets[index], output)
#     str_output = ''.join([str(x) for x in output])
#     return (idString, str_output)
    
    return (tuple(idString), tuple(output))


def generate_sequence(n, allowed_bases=['G', 'C', 'T', 'A'], base_probability=None):
    """ Generate a sequence of N base pairs.

    Bases are chosen from the allowed_bases [any sequence type], and
    according to the probability distribution base_probability - if
    none is specified, uses uniform distribution."""
    
    result = ""
    if base_probability is None:
        return result.join([random.choice(allowed_bases) for i in range(n)])
    else:
        def uniform_seq(r):
            """ This function returns a lambda to be used on a sequence of tuples of
            (p,item) e.g. via reduce(uniform_seq(.75), seq, (0.0,'none')).

            In this situation, the lambda computes (p',item'), where
            p' is the sum of consecutive probabilities in the sequence
            <= r, and item' is the corresponding item.

            It achieves this by updating the result with the new sum
            of probabilities and last item checked until the the
            summed probability has exceeded the target.

            Note: r should be in [0.0,1.0) as produced by random(),
            any input r >= 1.0 (assuming the sequence given sums to
            1.0) will give the final element as the result which is technically incorrect.
            """
            return lambda x, y: (x[0] + y[0], y[1]) if r >= x[0] else (x[0], x[1])

        return result.join([
                reduce(uniform_seq(random.random()),
                       zip(base_probability, allowed_bases),
                        # note this subscript [1] pulls out the item selected by
                        # the reduce since the result was a tuple.
                       (0.0, 'Invalid Probabilities'))[1]
                for _ in range(n)])
