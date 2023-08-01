# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

import os, random
from functools import reduce

import numpy as np

from ._objects.strand import Strand
from nupack import mfe


GAS_CONSTANT = 0.0019872036  # kcal / K mol


def meltingTemperature(seq, concentration=1.0e-9):
    """
    Returns the melting temperature in Kelvin for a duplex of the given sequence.
    Sequences should be at least 8 nt long for the SantaLucia model to reasonably apply.
    For shorter sequences, see notes on "Melting Temperature (Tm) Calculation" by biophp.org
    Specifically, look at basicTm vs Base-stacking Tm. FD mar 2018
    """
    strand = Strand(sequence=seq)

    energy20 = (float(mfe([strand.sequence, strand.C.sequence ], material='dna', T=20.0)[0][1])
                + GAS_CONSTANT * (273.15 + 20) * np.log(55.5))
    energy30 = (float(mfe([strand.sequence, strand.C.sequence ], material='dna', T=30.0)[0][1])
                + GAS_CONSTANT * (273.15 + 30) * np.log(55.5))

    dS = (energy20 - energy30) / 10.0  # kcal/ K mol
    dH = energy30 + (273.15 + 30.0) * dS  # kcal/mol
    return  (dH / (dS + GAS_CONSTANT * np.log(concentration / 4.0)))


def concentration_string(concentration):
    """
    An easy print function to format concentration in M
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


def uniqueStateID(idsList, structsList):
    """
    Takes a list of ids and structures, and computes the pairtype for each of
    the complexes. Then returns a list of pairtypes that is alphabetically
    ordered.
    """
    pairTypes = []
    
    for ids, struct in zip(idsList, structsList):
        myPairType = pairType(ids, struct)
        pairTypes.append(myPairType)

    # now sort the list of lists by the first element of each list (which is 
    mySortedList = sorted(pairTypes, key=lambda x: x[0])
    
    # to make this hashable, we make it into a tuple.
    return tuple(mySortedList)


def generatePairing(dotParen, stack, offset, output) -> None:
    """
    Utility function.
    """
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
            raise Warning('generatePairing: There is an error in the dot paren structure.')
        index += 1


def pairType(ids, structs):
    """
    Given identifiers and dot-parens for a complex,
    pairType returns a unique identifier for that secondary structure.
    """
    idList = ids.split(',')
    if all(id_.count(":") == 0 for id_ in idList):
        pass
    elif all(id_.count(":") == 1 for id_ in idList):
        nList = []
        for id_ in idList:
            id_ = id_.split(":")[1]
            nList.append(id_)
        idList = nList
    else:
        raise ValueError("Unsupported strand labelling.")
    
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
    for index in range(len(idList)):
        generatePairing(dotParens[index], myStack, offsets[index], output)
    return (tuple(idString), tuple(output))


def generate_sequence(n, allowed_bases=['G', 'C', 'T', 'A'], base_probability=None):
    """ Generate a sequence of N base pairs.

    Bases are chosen from the allowed_bases [any sequence type], and
    according to the probability distribution base_probability - if
    none is specified, uses uniform distribution."""
    
    result = ""
    if base_probability is None:
        return result.join([random.choice(allowed_bases) for _ in range(n)])
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


def dGC_feature(o, i: int):
    states = o.full_trajectory[i]
    dGC = 0.0
    for state in states:
        dGC += (state[5] - (
            o._temperature_kelvin * 0.0019872036
            * np.log(1.0 / o.join_concentration) * state[4].count("+")))
        # print("count is  " +  str(state[4].count("+")))
        # print("join conc is " + str(o.join_concentration))
        # print("dG-Complex is " + "%.2f" % dGC + " kcal/mol  for " + str(state[3]))
    return f"dGC={dGC:> 6.2f} kcal/mol"


def printTrajectory(o, timescale=(1e3, "ms"), feature=None, show_uid: bool = False):
    seqstring = " "
    for i in range(len(o.full_trajectory)):
        time = timescale[0] * o.full_trajectory_times[i]
        states = o.full_trajectory[i]

        ids = []
        newseqs = []
        structs = []
        dG = 0.0

        pairTypes = []
        for state in states:
            ids.append(str(state[2]))
            # extract the strand sequences in each complex
            # (joined by "+" for multistranded complexes)
            newseqs.append(state[3])
            # similarly extract the secondary structures for each complex
            structs.append(state[4])
            dG += state[5]

            if show_uid:
                uniqueID = pairType(state[2], state[4])
                pairTypes.append(
                    ''.join(uniqueID[0]) + '_' + ','.join(map(str, uniqueID[1])))

        # make a space-separated string of complexes, to represent the whole
        # tube system sequence
        newseqstring = ' '.join(newseqs)
        # give the dot-paren secondary structure for the whole test tube
        tubestruct = ' '.join(structs)
        if show_uid:
            identities = '+'.join(pairTypes)

        if not newseqstring == seqstring:
            print(newseqstring)
            # because strand order can change upon association of dissociation,
            # print it when it changes
            seqstring = newseqstring

        print(f"{tubestruct}   t={time:.6f} {timescale[1]}, dG={dG:> 6.2f} kcal/mol"
              + (f", {feature(o, i)}" if feature is not None else "")
              + (f", uID='{identities}'" if show_uid else ""))
