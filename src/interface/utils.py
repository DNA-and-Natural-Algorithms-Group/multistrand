import random
import os
    



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



# ## Pairtype util
def pairType(ids, structs):
    """Given identifiers and dot-parens for a complex, 
        pairType returns a unique identifier for that secondary structure
    """
    
    # utility function
    def generatePairing(dotParen, size, stack, offset, output):
    
        index = 0
    
        for c in dotParen:
            
            if c == '(':
                # pushing the first end of the basepair
                stack.append(offset + index)
            elif c == ')':
                # popping the stack, setting two locations
                currIndex = offset + index
                otherIndex = stack.pop()
                                
                output[currIndex] = otherIndex
                output[otherIndex] = currIndex
                
            elif not c == '.':
                raise Warning('There is an error in the dot paren structure -- PairType function') 
            
            index += 1
    
    
    
    idList = ids.split(',')
    dotParens = structs.split('+') 
    
    outputSize = sum(len(x) for x in dotParens)
    

    offset = 0;    
    myStack = []
    output = [0, ] * outputSize

    myEnum = enumerate(idList)
    
    for index, val in myEnum:      
        
        size = len(dotParens[index])        
        generatePairing(dotParens[index], size, myStack, offset, output)
        
        offset += size

        
        
    
    indices = sorted(range(len(idList)), key=idList.__getitem__)
        
    # now weave the new pairTypes in the right order.
    
    actualOutput = []    
        
    idString = ''
    for index in indices:
        
        idString += idList[index]        
        size = len(dotParens[index]) 
        
        offset = sum([len(dotParens[idx]) for idx in range(0, index - 1) ])  
        
        
        actualOutput += output[offset : (offset + size)]
        
        
     
    
    return  (idString, tuple(output))
    



def generate_sequence( n, allowed_bases = ['G','C','T','A'], base_probability = None ):
    """ Generate a sequence of N base pairs.

    Bases are chosen from the allowed_bases [any sequence type], and
    according to the probability distribution base_probability - if
    none is specified, uses uniform distribution."""

    
    result = ""
    if base_probability is None:
        return result.join([random.choice( allowed_bases ) for i in range(n)])
    else:
        def uniform_seq( r ):
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
            return lambda x,y: (x[0]+y[0],y[1]) if r>=x[0] else (x[0],x[1])
        return result.join( [reduce( uniform_seq( random.random() ),
                              zip(base_probability,allowed_bases),
                              (0.0,'Invalid Probabilities'))[1] 
                              # note this subscript [1] pulls out the item
                              # selected by the reduce since the result was a tuple.
                      for i in range(n)]
                     )

