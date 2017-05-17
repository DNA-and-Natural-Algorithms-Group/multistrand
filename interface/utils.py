import random
    
#FD, May 5th 2017
# Supplying rate options for Metropolis and Kawasaki methods,
# all using the dangles = some option. Also:  one general default,
# and one setting for Metropolis rates derived for DNA23. 

def JSDefault(options): 
    """ Default rates from Joseph Schaeffer's thesis  """
    
    options.unimolecular_scaling =1.50e+08;
    options.bimolecular_scaling = 1.38e+06;
    
    

def JSMetropolis25(options): 
    """ Default rates for Metropolis at 25 degree Celcius, from Joseph Schaeffer's thesis
    
    """
    
    options.unimolecular_scaling = 4.4e8;
    options.bimolecular_scaling = 1.26e6;
    

def JSKawasaki25(options): 
    """ Default rates for Kawasaki at 25 degree Celcius, from Joseph Schaeffer's thesis
    
    """
    
    options.unimolecular_scaling = 6.1e7;
    options.bimolecular_scaling = 1.29e6;
    
    

def JSKawasaki37(options):
    """ Default rates for Kawasaki at 37 degree Celcius, from Joseph Schaeffer's thesis
    """
    
    options.unimolecular_scaling = 1.5e8;
    options.bimolecular_scaling = 1.38e6;
    
     

def JSMetropolis37(options): 
    """ Default rates for Metropolis at 37 degree Celcius, from Joseph Schaeffer's thesis
    """
    
    options.unimolecular_scaling = 7.3e8;
    options.bimolecular_scaling = 1.40e6;
    
     

def DNA23Metropolis(options):
    """ A default rate for Metropolis at 25 degree Celcius, from the DNA23 conference
    """
        
    options.unimolecular_scaling = 5.0e6;
    options.bimolecular_scaling = 1.4e6;



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

