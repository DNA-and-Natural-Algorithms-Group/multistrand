import random


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
