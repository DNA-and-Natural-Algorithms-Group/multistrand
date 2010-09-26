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

class energy( object ):
    """ Gets the energy of a particular sequence and structure, or a list of them.

    If passed two non keyword args, they'll be considered the sequence and structure to compute.
    Keywords: 'sequence','structure': for sequence and structure of a single state.
              'sequence_list','structure_list': for sequence and structures of multiple states
              'model':  'Vienna', 'Nupack' or 'Multistrand'
                         (Note that Multistrand indicates the usage of one of the other
                         models, depending on what you initialized multistrand.system with.
                         The energy calculation is done by Multistrand natively.)
    """
    def __init__(self, sequence=False, structure=False,
                 sequence_list=False, structure_list=False, model=None ):
        if sequence and structure:
            # single energy, no mfe
            sequence_list=   [sequence]
            structure_list=  [structure]
        elif sequence:
            # single energy, mfe
            sequence_list = [sequence]
            structure_list = False

        if sequence_list and structure_list:
            # multiple energy, no mfe
            g = self.energy_generator( sequence_list, structure_list, model )
        elif sequence_list:
            g = self.mfe_generator( sequence_list, None, model )
            # multiple energy, mfe
        return [i for i in g]

    def energy_generator( seql, strucl, model ):
        if strucl == None:
            e_call = self.__getattribute__( 'energy_' + model.lower() + '_mfe' ) # mfe
        else:
            e_call = self.__getattribute__( 'energy_' + model.lower() )

        for seq,struc in zip(seql, strucl or [None]*len(seql)):
            yield e_call( seq, struc )

        
        return __energy_dispatch( proc=energy_proc, seq_l = kargs['sequence_list'], struc_l = kargs['strucure_list'], *args, **kargs )

    
    elif len(args)==2:
        return __energy_dispatch( seq = args[0], struc = args[1] )
    else:
        raise ValueError('Not sure how to use the passed parameters, please see the docstring.')

def energy_dispatch( seq = None, struc = None, seq_l = None, struc_l = None, *args, **kargs ):
    
            

def vienna_energy_setup( *args, **kargs ):
    proc = subprocess.Popen(["RNAeval", "-d0","-P","dna.par",'-logML'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    def energyfunc( sequence, structure ):
        input_str = "{seq}\n{struc}\n".format( seq=sequence, struc=structure )
        proc.write(input_str)
        proc.stdin.flush()

        l = self.vienna_proc.stdout.readline()
        while not (l.startswith(" energy") or (l.startswith( structure ) and len(l)>len(structure))):
            l = self.vienna_proc.stdout.readline()

        if l.startswith(" energy"):
            energy = float(l[l.rfind('=')+1:])
        else:
            energy = float(l[l.rfind('(')+1:l.rfind(')')])
        return round(energy,2)
    return energyfunc

def nupack_energy_setup( *args, **kargs ):
    def energyfunc( sequence, structure ):
        p = subprocess.Popen(["energy", "-material", "dna"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        input_str = "{seq}\n{struc}\n".format( seq=sequence, struc=structure )
        
        stdout_chars = p.communicate(input_str)[0]
        flag = False
        energy = None
        # value we want is on the line after "% Energy"
        result = stdout_chars.split("\n")
        for l in result:
            if flag:
                energy = float( l )
                flag = False
            if l.startswith("% Energy"):
                flag = True
        return round(energy,2)
    return energyfunc

def multistrand_energy_setup( *args, **kargs ):
    import multistrand.objects
    import multistrand.system
    def energyfunc( sequence, structure ):
        c = multistrand.objects.Complex( sequence=sequence, structure=structure )
        energy = multistrand.system.energy([c])
        return round(energy[0],2)    
    return energyfunc
