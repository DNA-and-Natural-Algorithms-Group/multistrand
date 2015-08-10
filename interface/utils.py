####################################################################
#                                                                  #
#  Copyright (c) 2010-2015 California Institute of Technology.     #
#  Distributed under the MIT License.                              #
#  (See accompanying file LICENSE or copy at                       #
#  http://opensource.org/licenses/MIT)                             #
#                                                                  #
####################################################################
#                                                                  #
#   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)         #
#                                                                  #
####################################################################


def generate_sequence( n, allowed_bases = ['G','C','T','A'], base_probability = None ):
    """ Generate a sequence of N base pairs.

    Bases are chosen from the allowed_bases [any sequence type], and
    according to the probability distribution base_probability - if
    none is specified, uses uniform distribution."""

    import random
    
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

def energy( sequence=False, structure=False,
                 sequence_list=False, structure_list=False, model=None ):
    """ Gets the energy of a particular sequence and structure, or a list of them.

    If passed two non keyword args, they'll be considered the sequence and structure to compute.
    Keywords: 'sequence','structure': for sequence and structure of a single state.
              'sequence_list','structure_list': for sequence and structures of multiple states
              'model':  'Vienna', 'Nupack' or 'Multistrand'
                         (Note that Multistrand indicates the usage of one of the other
                         models, depending on what you initialized multistrand.system with.
                         The energy calculation is done by Multistrand natively.)
    """
    def energy_generator( seql, strucl, model ):
        if strucl == None:
            e_call = __getattribute__( 'energy_' + model.lower() + '_mfe' ) # mfe
        else:
            e_call = __getattribute__( 'energy_' + model.lower() )

        for seq,struc in zip(seql, strucl or [None]*len(seql)):
            yield e_call( seq, struc )

    def energy_vienna( sequence, structure):
        if multiple_proc == False:
            proc = subprocess.Popen(["RNAeval", "-d0","-P","dna.par",'-logML'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        elif multiple_proc == True:
            multiple_proc = subprocess.Popen(["RNAeval", "-d0","-P","dna.par",'-logML'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            proc = multiple_proc
        else:
            proc = multiple_proc

        input_str = "{seq}\n{struc}\n".format( seq=sequence, struc=structure )
        proc.stdin.write(input_str)
        proc.stdin.flush()

        l = proc.stdout.readline()
        while not (l.startswith(" energy") or (l.startswith( structure ) and len(l)>len(structure))):
            l = proc.stdout.readline()

        if l.startswith(" energy"):
            energy = float(l[l.rfind('=')+1:])
        else:
            energy = float(l[l.rfind('(')+1:l.rfind(')')])
        return round(energy,2)

    def energy_vienna_mfe( sequence, structure):
        if multiple_proc == False:
            proc = subprocess.Popen(["RNAfold", "-d0","-P","dna.par"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        elif multiple_proc == True:
            mfe_proc = subprocess.Popen(["RNAfold", "-d0","-P","dna.par"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            proc = mfe_proc
        else:
            proc = mfe_proc

        input_str = "{seq}\n".format( seq=sequence )
        proc.stdin.write(input_str)
        proc.stdin.flush()

        l_seq = proc.stdout.readline()
        l_mfe_struc = proc.stdout.readline()
        # we then call energy because fold does not actually use the
        # -logML term, so just in case we have a multiloop we need to
        # use the logML part.

        return energy_vienna( l_seq, l_mfe_struc.split()[0] )

    def energy_nupack_mfe( sequence, structure):
        pass

    def energy_nupack( sequence, structure):
        pass

    def energy_multistrand( sequence, structure):
        pass



    if sequence_list or structure_list:
        multiple_proc = True
    else:
        multiple_proc = False

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
        g = energy_generator( sequence_list, structure_list, model )
    elif sequence_list:
        g = mfe_generator( sequence_list, None, model )
        # multiple energy, mfe
    return [i for i in g]



def vienna_energy_setup( *args, **kargs ):
    pass

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
