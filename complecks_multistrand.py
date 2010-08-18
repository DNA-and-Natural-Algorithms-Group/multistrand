##
#
################################################################################
#                                                                              #
#                                                                              #
#  ____                                ___                  __                 #
# /\  _`\                             /\_ \                /\ \                #
# \ \ \/\_\    ___     ___ ___   _____\//\ \      __    ___\ \ \/'\     ____   #
#  \ \ \/_/_  / __`\ /' __` __`\/\ '__`\\ \ \   /'__`\ /'___\ \ , <    /',__\  #
#   \ \ \L\ \/\ \L\ \/\ \/\ \/\ \ \ \L\ \\_\ \_/\  __//\ \__/\ \ \\`\ /\__, `\ #
#    \ \____/\ \____/\ \_\ \_\ \_\ \ ,__//\____\ \____\ \____\\ \_\ \_\/\____/ #
#     \/___/  \/___/  \/_/\/_/\/_/\ \ \/ \/____/\/____/\/____/ \/_/\/_/\/___/  #
#                                  \ \_\                                       #
#                                   \/_/                                       #
#                                                                              #
#                                                                              #
#                                                                              #
################################################################################

#
# BAM! No imports = good! 
#

class Complecks(object):
    ''' This class is designed to store the information about a single complex from Multistrands output.
    
    A complex consists of a number of connected strands along with dot-paren structure information'''
    def __init__(self, stra = [], sequ = '', struc= '', ener = 0.0, snames = dict()):
        super(Complecks,self).__init__()

        if stra is []:
            raise PendingDeprecationWarning("WARNING: A Complecks was instantiated with no strands in it. This shouldn't be necessary, I think!")
        self.strands = stra
        self.sequence = sequ
        self.structure = struc
        self.energy = ener
        self.strandnames = snames #Associates strand ids with strand names
        self.seqlist = []
        self.struclist = []
        self.regioninfo = []
        
        self._make_canonical()
        self._make_regioninfo()

    def __str__(self):
        res = "Multistrand Complex: \n" + \
        "  Strand List: {0}\n" +\
        "  Sequence   : {1}\n" +\
        "  Structure  : {2}\n" + \
        "  Energy     : {3}\n"
        
        return res.format( str(self.strands),
                    str(self.sequence),
                    str(self.structure),
                    str(self.energy))
    

    def __repr__(self):
        return "Complecks({0},{1},{2},{3})".format( repr(self.strands),
                                                    repr(self.sequence),
                                                    repr(self.structure),
                                                    repr(self.energy))
    
    def _make_canonical( self ):
        """ Re-orders the strands according to the 'canonical' ordering.

        Specifically, the 'canonical' ordering is the circular permutation of the strands
        that has the smallest strand ID (a number) at the beginning.

        As a side effect, this function also breaks up the
        sequences/structures into strand-level sequences/structures,
        accessible via <thisobject>.seqlist, <thisobject>.struclist,
        each of which is a list and in the same order as
        <thisobject>.strands.  (Note that <thisobject>.sequence
        possibly is rearranged after this step, as well.)
        """
        idx = self.strands.index( min( self.strands, key=int ))
        # important to use key=int, otherwise we get a lexicographic ordering
        # on the strings: e.g. '1'<'10'<'2'. Whoops.
        if idx is 0:
            self.seqlist = self.sequence.split("+")
            self.struclist = self.structure.split("+")
            return
        #else we need to shift circular permutations to get the
        #canonical ordering
        self.pivot( idx )
    def join(self, other, target_strand_order ):
        if self.strands[0] == target_strand_order[0]:
            begin = self
            insert = other
        else:
            begin = other
            insert = self
        idx = target_strand_order.index(insert.strands[0])
        # where is the inserted strand going in the target ordering?
        # idx is now the location where we can pivot begin to make the insertion a concatenation.
        # e.g. begin = 0,1,2,3,  insert = 7,8
        #  target = 0,1,2,7,8,3
        #  idx = 3 (7,8 starts at index 3 of target's strands)
        #  begin.pivot(idx) is now:
        #     3,0,1,2
        # so insert + begin is:
        #     7,8,3,0,1,2
        # and then the canonical ordering is:
        #     0,1,2,7,8,3. Done!
        begin.pivot(idx)
        return insert + begin
    
    def __add__(self, other ):
        tempstrands = self.strands + other.strands
        tempseq = "+".join([self.sequence,other.sequence])
        tempstruc = "+".join([self.structure, other.structure])
        tempenergy = str(float(self.energy) + float(other.energy))
        tempstrandnames = self.strandnames.copy()
        tempstrandnames.update( other.strandnames )
        return Complecks( tempstrands, tempseq, tempstruc, tempenergy, tempstrandnames )

    def make_abs_pairtable(self,structure=None,unpairedval = -1):
        """ returns a list of integers of at most len(structure), with
        the pairing info. Also returns a list of tuples of pairing
        indexes.

        NOTE: must always receive sane python objects, no Cstate stuff."""
        if structure == None:
            structure = "".join(self.struclist)
        paren_stack = []
        pairs = [-1]*len(structure)
        pairs[:] = structure[:]
        pair_tuples = []
        try:
            for i in range(len(pairs)):
                if pairs[i] == '(':
                    paren_stack.append(i)
                elif pairs[i] == ')':
                    idx = paren_stack.pop(-1) # pop from right
                    pairs[idx] = i
                    pairs[i] = idx
                    pair_tuples.append((idx,i))
                else:
                    pairs[i] = unpairedval
        except IndexError(e):
            raise IndexError("Unbalanced parentheses in structure [{0}]: (too few '(' at position {1}.".format(structure, i))
        if len(paren_stack) > 0:
            raise IndexError("Unbalanced parentheses in structure [{0}]: (too few ')' by end of structure, needed {1} more.".format(structure, len(paren_stack)))
            
        return pairs,pair_tuples




    def pivot(self, index ):
        if index is 0:
            return
        self.strands = self.strands[index:] + self.strands[:index]
        n = len(self.strands)
        
#        from IPython.Debugger import Pdb; Pdb().set_trace()
        # See appendix for an amusing digression on an earlier
        # solution to this segment.
        
        pivot_chr = len( self.sequence.rsplit("+",n-index)[0])
          # character at this position should be a +.
          # note that once again we fail if index=0.
        self.sequence = self.sequence[pivot_chr+1:] + "+" + \
                        self.sequence[:pivot_chr]
        t_struc       = self.structure[pivot_chr+1:] + "+" + \
                        self.structure[:pivot_chr]

        self.seqlist = self.sequence.split("+")

        # End setup of the flat lists in a rotated order.
        
        def magic_parens( struc_it ):
            """ This function abuses boolean operators.

            No, really, it does. 
            """
            tot = 0
            for i in struc_it:
                tot = tot + (i == '(' and 1 or \
                             i == ')' and -1 or \
                             0)
                res = tot <  0 and (i == ')' and '(') or \
                      tot <= 0 and (i == '(' and ')') or \
                                    i
                if tot<0:
                    tot=0
                yield res
        def paren_reverse( struc ):
            for i in struc:
                yield (i == ')' and '(') or \
                      (i == '(' and ')') or \
                      i
        def either( a,b,c):
            for i,j,k in zip(a,b,c):
                yield i != k and i or j != k and j or k
        temp_forward = "".join(magic_parens(t_struc))
        
        temp_rev = "".join(magic_parens(paren_reverse(reversed(t_struc))))
        temp_reverse = "".join( paren_reverse(reversed(temp_rev)))
        
        self.structure = "".join(either( temp_forward, temp_reverse, t_struc ))
        self.struclist = self.structure.split("+")
    def get_pair_list(self, **kargs):
        """ Returns a list of base pairing info, indexed across the entire complex, including + separators. Use a flag for strand-relative indices.

        For a 'prettier' version of the return value, try .sort() to at least put
        it in some sort of interesting order.

        Details:
          returns a list of (i,j) pairs, where i,j are indices into self.sequence
          such that base i is bound to base j in self.structure.

        Alternate form:
          x.get_pair_list( strand_relative=True )
        Details:
          returns a list of (u,v) pairs, where
              u ~ (strand_i, base_i) and v ~ (strand_j, base_j)
              such that base_i on strand_i is bound to base_j on strand_j
              if self.structure.

              the strand_i are integer indices into self.strands, and
              the base_i are integer indices into self.seqlist[strand_i]
        """
        stack = []
        pairs = []
        if kargs.has_key("strand_relative") and kargs["strand_relative"] is True:
            pair_item = \
               lambda *args: [(\
                               (args[0], \
                                args[1] \
                                - int(sum(map(len,self.seqlist[:args[0]])))\
                                - args[0] ),\
                               (args[2], \
                                args[3] \
                                - int(sum(map(len,self.seqlist[:args[2]])))\
                                - args[2] )\
                              )]
            s_idx = lambda base_idx: \
                               int(sum(map(lambda x: x=='+' and 1 or 0, \
                                       self.sequence[:base_idx])))
        else:
            pair_item = lambda *args: [(args[1],args[3])]
            s_idx = lambda base_idx: 0
            
        """ s_idx: Returns the strand index associated with a particular base index into the complex's sequence (including +'s).
            """

        for i in range(len(self.sequence)):
            ch = self.structure[i]
            if ch == '(':
                stack = stack + [i]
            elif ch == ')':
                item = stack.pop()
                pairs = pair_item( s_idx(item), item, \
                                   s_idx(i), i) \
                      + pairs
        return pairs
    def _make_regioninfo(self):
        idx = 0
        self.regioninfo = [None] * len(self.strands)
        for i in range(len(self.strands)):
            length = len(self.seqlist[i])
            self.regioninfo[i] = (idx,idx+length-1)
            idx = idx + length
    def __len__(self):
    #Overloaded len() function
        return len(self.sequence)
    def __eq__(a, b):
    #Overloaded equivalence operator. If any of the member variables are not equal, the Complecks are not equal
        if (a.strands != b.strands):
            return False
        elif (a.sequence != b.sequence):
            return False
        elif (a.structure != b.structure):
            return False
        elif (a.energy != b.energy):
            class MultistrandError(Exception): pass
            raise MultistrandError("ERROR: Computed energies for identical structures were different! Either there's a floating point precision issue or Multistrand is seriously borked. Please let JS know!")
            # JS: if this happens, there are more serious
            #     issues. Flagging with an exception now - it's
            #     possible it will be triggered due to floating point
            #     arithmetic, if that is the case we should stop using
            #     != for this case.
        else:         
            return True 
    def __ne__(a, b):
            if (a == b): return False
            else: return True
    def combine(self, com):
        raise NotImplementedError("ERROR: This method is obsolete and should be rewritten if it's needed. In the current form it's trying to return a new object, that could be easily implemented by e.g.\n return Complecks(self.strands+com.strands, self.sequence + '+' + com.sequence, ...)\n This would probably be good to just have as the default action for the + operator instead.")
        temp = Complecks()
        temp.strands.extend(self.strands)
        temp.strands.extend(com.strands)
        temp.sequence = self.sequence + '+' + com.sequence
        temp.structure = self.structure + '+' + com.sequence
        return temp
    def clear(self):
        raise PendingDeprecationWarning("WARNING: Shouldn't need to clear data out of a Complecks anymore.")
        self.strands = []
        self.sequence = ''
        self.structure = ''
        self.energy = 0.0
    def to_Cstate(self):
        raise PendingDeprecationWarning("This function will be removed shortly! Please don't use it - it breaks the information encapsulation ideas. If you need a list of basepairs, try the new function get_pair_list: \n help(Complecks.get_pair_list)")
        #####returns the information in the complex stored as a C++ state struct
        ###to_return = naview.Cstate()
        ###to_return = make_pairtable( self.structure, self.sequence)
        ###return to_return
        return None
    def write(self):
        raise PendingDeprecationWarning("WARNING: Complecks.write() is obsolete. To print it to stdout, try: \n print myComplecksObject\n To print it to another stream, just send str(myComplecksObject) to that stream!")
        #A method that prints the contents of the Complecks to stdout.
        print self.strands;print self.sequence; print self.structure; print "Energy: %f"%self.energy
    def output(self, outs):
        raise PendingDeprecationWarning("WARNING: Complecks.output() is obsolete. To print it to stdout, try: \n print myComplecksObject\n To print it to another stream, just send str(myComplecksObject) to that stream!")
        #A method analogous to write that outputs to the filestream outs
        outs.write(str(self.strands));outs.write('\n' + self.sequence+'\n');outs.write (self.structure + '\n'); outs.write("Energy: %f \n"%self.energy)






# APPENDIX 1:
#   a digression on the pivot algorithm, also known as an interesting
#    algorithm that ended up leading to a much simpler one, but only
#    after pondering this one for a while.
#
# V1 algorithm for setting up the rotated flat lists:
#
# self.seqlist = self.sequence.rsplit("+",n-index)[1:] + \
#                self.sequence.split("+",index)[:index]
# t_strucs = self.structure.rsplit("+",n-index)[1:] + \
#            self.structure.split("+",index)[:index]
#        
# self.sequence = "+".join(self.seqlist)
#
# t_struc = "+".join(t_strucs)
#
# Note: This fails if index = 0.
#  Why? Because we use rsplit and split to chop up the strands
#       in what should be a symmetrical way, but fails at an edge
#       case due to hard coding [1:]...
#       idx = 3   =>
#       A B C D   => D A B C
#       rsplit(4-idx)[1:] => D     # 1 splitting performed.
#        split(3)[:3]     => A B C # 3 splits.
#       # of segments  = 2 + 4  (ABC,D for rsplit and A,B,C,D for split)
#       # of splits seleted = (2-1) + (4-1) = 4
#
#
#       idx = 0   =>
#       A B C D   => A B C D
#       rsplit(4-idx)[1:] => B C D # 4 splits requested, 3 performed
#        split(0)[:0]     =>       # 0 splits requested.
#       # of segments  = 4 + 1  (A,B,C,D for rsplit and ABCD for split)
#       # of splits seleted = (4-1) + (1-1) = 3  #doh!
#
#  I see you note that we didn't have to use rsplit(...)[1:],
#  we could make it /slightly/ more complex to cover this edge
#  case. Well, yes. It'd become this:
#        rsplit(...)[index and 1 or (not index) and 0:]
#  Still, it's faster to just special case the index anyways.
#
#  You may note that this exact issue shows up again in the new algorithm:
#   because there are 4 splits requested and only 3 performed, the pivot_char
#   position (which should be at -1) is actually the same as if index=1.
#   The same special casing fix can solve this, but as that's also a fencepost
#   issue (since there was no + in the pivot_char spot!) and it's a trivial
#   fix to special case, we do that.
#

#
# Footnote: http://www.network-science.de/ascii/
#     font: larry3d
#
