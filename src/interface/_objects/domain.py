from ..utils import generate_sequence
from .strand import Strand


class Domain(object):
    """
    Represents a sequence domain, for use in defining strands and stop conditions for Multistrand.
    """
    _domain_unique_id = 0

    def __init__(self, *args, **kargs ):
        """
        Initialization:
        
        Keyword Arguments:
        name [type=str,required]       -- Name of this domain.
        sequence [type=str,default=""] -- Sequence of this domain, e.g. 'AGGCAGATA'
        length [default=0]             -- Length of this domain. A 0-length domain
                                          may be useful in some rare cases. If
                                          sequence is set, length should ALWAYS be
                                          == len(sequence), but this is not strictly
                                          enforced.
        """
        if len(args)==4 or (len(args)==3 and 'is_complement' in kargs):
            self.id = args[0]
            self.name = str(args[1])
            self.length = int(args[2])
            self._sequence = ""
        else:
            self.id = Domain._domain_unique_id
            Domain._domain_unique_id += 1
            self._sequence = ""
            self.length = 0
            for k,v in kargs.items():
                self.__setattr__( k, v )
            if 'name' not in kargs:
                raise ValueError("Must pass a 'name' keyword argument.")

    @property
    def sequence( self ):
        """
        The sequence associated with this domain.
        """
        return self._sequence

    @sequence.setter
    def sequence( self, value ):
        """
        Check the passed in value to make sure it's sane and convert from unicode. Also automatically set the length.
        """
        value = str(value)  # convert from unicode if needed
        try:
            if not all([i.upper() in 'AGCT' for i in value]):
                raise ValueError("At least one of the bases in sequence [{0}] was not a valid base; The first offending character was '{1}', at position {2}.".format( value, value.lstrip('agctAGCT')[0], value.index( value.lstrip('agctAGCT')[0] ) ))
        except TypeError:
            raise ValueError("A strand may only be set to a string of valid bases.")
        
        self._sequence = value
        if len(self._sequence) > 0:
            self.length = len(self._sequence)

    def gen_sequence( self, *args, **kargs ):
        """ Uses the same parameters as 'multistrand.utils.generate_sequence', but sets the length to the domain's length."""
        if 'n' in kargs:
            del kargs['n']
        self.sequence = generate_sequence(n = self.length, *args, **kargs )

    @property
    def C(self):
        """
        The complementary domain, specifically the one whose bases in
        the standard 5' to 3' ordering are exactly matching base pairs to
        the original.
        
        >>> a = Domain( sequence='AGGACCATT')
        >>> a.sequence
        AGGACCATT
        >>> b = a.C
        >>> b.sequence
        AATCCTCCT
        >>> b.name
        >>> b.C == a
        True
        """
        return ComplementaryDomain( self )
    
    def __add__( self, other ):
        if isinstance(other,Domain):
            return Strand( domains = [self,other] )
        else:
            return NotImplemented
    
    def __str__( self ):
        return ("Domain : {fieldnames[0]:>9}: '{0.name}'\n : {fieldnames[1]:>9}: {0.length}\n".format( self, fieldnames=['Name','Length'] ) + (self.sequence and 
         '       : {fieldname:>9}: {0}\n'.format( self.sequence, fieldname='Seq' )
         or ''))
    
    
    
class ComplementaryDomain( Domain):
    """
    Represents a complemented domain. Note that this is always
    defined in terms of an original domain and does not have the same
    data members, instead providing an interface to the complementary
    members.
    
    """
    complement = {'G':'C',
                  'C':'G',
                  'A':'T',
                  'T':'A'}
      
    def __init__(self, complemented_domain ):
        self.id = complemented_domain.id
        
        self._domain = complemented_domain
    
    @property
    def length(self):
        return self._domain.length
    
    @property
    def name(self):
        if self._domain.name.endswith("*") or self._domain.name.endswith("'"):
            return self._domain.name.rstrip("*'")
        else:
            return self._domain.name + "*"
    
    @property
    def sequence( self ):
        if self._domain.sequence == None:
            raise ValueError
        else:
            return "".join([ComplementaryDomain.complement[i] for i in reversed(self._domain.sequence.upper())])
    
    @sequence.setter
    def sequence( self, value ):
        self._domain.sequence = "".join([ComplementaryDomain.complement[i] for i in reversed(value.upper())])
    
    def gen_sequence( self, *args, **kargs ):
        """ Uses the same parameters as 'multistrand.utils.generate_sequence', but sets the length to the domain's length."""
        self._domain.gen_sequence( *args, **kargs )
    
    @property
    def C(self):
        """
        The complementary domain, specifically the one whose bases in
        the standard 5' to 3' ordering are exactly matching base pairs to
        the original.
        
        >>> a = Domain( sequence='AGGACCATT')
        >>> a.sequence
        AGGACCATT
        >>> b = a.C
        >>> b.sequence
        AATCCTCCT
        >>> b.name
        >>> b.C == a
        True
        """
        return self._domain
