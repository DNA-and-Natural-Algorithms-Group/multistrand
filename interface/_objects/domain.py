from ..utils import generate_sequence

class Domain(object):
  """Represents a sequence domain, for use in defining strands and stop conditions for Multistrand.

  Data Members:
  id [default=automatic]: unique id representing this domain, automatically
                          set as they are created.
  name [required]: Name of this domain.
  sequence [default=None]: Sequence of this domain, e.g. 'AGGCAGATA'
  length [default=0]: Length of this domain. A 0-length domain may be useful
                      in some rare cases. If sequence is set, length should
                      ALWAYS be == len(sequence), but this is not strictly
                      enforced.
  is_complement [deprecated]: Whether this object is a complement of a different
                              domain object. Should no longer be used, see other
                              notations for handling this in Strand and higher
                              level primitives.
  """
  unique_id = 0
  
  def __init__(self, *args, **kargs ):
    if len(args)==4 or (len(args)==3 and 'is_complement' in kargs):
      self.id, self.name, self.length = args
      self.is_complement = (len(args)==4 and args[3]) or kargs['is_complement']
      self.sequence = None
    else:
      self.id = Domain.unique_id
      Domain.unique_id += 1
      self.sequence = None
      self.length = 0
      for k,v in kargs.iteritems():
        self.__setattr__( k, v )
      if 'name' not in kargs:
        raise ValueError("Must pass a 'name' keyword argument.")

  def gen_sequence( self, *args, **kargs ):
    """ Uses the same parameters as 'multistrand.utils.generate_sequence', but sets the length to the domain's length."""
    if 'n' in kargs:
      del kargs['n']
    self.sequence = generate_sequence(n = self.length, *args, **kargs )

  def __str__( self ):
    return ("\
Domain : {fieldnames[0]:>9}: '{0.name}'\n\
       : {fieldnames[1]:>9}: {0.length}\n".format( self, fieldnames=['Name','Length'] ) +
  (self.sequence and
     '       : {fieldname:>9}: {0}\n'.format( self.sequence, fieldname='Seq' )
     or ''))
