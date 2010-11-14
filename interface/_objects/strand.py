import warnings

class Strand(object):
  """
  Represents a Multistrand Strand object.
  """
  unique_id = 0

  def __init__(self, *args, **kargs ):
    """ Initializes a new Strand object. """
    if len(args) == 4:
      warnings.warn( DeprecationWarning("Passing the strand ID is deprecated, it's a private matter internal to Multistrand.") )
      self.name,self._sequence,self.domain_list = args[1:4]
    elif len(args) == 3:
      self.name,self._sequence,self.domain_list = args[0:3]
    else:
      try:    self._sequence = kargs['sequence']
      except: self._sequence = ""

      try:    self.name = kargs['name']
      except: self.name = "Automatic_" + str(Strand.unique_id)

      try:    self.domain_list = kargs['domains']
      except: self.domain_list = []

    self.id = Strand.unique_id
    Strand.unique_id += 1

  @property
  def sequence( self ):
    """
    The sequence associated with this strand, computing it as
    necessary from the domains.
    
    Raises:
    ValueError -- Is raised when this attribute is accessed and there was
                  no sequence, either from Domains or by explicit assignment. 

    .. warning::
       While you may set the sequence of a Strand after creation,
       if domains are defined and the given sequence does not match that
       domain specification, result is undefined. In the future this will
       have error checking.

       The following are examples of things you should not do:
       
       >>> d = Domain(name="d", length=5)
       >>> s = Strand(name="s", domains=[d, d])
       >>> s.sequence = "ACTG"
       >>> print len(d.sequence) == d.length
       False
       >>> s.sequence = "AAAAATTTTT"
       >>> print d.sequence
       TTTTT
    """
    if len(self.domain_list) == 0 and len(self._sequence) == 0:
      raise ValueError("ERROR: Strand was queried for a sequence, but it has no domains and no explicitly set sequence.")
    if len(self._sequence) > 0 and len(self.domain_list) == 0:
      return self._sequence
    if len(min( self.domain_list, key = lambda x: len(x.sequence)-int(x.length) ).sequence) < 0:
      raise ValueError("ERROR: Strand was queried for a sequence, but at least one contained domain had a length longer than its current sequence.")

    self._sequence = "".join( [x.sequence for x in self.domain_list] )
    return self._sequence

  @sequence.setter
  def sequence( self, value ):
    """
    Check the passed in value to make sure it's sane.
    """
    try:
      if not all([i.upper() in 'AGCT' for i in value]):
        raise ValueError("At least one of the bases in sequence [{0}] was not a valid base; The first offending character was '{1}', at position {2}.".format( value, value.lstrip('agctAGCT')[0], value.index( value.lstrip('agctAGCT')[0] ) ))
    except TypeError:
      raise ValueError("A strand may only be set to a string of valid bases.")

    if not self.domain_list:
      self._sequence = value
    else:
      # Do some error checking (lengths, domain assumptions, etc.)
      current = 0
      for d in self.domain_list:
        d.sequence = value[current:current + d.length]
        current += d.length

  def __str__(self):
    try:
      temp_seq = self.sequence
    except ValueError:
      temp_seq = self._sequence
    # We want to ignore a possible ValueError here, as it just means
    # the base components haven't been set up yet, but someone may
    # want to actually print a strand at any point. Other types of
    # exceptions should still be passed on though, so we need to be
    # specific here.

    if len(temp_seq) > 30:
      return "\
Strand: {b}       Name:'{0.name}'\n\
      : Sequence [{2}]:{0.sequence}\n\
      : {b}    Domains:{1}".format( self,
                                    [i.name for i in self.domain_list],
                                    len(temp_seq),
                                    b=" "*(len(str(len(temp_seq)))))
    elif len(temp_seq) > 0:
      return "Strand: Name:'{0.name}' Sequence [{2}]:{0.sequence} Domains:{1}".format( self, [i.name for i in self.domain_list], len(temp_seq) )
    else:
      return "Strand: Name:'{0.name}' Domains:{1}".format( self, [i.name for i in self.domain_list] )

  def __add__(self, other ):
    """ Addition of two strands results in a strand composed of each piece. Addition of a strand and a domain adds the domain onto the strand."""
    if isinstance(other, Strand):
      try:
        fullseq = self.sequence + other.sequence
      except ValueError:
        fullseq = ""
      return Strand( name = self.name + '+' + other.name,
                     domains = self.domain_list + other.domain_list,
                     sequence = fullseq )
    try:
      try:
        fullseq = self.sequence + other.sequence
      except ValueError:
        fullseq = ""
      return Strand( name = self.name,
                     domains = self.domain_list + [other],
                     sequence = fullseq)
    except AttributeError:
      return NotImplemented

  def __radd__( self, other ):
    try:
      try:
        fullseq = other.sequence + self.sequence
      except ValueError:
        fullseq = ""
      return Strand( name = self.name,
                     domains =  [other] + self.domain_list,
                     sequence = fullseq )
    except AttributeError:
      return NotImplemented
    #radd is used when the other operand does not support the op.


  @property
  def C(self):
    """
    Returns a Strand object that is complementary to this one.
    """

    return ComplementaryStrand(self)


class ComplementaryStrand( Strand ):
  """
  Represents a complemented strand. This is always defined in
  terms of an original strand, so that it reflects any change in the
  original. It provides the same interfaces as a strand.
  """

  complement = {'G':'C',
                'C':'G',
                'A':'T',
                'T':'A'}

  unique_id = 0

  def __init__( self, complemented_strand ):
    self.id = ComplementaryStrand.unique_id
    ComplementaryStrand.unique_id += 1
    self._strand = complemented_strand
    self._sequence = ""

  @property
  def name( self ):
    if self._strand.name.endswith("*") or \
       self._strand.name.endswith("'"):
      return self._strand.name.rstrip("*'")
    else:
      return self._strand.name + "*"

  @property
  def sequence(self):
    self._sequence = "".join(
        [ComplementaryStrand.complement[i] for i in reversed(self._strand.sequence)])
    # Note that if Strand.sequence was 0 (e.g. since one domain didn't
    # have a sequence yet), the Strand.sequence call will raise an
    # exception, so we never reach these later lines.
    return self._sequence

  @property
  def domain_list(self):
    return [d.C for d in reversed(self._strand.domain_list)]

  @property
  def C(self):
    """
    Returns a Strand object that is complementary to this one.
    """

    return self._strand

