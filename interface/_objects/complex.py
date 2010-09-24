from objects import Strand

class Complex(object):
  """Represents a Multistrand Complex object."""
  unique_id = 0
  
  def __init__(self, *args, **kargs):
    """ Initialize a Complex object.

    For the old style initialization function, see old_init - using
    that type of parameters when creating a Complex will fall back to
    that function.

    New parameters [all keywords]:
    sequence [required]:  Flat sequence to use for this complex.
    structure [required]: Flat structure to use for this complex.
    name [default=automatic]: Name to use for this complex. Defaults to 'automatic' + a unique integer.
    boltzmann_sample [default=False]: Whether we should boltzmann sample this complex.

    You must include both of the required keyword arguments to create a Complex with the new style init. 
    """
    if len( args ) == 4 or len( args ) == 5:
      self.old_init( *args, **kargs )
    elif 'sequence' in kargs and 'structure' in kargs:
      self.id = Complex.unique_id
      self.name = kargs.get('name') or "automatic" + str(Complex.unique_id)
      self.strand_list = [Strand(sequence=i) for i in kargs['sequence'].split("+")]
      self._fixed_structure = kargs['structure']
      self.boltzmann_sample = kargs.get('boltzmann_sample', False)
      self._last_boltzmann_structure = False
      self._boltzmann_sizehint = 1
      self._boltzmann_queue = []
      Complex.unique_id += 1
      
  def __str__( self ):
    return "\
Complex: {fieldnames[0]:>9}: '{0.name}'\n\
       : {fieldnames[1]:>9}: {0.sequence}\n\
       : {fieldnames[2]:>9}: {0.structure}\n\
       : {fieldnames[3]:>9}: {1}\n\
       : {fieldnames[4]:>9}: {0.boltzmann_sample}".format(
       self,
       [i.name for i in self.strand_list],
       fieldnames =('Name','Sequence','Structure','Strands','Boltzmann') )

  def old_init( self, id, name, strand_list, structure, boltzmann_sample = False):
    """ Old style init function, uses fixed argument list.

    id:  unique identifier for this complex [Not currently used]
    name:  name to refer to this complex by. [Vaguely used]
    strand_list: list of multistrand.object.Strand objects to retreive the sequence from
    structure: flat dot-paren structure of this scomplex
    boltzmann_sample [default=False]: whether we should call a sampler when retreiving a structure from this complex.
    """
    self.id = id
    # what is this id?
    self.name = name
    # what is this name?
    self.strand_list = strand_list
    self._fixed_structure = structure
    self._last_boltzmann_structure = False
    self._boltzmann_sizehint = 1
    self._boltzmann_queue = []
    self.boltzmann_sample = boltzmann_sample

  def get_unique_ids( self ):
    return set([i.id for i in self.strand_list])
  
  def __len__(self):
    """ This may not be a very good definition of length. See notes elsewhere."""
    return len(self.strand_list)

  @property
  def structure(self):
    """ If this complex is set to use boltzmann sampling, this property returns a newly sampled structure. Otherwise it gets the fixed structure for the complex."""
    if self.boltzmann_sample:
      self.generate_boltzmann_structure()
      # puts the generated structure in self._last_boltzmann_structure
      # for use by other properties as well.
      return self._last_boltzmann_structure
    else:
      return self._fixed_structure

  @structure.setter
  def structure(self,value):
    # I include the following due to the 'normal' use cases of Complex
    # involve it immediately being deepcopy'd (e.g. in starting
    # states, etc). So allowing someone to set the fixed structure
    # member is probably not a good idea - we could remove this
    # completely so the structure component is pretty much immutable.
    #
    # In practice, I can see a case where it may be useful to modify it if you're using an interactive mode, but even in those cases it may not do what the user wants, if they're relying on that changing previous parts. For example:
    # c1 = Complex("c1","c1", "......((.......)).....")
    # o = Options()           
    # o.start_state = [c1]
    # ... user runs a simulation and gets a syntax warning about mismatched parens...
    # c1.structure = "......((........))...."
    #
    # # This did NOT actually change o at all! So the warning used
    # # here mentions that, and just in case, returns a new object
    # # anyways!
    #
    import warnings
    warnings.warn("Setting a Complex's structure does not [usually] change existing uses of this Complex, so the object returned is a NEW object to avoid any confusion as to how it may affect previous usages.",SyntaxWarning)
    import copy
    retval = copy.deepcopy(self)
    retval._fixed_structure = value
    return retval
    
  @property
  def fixed_structure(self):
    """ The structure used to create this object. """
    return self._fixed_structure

  @property
  def current_boltzmann_structure(self):
    """ What structure was used for the last sample. Is the value False if no sampling has occurred."""
    return self._last_boltzmann_structure

  @property
  def boltzmann_count(self):
    """ Used to provide a hint as to how many times we're going to need boltzmann samples from this complex, so that get_boltzmann_structure can ask nupack for more in one shot.

    Default: 1"""
    return self._boltzmann_sizehint
  
  @boltzmann_count.setter
  def boltzmann_count(self,value):
    if value >= 1:
      self._boltzmann_sizehint = value
    else:
      self._boltzmann_sizehint = 1
      
  @property
  def sequence(self):
    """ This property is actually the 'flat' sequence for the complex."""
    return "+".join([strand.sequence for strand in self.strand_list])
  
  def generate_boltzmann_structure(self):
    """ Creates a new boltzmann sampled structure for this complex.

    Can be retrieved from other objects via the .current_boltzman_structure property, or
    internal functions can access it directly."""

    if len(self._boltzmann_queue) > 0:
      self._pop_boltzmann()
      return

    import subprocess, tempfile, os
    
    tmp = tempfile.NamedTemporaryFile(delete=False,suffix=".sample")
    tmp.close()
    # we close it here as some OS's have issues opening the same file
    # simultaneous. Will reopen it later to get the data back out.

    # set up the # of structures to grab from the file, max out at
    # ~100 so we don't use too much CPU on this step. JS's testing
    # timed a 100 count at ~ .1s and 10 and 1 counts were almost
    # always around .08s, so at least in this range there's a lot more
    # call overhead than generation time being used.
    if self._boltzmann_sizehint > 100:
      count = 100
    elif self._boltzmann_sizehint >= 1:
      count = self._boltzmann_sizehint
    else:
      count = 1
    # I have changed the below call to just use "sample". It should be
    # in your path, and if not it's probably a 'better' idea to
    # ln -s pathtosample/sample ~/bin/sample
    # rather than do a better path search here. Once it is a part of
    # the NUPACK package we can search in NUPACKHOME. This is preferred as
    # Mac OS X [BSD] has a inbuilt tool called 'sample' as well.
    #
    # Notes:
    # 1) ~/bin may not be in your path. Use the path to your preferred user
    # binaries directory that is in your path.
    # 2) the executable generated by nupack is called 'sample', however there is a
    # standard BSD tool with that name: if you are using OS X, make sure that your path
    # to the nupack 'sample' occurs before /usr/bin or it may not find it correctly.
    #

    p = subprocess.Popen(["sample", "-multi", "-material", "dna", "-count", str(count), tmp.name[:-7]],stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    input_str = "{0}\n{1}\n{2}\n".format( len(self.strand_list),
                                        "\n".join( [i.sequence for i in self.strand_list] ),
                                        " ".join( [str(i+1) for i in range( len( self.strand_list ))])
                                          )
    result = p.communicate(input_str)[0]
    # note we toss the result as it's mostly just spam from the subprocess
    
    f = open(tmp.name, "rt")
    lines = f.readlines()
    f.close()
    os.remove(tmp.name) # was created by us [NamedTemporaryFile] and
                          # used by the sampler, thus we need to clean
                          # it up.
    self._boltzmann_queue = lines[10:]
    if len(self._boltzmann_queue) < 1:
      raise IOError("Did not get any results back from the boltzmann sampler function.")

    self._pop_boltzmann()
    
  def _pop_boltzmann(self):
    """ Pops a structure off our waiting queue, putting it in the correct internal.

    Does not check for the queue being empty: any caller must ensure
    there is something in the queue, or catch the exception raised by
    pop.

    Note also that this implicitly decrements the size hint, so if you
    use more requests than you noted in the size hint, the later ones
    get pulled in much smaller amounts. Theoretically the user should
    poke the complexes and reset the size hint back upwards if they
    need to use more, rather than making this pop smart about dynamic
    resizing of the requested amounts."""
    self._last_boltzmann_structure = self._boltzmann_queue.pop().strip()
    self._boltzmann_sizehint -= 1
