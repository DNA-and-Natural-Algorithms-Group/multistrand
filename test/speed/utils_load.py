import numpy
import os,os.path
import re

import cPickle

from speed_general import LengthResult
Length_Result = LengthResult

class SpeedTestData(object):
  def __init__(self, data_name, directory ):
    self.data_name = data_name
    if isinstance(directory,list):
      self.directories = directory
    else:
      self.directories = [directory]
    res = {}
    if os.path.isfile(self.data_name + '.dat'):
      f = open(self.data_name + '.dat','rb')
      res = cPickle.load(f)
      f.close()

    read_count = 0
    for k in res.iterkeys():
      read_count += len([i for i in res[k] if i != None])

    directory_count = 0
    for search_dir in self.directories:
      directory_count += len([i for i in os.listdir(search_dir) if i.endswith('.out')])

    if directory_count != read_count:
      res = {}
      print("Loading '{0}' from individual files.".format(self.data_name))
      for load_dir in self.directories:
        res.update(self._load_numbered(load_dir, seq_count=100))
    else:
        print("Loaded '{0}' from summary file.".format(self.data_name))

    self.raw_data = res
    self._prep_data()
    
  def to_file(self):
    f = open(self.data_name + '.dat','wb')
    cPickle.dump( self.raw_data, f, protocol=-1 )
    f.close()

  def _load_numbered( self, prefix, seq_range = False, seq_count = False ):
    c = re.compile(r'len_(\d*)_sequence_(\d*)[.]out$')
    _,_,files = os.walk(prefix).next()
    results = {}
    for f in files:
      m = c.match( f )
      try:
        seqlen = int(m.group(1))
        seqcnt = int(m.group(2))
      except AttributeError:
        continue
    
      if seq_range and seqlen not in seq_range or \
             not (seq_count and 0 <= seqcnt < seq_count):
        continue


      input_file = open( os.path.join( prefix, f ) )
      data = eval( input_file.readline() )
      input_file.close()
      if isinstance(data, tuple ):
        def initval( a, b ):
          return tuple( [i-j for i,j in zip(a,b)] )

        res = LengthResult({'Kinfold':data[0][0],
                             'Multistrand':data[1][0],
                             'Kinfold_init':initval(data[0][1],data[0][0]),
                             'Multistrand_init':initval(data[1][1],data[1][0]),
                             'length':seqlen,
                             'maxtime': seqlen <= 40 and 5000.0 or 1000.0})
      elif isinstance(data, dict):
        res = LengthResult(data)
      elif isinstance(data,LengthResult):
        res = data
      else:
        raise ValueError("Did not get a tuple or dict input from file [{0}]".format(f))

      if seqlen not in results and seq_count:
        results[seqlen] = [None]*seq_count
      elif seqlen not in results:
        results[seqlen] = []
      if seqcnt >= len(results[seqlen]):
        results[seqlen] = results[seqlen] + [None]*(len(results[seqlen])-seqcnt-1) + [res]
      else:
        results[seqlen][seqcnt] = res
    return results

  @property
  def names( self ):
    return ['Kinfold','Multistrand']

  def normalize(self, base_data=None, normalizer = None, base_data_key=None ):
    if normalizer == None:
      normalizer = Normalizer( base_data=base_data, sampled_data=self, keyflag=base_data_key )
    poly_normalize = lambda x: x * normalizer.linear + normalizer.constant
    self.data_kinfold[0:4,:,:] = poly_normalize(self.data_kinfold[0:4,:,:])
    self.data_kinfold_init[0:4,:,:] = poly_normalize(self.data_kinfold_init[0:4,:,:])
    self.data_multistrand[0:4,:,:] = poly_normalize(self.data_multistrand[0:4,:,:])
    self.data_multistrand_init[0:4,:,:] = poly_normalize(self.data_multistrand_init[0:4,:,:])
    self.data_full[0:16,:,:] = poly_normalize(self.data_full[0:16,:,:])

  def data_series( self, name, avg =  False, interval = None, scaling=None):
    if name == 'Kinfold':
      data_used = self.data_kinfold
    else:
      data_used = self.data_multistrand

    if avg:
      avged_data = numpy.mean( data_used[0,:,:] , 1 )
      xvalues = numpy.max( data_used[4,:,:], 1 )
      for i in range(len(avged_data)):
        min_count = data_used[6,i,:][data_used[6,i,:].nonzero()].min()
        if min_count < 100.0:
          avged_data[i] = 100.0 * avged_data[i] / min_count

      return numpy.vstack( [xvalues, avged_data]).T
    
    if interval != None:
      sorted_data = numpy.sort( data_used[0,:,:], 1 )
      # Note that since we have selected only the timing data
      # 'column', there are only two axes; length and sequence. We
      # want to sort over the sequence results. This copies rather than
      # modifies the internal data.
      res = numpy.zeros( shape=(sorted_data.shape[0],3) )
      # array of x, low_y, high_y for use by the matplotlib.pylab.fill_between
      for i in range(sorted_data.shape[0]):
        res[i,0] = numpy.max( data_used[4,i,:], 0 )
        count = len([j for j in sorted_data[i,:] if not numpy.allclose(j,0.0)]) 
        res[i,1] = sorted_data[i,::-1][int(round( count * interval))]
        res[i,2] = sorted_data[i,::-1][int(round( count * (1.0 -interval)))]
      return res

    # not avg or interval
    xvals = data_used[4,:,:].flatten()
    yvals = data_used[0,:,:].flatten()
    return numpy.vstack([xvals[xvals.nonzero()],yvals[yvals.nonzero()]]).T

  def _prep_data(self):
    """ Separate the dictionary-keyed data into actual numpy.ndarrays to use for plotting and statistics.

    These data attributes are named SpeedTestData.data_[x], where [x] is replaced by
    'kinfold', 'multistrand', 'kinfold_init', 'multistrand_init'. The full dataset is then
    SpeedTestData.data_full.

    The sets have three dimensions, the first is the 'column' type,
    second is the length key index, third is the sequence key
    index. For the individual sets, columns 0-3 have the real,
    user+sys, user, sys times, column 4, 5, 6 are the sequence raw
    length, max simulation time, number of sequences tested, and columns 7-9 are the load
    averages at the end of the run. For the combined set, columns 0-3,
    4-7, 8-11, 12-15 are the times for kinfold, kinfold_init,
    multistrand, multistrand_init in the same format as above, and
    columns 16-21 are the same as 4-9 in the individual format:
    length, time, sequence count,  load averages.

    Note that any timing data has had the max time and # of
    simulations divided out. One can get the actual reported time by
    examining the raw dictionary entries, or by using the maxtime
    field and 100 trajectory counts.
    """
    res = []
    keys = self.raw_data.keys()
    keys.sort()
    l = len(keys)
    # l is the number of x locations we have sampled.

    for value in ['kinfold','kinfold_init','multistrand','multistrand_init']:
      self.__setattr__('data_' + value,
                       numpy.zeros( shape=(10,l,100) )
                       )
      d = self.__getattribute__('data_' + value)
      # note that this does not work if the attribute retreived is not mutable.

      for i in range(l):
        for j in range(len(self.raw_data[keys[i]])):
          if self.raw_data[keys[i]][j] != None:
            data_tuple = self.raw_data[keys[i]][j][value.capitalize()]
            d[0:4,i,j] = (max(data_tuple[0],data_tuple[2])+
                          max(data_tuple[1],data_tuple[3]),
                          data_tuple[4],
                          max(data_tuple[0],data_tuple[2]),
                          max(data_tuple[1],data_tuple[3]))
            ## user+sys, real, user, sys.
            d[4,i,j] = self.raw_data[keys[i]][j].get('length',keys[i])
            d[5,i,j] = self.raw_data[keys[i]][j].get('maxtime',1000.0)
            d[6,i,j] = len([k for k in self.raw_data[keys[i]] if k != None])
            d[7:10,i,j] = self.raw_data[keys[i]][j].get('load',(-1.0,-1.0,-1.0))

            d[0:4,i,j] /= (d[5,i,j] * 100.0)

    self.data_full = numpy.zeros( shape=(22,l,100) )
    self.data_full[0:4,:,:] = self.data_kinfold[0:4,:,:]
    self.data_full[4:8,:,:] = self.data_kinfold_init[0:4,:,:]
    self.data_full[8:12,:,:] = self.data_multistrand[0:4,:,:]
    self.data_full[12:16,:,:] = self.data_multistrand_init[0:4,:,:]
    self.data_full[16:22,:,:] = self.data_kinfold[4:10,:,:]

class Normalizer( object ):
  def __init__(self, base_data, sampled_data, keyflag=None):
    if keyflag == None:
      keyflag = 'Kinfold'
    base_kin_data = base_data.data_series('Kinfold',avg=True)
    sample_kin_data = sampled_data.data_series('Kinfold',avg=True)

    base_ms_data = base_data.data_series('Multistrand',avg=True)
    sample_ms_data = sampled_data.data_series('Multistrand',avg=True)

    selection = [i for i in range(base_kin_data.shape[0]) if base_kin_data[i,0] in sample_kin_data[:,0]]
    
    if keyflag == 'Both':
      base_data = numpy.hstack( [base_kin_data[:,1][selection], base_ms_data[:,1][selection]])
      sample_data = numpy.hstack( [sample_kin_data[:,1], sample_ms_data[:,1]])
    elif keyflag.startswith('Kinfold'):
      base_data = base_kin_data[:,1][selection]
      sample_data = sample_kin_data[:,1]
    elif keyflag.startswith('Multistrand'):
      base_data = base_ms_data[:,1][selection]
      sample_data = sample_ms_data[:,1]
      

    self.fitvalue = numpy.polyfit( sample_data, base_data, 1 )

  @property
  def linear(self):
    return self.fitvalue[0]

  @property
  def constant(self):
    return self.fitvalue[1]
    
             
