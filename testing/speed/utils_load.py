import numpy
import os,os.path
import re

#import matplotlib
#matplotlib.use('macosx')
#import matplotlib.pylab as pyp

from speed_general import Length_Result

def initval( a, b ):
  return tuple( [i-j for i,j in zip(a,b)] )

def load_numbered( prefix, seq_range = False, seq_count = False ):
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
      res = Length_Result({'Kinfold':data[0][0],
                           'Multistrand':data[1][0],
                           'Kinfold_init':initval(data[0][1],data[0][0]),
                           'Multistrand_init':initval(data[1][1],data[1][0]),
                           'length':seqlen,
                           'maxtime': seqlen <= 40 and 5000.0 or 1000.0})
    elif isinstance(data, dict):
      res = Length_Result(data)
    elif isinstance(data,Length_Result):
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

def data_series( data, name, avg=False, interval = None ):
  res = []
  keys = data.keys()
  keys.sort()
  if avg or interval != None:
    res = numpy.ndarray( shape=(len(keys),avg and 2 or 3) )
    for k in range(len(keys)):
      scale = max(data[keys[k]], key=lambda x:x==None and 0 or 1)['maxtime'] * 100.0
      if avg:
        mn = numpy.mean( [i[name] for i in data[keys[k]] if i != None] )
        mn /= scale
        # scale for total time units simulated.
        res[k,:] = [keys[k],mn]
      else:
        vals = [i[name] / scale for i in data[keys[k]] if i != None]
        vals.sort()
        res[k,:] = [keys[k], vals[int(interval * len(vals))], vals[int((1.0-interval) * len(vals))]]
    return res
  else:
    # not avg or interval
    for k in keys:
      for l in data[k]:
        if l is not None:
          l_scale = l[name] / (l['maxtime']*100.0)
          res.append( (k, l_scale))
    return numpy.array(res)


def close_check_kin( data ):
  def not_close(i, idx):
    return [i['Kinfold'],i['Kinfold_user_sys'],idx] if i != None and not i['Kinfold_close'] else []

  def checkall( k, v ):
    return [i + [k] for i in v if len(i) > 0]

  res = []
  for k,v in data.iteritems():
    res += checkall(k, map( not_close, data[k], range(len(data[k]))))
  return res

def close_check_ms( data ):
  def not_close(i, idx):
    return [i['Multistrand'],i['Multistrand_user_sys'],idx] if i != None and not i['Multistrand_close'] else []

  def checkall( k, v ):
    return [i + [k] for i in v if len(i) > 0]

  res = []
  for k,v in data.iteritems():
    res += checkall(k, map( not_close, data[k], range(len(data[k]))))
  return res
