import matplotlib
matplotlib.use('macosx')
import matplotlib.pylab as pyp

import utils_load

def figuresetup( num, xlabel=None, ylabel=None, title = None, yscale = None,xlim=None,ylim=None ):
    def actual_decorator( funct ):
        def figuresetup_return( self=None, show=False, *args, **kargs ):
            f = pyp.figure( num )

            if 'hold' in kargs and not kargs['hold']:
                f.clear()

            if xlabel == None:
                pyp.xlabel("Sequence length (nt)")
            else:
                pyp.xlabel( xlabel )

            if ylabel == None:
                pyp.ylabel("Real Clock Time (RCT) per Simulated Time Unit (STU)")
            else:
                pyp.ylabel( ylabel )

            if title == None:
                pyp.title("Kinfold vs Multistrand comparison on random sequences", position=(0.5,1.05) )
            else:
                pyp.title( title, position=(0.5,1.05) )
                # position it up slightly to make room for scientific notation labels. 
            if yscale != None:
                pyp.yscale( yscale )
            else:
                pyp.yscale('linear')
                pyp.ticklabel_format( style ='sci', scilimits = (0,0), axis='y' )

            if xlim or ylim:
                pass
            else:
                ax = pyp.axes()
                ax.autoscale(enable=True, axis='x',tight=True)
                ax.autoscale(enable=True, axis='y',tight=False)
            
            funct(self=self, *args, **kargs)
            pyp.legend( loc=0)
            if xlim:
                pyp.xlim( xlim )
            if ylim:
                pyp.ylim( ylim )
            else:
                pyp.ylim( ymin = 0.0 )

            if show:
                pyp.show()
        figuresetup_return.num = num
        return figuresetup_return
    return actual_decorator

def findnames( data ):
    """ Retreives the major names for a particular data point loaded from a file. """
    testval = max( data.itervalues().next(), key=lambda x:x == None and -1 or +1 )
    res = []
    for keyname in testval.iterkeys():
        if keyname[0].isupper() and '_' not in keyname:
            res.append(keyname)
    return res

def plotdata( self=None, data=None, label="", namecolors=None, avg = False, interval = None, user=False, offset=None, x_range=False, *args, **kargs ):
    """ Plots the data set, using the label name and optional args as given to
    select the plotting function style.

    data: should be a SpeedTestData object (see utils_load.py)

    label: Must be a string which has exactly one positional specifier, it will be replaced
           with the successive data series names.
    avg: plot the average for each length n.
    interval: for each length n, plot the range [y_low,y_high]
              where y_low  =y[interval*k] and
                    y_high =y[(1.0-interval)*k],
              k = # of sequences at that position, and the y array is sorted.
              NOTE: alpha is really handy for this plot. It's defaulted to .2 for this style.
              """
    plotfunc = (avg and pyp.plot) or \
               (interval and pyp.fill_between) or \
               pyp.scatter
    plotkargs = (avg and [('linewidth',1.0)] ) or \
                (interval and [('alpha',0.2)] ) or \
                 [('s',0.5)]
    for k,v in plotkargs:
        kargs.setdefault(k,v)
    # this prevents the defaults from overwriting any passed in kargs for those names.
    
    names = data.names

    if 'colors' in kargs:
        colordict = kargs['colors']
        del kargs['colors']
    elif 'color' not in kargs:
        if 'colors' in globals():
            colordict = colors
        elif self != None and hasattr(self,'namecolors'):
            colordict = self.namecolors
        else:
            colordict = {}

    if x_range:
        keyrange = [i for i in keyrange if x_range[0] <= i and x_range[1] > i]

    for n in names:
        kargs['label'] = label.format(n)
        kargs['color'] = colordict.get(n,'g')
                
        ds = data.data_series( name=n, avg=avg, interval=interval )
        if offset != None:
            ds[:,0] += offset
        if names.index(n) == 0:
            ds[:,0] += 0.5

        plotfunc( *(zip(*ds) + list(args)), **kargs )
        
    
def prepdata( rang = None, m_f = None):
    def actual_decorator( f ):
        data_x = rang or (range(10,192,2)*5)
        data_x.sort()
        if m_f:
            data_y = [m_f(i) for i in range( len( data_x ))]
        else:
            data_y = [random.random() * (1000.0 / (i+1)) * random.random() / 1000.0 for i in range(len( data_x ))]
        data_kin = (data_x, data_y)
        data_x = rang or (range(10,192,2)*5)
        data_x.sort()
        if m_f:
            data_y = [m_f(i) for i in range( len( data_x ))]
        else:
            data_y = [random.random() * (1000.0 / (i+1)) * random.random() / 1000.0 for i in range(len( data_x ))]
        data_ms = (data_x, data_y)
        def res():
            f(data=(data_kin,data_ms))
        return res
    return actual_decorator
            



