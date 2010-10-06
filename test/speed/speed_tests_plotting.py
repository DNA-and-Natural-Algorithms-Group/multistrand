# import matplotlib
# matplotlib.use('macosx')

# import matplotlib.pylab as pyp

from utils_plot import *
from utils_load import *

import cPickle

#utils_plot already grabs the backends for us.

# @figuresetup(10, title = "Kinfold vs Multistrand Comparison, Weak MFE sequences")
# @prepdata( None )
# def lowmfe_from_dataset( data = None):
#     pyp.scatter( data[0][0], data[0][1], label = "Weak MFE [Kinfold]", color='r' )
#     pyp.scatter( data[1][0], data[1][1], label = "Weak MFE [Multistrand]", color='m' )

    
# @figuresetup(11, title = "Kinfold vs Multistrand Comparison, Strong MFE sequences")
# def highmfe_from_dataset():
#     pyp.scatter( highmfe_xy[0], highmfe_xy[1], label = "Strong MFE [Kinfold]", color='b' )
#     pyp.scatter( highmfe_xy[0], highmfe_xy[1], label = "Strong MFE [Multistrand]", color='c' )



# @figuresetup(12, title = "Kinfold vs Multistrand Comparison, Strong vs Weak MFE sequences")
# @prepdata( None )
# def bothmfe_from_dataset( data = None):
#     pyp.scatter( data[0][0], data[0][1], label = "Weak MFE [Kinfold]", color='r' )
#     pyp.scatter( data[1][0], data[1][1], label = "Weak MFE [Multistrand]", color='m' )

#     pyp.scatter( highmfe_xy[0], highmfe_xy[1], label = "Strong MFE [Kinfold]", color='b' )
#     pyp.scatter( highmfe_xy[0], highmfe_xy[1], label = "Strong MFE [Multistrand]", color='c' )

# @figuresetup(13, title = "Kinfold vs Multistrand Comparison, 4-way Branch Migration")
# @prepdata( range(70,120,5)*5, lambda x: x * random.random() / 1000 )
# def fourway_bm( data = None ):
#     pyp.scatter( data[0][0], data[0][1], label = "4-way BM [Kinfold]", color = 'r' )
#     pyp.scatter( data[1][0], data[1][1], label = "4-way BM [Multistrand]", color = 'g' )


# @figuresetup(14, title = "Kinfold vs Multistrand Comparison, 3-way Branch Migration")
# @prepdata( range(70,120,5)*5, lambda x: x * random.random() / 1000 )
# def threeway_bm( data = None ):
#     pyp.scatter( data[0][0], data[0][1], label = "3-way BM [Kinfold]", color = 'r' )
#     pyp.scatter( data[1][0], data[1][1], label = "3-way BM [Multistrand]", color = 'g' )


# @figuresetup(15, title = "Kinfold vs Multistrand Comparison, Hairpin Folding")
# @prepdata( range(10,80,5)*5, lambda x: x * random.random() / 600 )
# def hairpin_data( data = None ):
#     pyp.scatter( data[0][0], data[0][1], label = "Hairpin [Kinfold]", color = 'r' )
#     pyp.scatter( data[1][0], data[1][1], label = "Hairpin [Multistrand]", color = 'g' )

class FigureObject( object ):
    def __init__(self):
        self.figures = []
        self.namecolors = {'Kinfold':'r','Multistrand':'b'}
        for k,v in type(self).__dict__.iteritems():
            if callable(v) and v.__name__  == 'figuresetup_return':
                self.figures.append( (v.num, k) )
    def __len__(self):
        return len(self.figures)
    def showAll(self):
        for num,fname in self.figures:
            self.__getattribute__(fname)()
        
    def loaddata(self, attrname, directory ):
        res = {}
        if os.path.isfile(attrname + '.dat'):
            f = open(attrname + '.dat','rb')
            res = cPickle.load(f)
            f.close()

        read_count = 0
        for k in res.iterkeys():
            read_count += len([i for i in res[k] if i != None])

        if read_count != len([i for i in os.listdir(directory) if i.endswith('.out')]):
            print("Loading '{0}' from individual files.".format(attrname))
            res = load_numbered(directory, seq_count=100)
        else:
            print("Loaded '{0}' from summary file.".format(attrname))
        return res
        
    def dump_data(self, obj, attrname):
        f = open(attrname + '.dat','wb')
        cPickle.dump( obj, f, protocol=-1 )
        f.close()

class FullyRandomSequences( FigureObject ):
    _unique_id = 0
    
    @property
    def unique_id():
        FullyRandomSequences._unique_id += 1
        return FullyRandomSequences._unique_id - 1
    
    def __init__(self):
        FigureObject.__init__(self)
        self.data_short = self.loaddata( 'data_short', 'length_short')
        self.data_bm = self.loaddata( 'data_bm', 'short_4way')
        # self.data_long = self.loaddata( 'data_long', 'length_longs' )
        # self.data_vlong = self.loaddata( 'data_vlong', 'length_very_longs')
        self.data_full = self.data_short.copy()
        # self.data_full.update( self.data_long )
        # self.data_quitefull = self.data_full.copy()
        # self.data_quitefull.update( self.data_vlong )

    def dump_all_data(self):
        self.dump_data( self.data_short, 'data_short')
        self.dump_data( self.data_vlong, 'data_vlong')
        self.dump_data( self.data_long, 'data_long')
        
    @figuresetup( 1, yscale='linear', xlim=(0,205) )
    def plain_full_no_log(self):
        plotdata( self=self, data=self.data_full, label="{0} (Random Sequences)" )
        plotdata( self=self, data=self.data_full, avg=True )

    @figuresetup(2,  yscale='log', xlim=(0,105) )
    def plain_full_log(self):
        plotdata( self=self, data=self.data_full, label="{0} (Random Sequences)" )
        plotdata( self=self, data=self.data_full, avg=True )
        
    @figuresetup(3,  yscale='log', xlim=(0,205) )
    def shading_full(self):
        plotdata( self=self, data=self.data_full, label="{0} (Random Sequences)", alpha=.8 )
        plotdata( self=self, data=self.data_full, interval=.4, alpha=.3 )
        plotdata( self=self, data=self.data_full, interval=.05 )
        plotdata( self=self, data=self.data_full, interval=.2, alpha=.4)

        
    @figuresetup(4,  yscale='log', xlim=(0,205) )
    def user_sys_comparison(self):
        plotdata( self=self, data=self.data_full, label="{0} (Real)" )
        plotdata( self=self, data=self.data_full, avg=True )
        plotdata( self=self, data=self.data_full, label="{0} (User+Sys)", user=True, offset=1.0, colors={'Kinfold':'m','Multistrand':'c'} )
        plotdata( self=self, data=self.data_full, avg=True, user=True, offset=1.0, colors={'Kinfold':'m','Multistrand':'c'} )

    @figuresetup(5,  yscale='log', xlim=(0,105) )
    def user_sys_comparison_small(self):
        plotdata( self=self, data=self.data_full, xrange=(0,105), s=.2, label="{0} (Real)" )
        plotdata( self=self, data=self.data_full, xrange=(0,105), avg=True )
        plotdata( self=self, data=self.data_full, xrange=(0,105), s=.2, label="{0} (User+Sys)", user=True, offset=1.0, colors={'Kinfold':'m','Multistrand':'c'} )
        plotdata( self=self, data=self.data_full, xrange=(0,105), avg=True, user=True, offset=1.0, colors={'Kinfold':'m','Multistrand':'c'} )

    @figuresetup(6,  yscale='log', xlim=(0,105) )
    def user_sys_comparison_areas(self):
        plotdata( self=self, data=self.data_full, xrange=(0,105), interval=.2)
        plotdata( self=self, data=self.data_full, xrange=(0,105), interval=.2,  user=True, alpha=.3)
        #        plotdata( self=self, data=self.data_full, xrange=(0,105), s=.2, label="{0} (Real)", alpha=.5 )
        #       plotdata( self=self, data=self.data_full, xrange=(0,105), s=.2, label="{0} (User+Sys)", user=True, offset=1.0, colors={'Kinfold':'g','Multistrand':'g'}, alpha = 0.5 )
        
    @figuresetup(7,  yscale='log', xlim=(0,305) )
    def longer_full(self):
        plotdata( self=self, data=self.data_quitefull, label="{0} (Random Sequences)", alpha=.7 )
        plotdata( self=self, data=self.data_quitefull, avg=True, alpha=.7 )
        plotdata( self=self, data=self.data_quitefull, avg=True, user=True, colors={'Kinfold':'m','Multistrand':'c'} )

    @figuresetup(8, yscale='log', xlim=(18,102) )
    def alternate_colors_log( self, primary=False ):
        if primary:
            plotdata( self=self, data=self.data_full,
                      label="{0} (No cheats)")
            plotdata( self=self, data=self.data_full, avg=True )
        else:
            plotdata( self=self, data=self.data_full,
                      label="{0} (Kinfold cheats)",offset=1.0,
                      colors={'Kinfold':'m','Multistrand':'c'} )
            plotdata( self=self, data=self.data_full, avg=True,
                      offset=1.0, colors={'Kinfold':'m','Multistrand':'c'} )

    @figuresetup(9,yscale='log',xlim=(40,202))
    def fourway_data( self ):
        plotdata(self=self, data=self.data_bm, label="{0} (4-way)", alpha=.8)
        plotdata(self=self, data=self.data_bm, avg=True)
        plotdata(self=self, data=self.data_bm, interval=.2, alpha=.4)
        plotdata(self=self, data=self.data_bm, interval=.05)
        plotdata(self=self, data=self.data_bm, interval=.4, alpha=.3)

if __name__ == '__main__':
    Figures_FullRandom = FullyRandomSequences()



