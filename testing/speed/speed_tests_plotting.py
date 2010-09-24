# import matplotlib
# matplotlib.use('macosx')

# import matplotlib.pylab as pyp

from utils_plot import *
from utils_load import *

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


class FullyRandomSequences( FigureObject ):
    _unique_id = 0
    
    @property
    def unique_id():
        FullyRandomSequences._unique_id += 1
        return FullyRandomSequences._unique_id - 1
    
    def __init__(self):
        FigureObject.__init__(self)
        
        self.data_short = load_numbered('length_short', seq_count=100)
        self.data_long = load_numbered('length_longs',seq_count=100)
        self.data_full = dict( self.data_long.items() + self.data_short.items() )
        
    @figuresetup( 1, yscale='linear', xlim=(0,205) )
    def plain_full_no_log(self):
        plotdata( self=self, data=self.data_full, label="{0} (Random Sequences)" )
        plotdata( self=self, data=self.data_full, avg=True )

    @figuresetup(2,  yscale='log', xlim=(0,205) )
    def plain_full_log(self):
        plotdata( self=self, data=self.data_full, label="{0} (Random Sequences)" )
        plotdata( self=self, data=self.data_full, avg=True )

    @figuresetup(3,  yscale='log', xlim=(0,205) )
    def shading_full(self):
        plotdata( self=self, data=self.data_full, label="{0} (Random Sequences)", alpha=.8 )
        plotdata( self=self, data=self.data_full, interval=.4, alpha=.3 )
        plotdata( self=self, data=self.data_full, interval=.05 )
        plotdata( self=self, data=self.data_full, interval=.2, alpha=.4)
        


if __name__ == '__main__':
    Figures_FullRandom = FullyRandomSequences()



