# import matplotlib
# matplotlib.use('macosx')

# import matplotlib.pylab as pyp
import os
import sys

if 'MULTISTRANDHOME' in os.environ:
    if not os.path.isfile( os.path.join( os.environ['MULTISTRANDHOME'], 'setup.py') ):
        warnings.warn( ImportWarning("Could not find the file 'setup.py' in your MULTISTRANDHOME [{0}]; this environment variable is possibly out of date or not referring to the new Mulistrand distribution."))
        multihome=None
    else:
        if os.environ['MULTISTRANDHOME'] not in sys.path:
            multihome= os.environ['MULTISTRANDHOME']     

idx = os.getcwd().find( os.path.join('test','speed') )
if idx > -1:
    rootpath = os.path.abspath(os.getcwd()[:idx])
    #abspath cleans up the tail of the pathname as needed.
    if rootpath not in sys.path:
        sys.path.append( rootpath )
elif multihome != None:
    sys.path.append(multihome)



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
        self.data_short = SpeedTestData( 'data_short', 'length_short')
        self.data_long = SpeedTestData( 'data_long', 'length_longs' )
        self.data_combined = SpeedTestData('data_combined', ['length_short','length_longs'] )
        self.data_vlong = SpeedTestData( 'data_vlong', 'length_very_longs')
        self.data_bm4 = SpeedTestData( 'data_bm4', 'short_4way')
        self.data_bm3 = SpeedTestData( 'data_bm3', 'short_3way')
        self.data_bm4_tc = SpeedTestData( 'data_bm4_tc', 'new_fourway')
        self.data_bm4_tc.normalize(base_data = self.data_bm4, base_data_key='Kinfold')
        self.data_combined_tc = SpeedTestData( 'data_combined_tc', 'length_combined')
                
        #self.data_combined = self.data_short.copy()
        # self.data_combined.update( self.data_long )
        # self.data_quitefull = self.data_combined.copy()
        # self.data_quitefull.update( self.data_vlong )

    def dump_all_data(self):
        for v in ['data_short', 'data_vlong', 'data_long', 'data_bm3',
                  'data_bm4', 'data_bm4_tc']:
            self.__getattribute__(v).to_file()

        
    @figuresetup( 1, yscale='linear', xlim=(0,205) )
    def plain_full_no_log(self):
        plotdata( self=self, data=self.data_combined, label="{0} (Random Sequences)" )
        plotdata( self=self, data=self.data_combined, avg=True )

    @figuresetup(2,  yscale='log', xlim=(0,105) )
    def plain_full_log(self):
        plotdata( self=self, data=self.data_combined, label="{0} (Random Sequences)" )
        plotdata( self=self, data=self.data_combined, avg=True )
        
    @figuresetup(3,  yscale='log', xlim=(0,205) )
    def shading_full(self):
        plotdata( self=self, data=self.data_combined, label="{0} (Random Sequences)", alpha=.8 )
        plotdata( self=self, data=self.data_combined, interval=.4, alpha=.3 )
        plotdata( self=self, data=self.data_combined, interval=.05 )
        plotdata( self=self, data=self.data_combined, interval=.2, alpha=.4)

        
    @figuresetup(4,  yscale='log', xlim=(0,205) )
    def user_sys_comparison(self):
        plotdata( self=self, data=self.data_combined, label="{0} (Real)" )
        plotdata( self=self, data=self.data_combined, avg=True )
        plotdata( self=self, data=self.data_combined, label="{0} (User+Sys)", user=True, offset=1.0, colors={'Kinfold':'m','Multistrand':'c'} )
        plotdata( self=self, data=self.data_combined, avg=True, user=True, offset=1.0, colors={'Kinfold':'m','Multistrand':'c'} )

    @figuresetup(5,  yscale='log', xlim=(0,105) )
    def user_sys_comparison_small(self):
        plotdata( self=self, data=self.data_combined, xrange=(0,105), s=.2, label="{0} (Real)" )
        plotdata( self=self, data=self.data_combined, xrange=(0,105), avg=True )
        plotdata( self=self, data=self.data_combined, xrange=(0,105), s=.2, label="{0} (User+Sys)", user=True, offset=1.0, colors={'Kinfold':'m','Multistrand':'c'} )
        plotdata( self=self, data=self.data_combined, xrange=(0,105), avg=True, user=True, offset=1.0, colors={'Kinfold':'m','Multistrand':'c'} )

    @figuresetup(6,  yscale='log', xlim=(0,105) )
    def user_sys_comparison_areas(self):
        plotdata( self=self, data=self.data_combined, xrange=(0,105), interval=.2)
        plotdata( self=self, data=self.data_combined, xrange=(0,105), interval=.2,  user=True, alpha=.3)
        #        plotdata( self=self, data=self.data_combined, xrange=(0,105), s=.2, label="{0} (Real)", alpha=.5 )
        #       plotdata( self=self, data=self.data_combined, xrange=(0,105), s=.2, label="{0} (User+Sys)", user=True, offset=1.0, colors={'Kinfold':'g','Multistrand':'g'}, alpha = 0.5 )
        
    @figuresetup(7,  yscale='log', xlim=(0,305) )
    def longer_full(self):
        plotdata( self=self, data=self.data_quitefull, label="{0} (Random Sequences)", alpha=.7 )
        plotdata( self=self, data=self.data_quitefull, avg=True, alpha=.7 )
        plotdata( self=self, data=self.data_quitefull, avg=True, user=True, colors={'Kinfold':'m','Multistrand':'c'} )

    @figuresetup(9,yscale='log',xlim=(40,300))
    def fourway_data( self ):
        plotdata(self=self, data=self.data_bm4, label="{0} (4-way)", alpha=.8)
        plotdata(self=self, data=self.data_bm4, avg=True,colors={'Kinfold':'k','Multistrand':'k'}, alpha=.5)
        plotdata(self=self, data=self.data_bm4, interval=.2, alpha=.4)
        plotdata(self=self, data=self.data_bm4, interval=.05)
        plotdata(self=self, data=self.data_bm4, interval=.4, alpha=.3)

    @figuresetup(10,yscale='log',xlim=(40,202))
    def fourway_data_tcmalloc( self ):
        plotdata(self=self, data=self.data_bm4_tc, avg=True, colors={'Kinfold':'k','Multistrand':'k'}, alpha=.5)
        plotdata(self=self, data=self.data_bm4_tc, label="{0} (4-way)", alpha=.8)
        plotdata(self=self, data=self.data_bm4_tc, interval=.2, alpha=.4)
        plotdata(self=self, data=self.data_bm4_tc, interval=.05)
        plotdata(self=self, data=self.data_bm4_tc, interval=.4, alpha=.3)

    @figuresetup(11,yscale='log',xlim=(40,202))
    def fourway_combined(self):
        plotdata(self=self, data=self.data_bm4   , avg=True, colors={'Kinfold':'k','Multistrand':'k'}, alpha=.5)
        plotdata(self=self, data=self.data_bm4_tc, avg=True, colors={'Kinfold':'k','Multistrand':'k'}, offset=1.0, alpha=0.5 )

        plotdata(self=self, data=self.data_bm4   , label="{0} (4-way; old malloc)", alpha=.7)
        plotdata(self=self, data=self.data_bm4_tc, label="{0} (4-way; tcmalloc)"  , alpha=.7, colors={'Kinfold':'m','Multistrand':'c'}, offset=1.0 )

        plotdata(self=self, data=self.data_bm4, interval=.3, alpha=.3)
        plotdata(self=self, data=self.data_bm4_tc, interval=.3, alpha=.3, colors={'Kinfold':'m','Multistrand':'c'}, offset=1.0 )

    @figuresetup(12,yscale='log',xlim=(40,202))
    def fourway_normalized(self):
        plotdata(self=self, data=self.data_bm4   , avg=True, alpha=0.5)

        plotdata(self=self, data=self.data_bm4   , label="{0} (4-way; old malloc)", alpha=.7)
        plotdata(self=self, data=self.data_bm4_tc, label="{0} (4-way; tcmalloc)"  , alpha=.7, colors={'Kinfold':'m','Multistrand':'c'}, offset=1.0 )

        plotdata(self=self, data=self.data_bm4, interval=.3, alpha=.3)
        plotdata(self=self, data=self.data_bm4_tc, interval=.3, alpha=.3, colors={'Kinfold':'m','Multistrand':'c'}, offset=1.0 )

    @figuresetup(13,yscale='log',xlim=(0,205))
    def combined_tc( self ):
        plotdata( self=self, data=self.data_combined_tc, avg=True, colors={'Kinfold':'k','Multistrand':'k'}, alpha=.5 )
        
        plotdata( self=self, data=self.data_combined_tc, label="{0}", alpha=.8 )
        plotdata( self=self, data=self.data_combined_tc, interval=.4, alpha=.3 )
        plotdata( self=self, data=self.data_combined_tc, interval=.05 )
        plotdata( self=self, data=self.data_combined_tc, interval=.2, alpha=.4)

    @figuresetup(15,yscale='log',xlim=(22,200))
    def threeway_data( self ):
        plotdata(self=self, data=self.data_bm3, label="{0} (3-way)", alpha=.8)
        plotdata(self=self, data=self.data_bm3, avg=True)
        plotdata(self=self, data=self.data_bm3, interval=.2, alpha=.4)
        plotdata(self=self, data=self.data_bm3, interval=.05)
        plotdata(self=self, data=self.data_bm3, interval=.4, alpha=.3)

    # @figuresetup(8, yscale='log', xlim=(18,102) )
    # def alternate_colors_log( self, primary=False ):
    #     if primary:
    #         plotdata( self=self, data=self.data_combined,
    #                   label="{0} (No cheats)")
    #         plotdata( self=self, data=self.data_combined, avg=True )
    #     else:
    #         plotdata( self=self, data=self.data_combined,
    #                   label="{0} (Kinfold cheats)",offset=1.0,
    #                   colors={'Kinfold':'m','Multistrand':'c'} )
    #         plotdata( self=self, data=self.data_combined, avg=True,
    #                   offset=1.0, colors={'Kinfold':'m','Multistrand':'c'} )


if __name__ == '__main__':
    fr = FullyRandomSequences()


