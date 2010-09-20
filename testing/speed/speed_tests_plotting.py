
import matplotlib
matplotlib.use('macosx')

import matplotlib.pylab as pyp


f = open('results_python.txt')

data = eval( f.readline() )
f.close()

n100_t1000_s20 = [i for i in data if i[0] is not None and len(i[0])==20]
times_s20 = [i[4] for i in n100_t1000_s20]

times_s20_kin = [i[0] for i in times_s20]
times_s20_ms = [i[2] for i in times_s20]

def show_20_plot(show = False):
    pyp.figure(1)
    kin = pyp.plot( times_s20_kin[:-1] )
    ms = pyp.plot( times_s20_ms[:-1] )
    pyp.xlabel("Sequence Composition (indexed)")
    pyp.ylabel("Real Clock Time (s)")
    if show:
        pyp.show()


n100_t1000_s40 = [i for i in data if i[0] is not None and len(i[0])==40]
times_s40 = [i[4] for i in n100_t1000_s40]

times_s40_kin = [i[0] for i in times_s40]
times_s40_ms = [i[2] for i in times_s40]

def show_40_plot(show=False):
    pyp.figure(2)
    kin = pyp.plot( times_s40_kin[:-1] )
    ms = pyp.plot( times_s40_ms[:-1] )
    pyp.xlabel("Sequence length (indexed)")
    pyp.ylabel("Real Clock Time (s)")
    if show:
        pyp.show()



length_run = [i for i in data if i[2] == 10]
times_length = [i[4] for i in length_run]
times_length_kin = [i[0] for i in times_length]
times_length_ms = [i[2] for i in times_length]
lengths = [20,40,60,80,100]

def show_length_plot(show=False):
    pyp.figure(3)
    kin = pyp.plot( lengths, times_length_kin )
    ms = pyp.plot( lengths, times_length_ms )
    pyp.xlabel("Sequence length (nt)")
    pyp.ylabel("Real Clock Time (s)")
    if show: pyp.show()

Figures = []

def figuresetup( num, xlabel=None, ylabel=None, title = None ):
    def actual_decorator( funct ):
        def res( show=False ):
            f = pyp.figure( num )
            f.clear()

            if xlabel == None:
                pyp.xlabel("Sequence length (nt)")
            else:
                pyp.xlabel( xlabel )

            if ylabel == None:
                pyp.ylabel("Real Clock Time (RCT) per Simulated Time Unit (STU)")
            else:
                pyp.ylabel( ylabel )

            pyp.ticklabel_format( style ='sci', scilimits = (0,0), axis='y' )
            if title == None:
                pyp.title("Kinfold vs Multistrand comparison on random sequences", position=(0.5,1.05) )
            else:
                pyp.title( title, position=(0.5,1.05) )
                # position it up slightly to make room for scientific notation labels. 

            funct()
            pyp.legend( loc=0)
            if show:
                pyp.show()
        Figures.append( (res, num ) )
        return res
    return actual_decorator


    
@figuresetup( 6 )
def decor_random_length_plot():
    length_run = [i for i in data if i[2] == 10]
    times_scaling = [i[1] * i[2] * i[3] for i in length_run]
       # i[1] is trajectory count
       # i[2] is number of sims
       # i[3] is simulated time units
    times_length = [i[4] for i in length_run]
    times_length_kin_scaled = [i[0]/j for i,j in zip(times_length,times_scaling)]
    times_length_ms_scaled = [i[2]/j for i,j in zip(times_length, times_scaling)]

    lengths = [20,40,60,80,100]
    
    pyp.plot( lengths, times_length_kin_scaled, label="Kinfold")
    pyp.plot( lengths, times_length_ms_scaled, label="Multistrand")

    times_s40 = [i[4] for i in n100_t1000_s40]
    times_s40_kin_scaled = [i[0] / (100 * 1000.0 ) for i in times_s40[:-1]]
    times_s40_ms_scaled = [i[2] / (100 * 1000.0 ) for i in times_s40[:-1]]

    times_s20 = [i[4] for i in n100_t1000_s20]
    times_s20_kin_scaled = [i[0] / (100 * 1000.0 ) for i in times_s20[:-1]]
    times_s20_ms_scaled = [i[2] / (100 * 1000.0 ) for i in times_s20[:-1]]

    pyp.scatter( [20]*(len(times_s20_kin)-1), times_s20_kin_scaled, label="Kinfold [20]", color='b' )
    pyp.scatter( [20]*(len(times_s20_ms)-1), times_s20_ms_scaled, label="Multistrand [20]",marker='+', color='g' ) 
    pyp.scatter( [40]*(len(times_s40_kin)-1) , times_s40_kin_scaled, label="Kinfold [40]", color='b')
    pyp.scatter( [40]*(len(times_s40_ms)-1) , times_s40_ms_scaled, label="Multistrand [40]",marker='+', color = 'g' )

    k_times = [(.1 / 1000.0) for i in range(10,192,2)]
    k_times_2 = [(.2 / 1000.0) for i in range(10,192,2)]
    pyp.plot( range(10,192,2), k_times, label="PlaceholderK")
    pyp.plot( range(10,192,2),k_times_2, label="PlaceholderM")

@figuresetup(4)
def show_20_40():
    pyp.scatter( [20]*(len(times_s20_kin)-1), times_s20_kin[:-1], label = "Kinfold [20]" )
    pyp.scatter( [20]*(len(times_s20_kin)-1), times_s20_ms[:-1], label = "Multistrand [20]" ) 
    pyp.scatter( [40]*(len(times_s40_kin)-1) , times_s40_kin[:-1], label = "Kinfold [40]")
    pyp.scatter( [40]*(len(times_s40_kin)-1) , times_s40_ms[:-1], label = "Multistrand [40]" )

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
            

@figuresetup(6, title = "Kinfold vs Multistrand Comparison, Weak MFE sequences")
@prepdata( None )
def lowmfe_from_dataset( data = None):
    pyp.scatter( data[0][0], data[0][1], label = "Weak MFE [Kinfold]", color='r' )
    pyp.scatter( data[1][0], data[1][1], label = "Weak MFE [Multistrand]", color='m' )

    
@figuresetup(5, title = "Kinfold vs Multistrand Comparison, Strong MFE sequences")
def highmfe_from_dataset():
    highmfe_x_range = (range(10,192,2) * 5)
    highmfe_x_range.sort()
    highmfe_data = [random.random() * (i / 1000.0) * random.random() / 1000.0 for i in range(len( highmfe_x_range))]
    highmfe_xy = ( highmfe_x_range, highmfe_data )
    pyp.scatter( highmfe_xy[0], highmfe_xy[1], label = "Strong MFE [Kinfold]", color='b' )

    highmfe_x_range = (range(10,192,2) * 5)
    highmfe_x_range.sort()
    highmfe_data = [random.random() * (i / 1000.0 ) * random.random() / 1000.0 for i in range(len( highmfe_x_range))]
    highmfe_xy = ( highmfe_x_range, highmfe_data )
    pyp.scatter( highmfe_xy[0], highmfe_xy[1], label = "Strong MFE [Multistrand]", color='c' )



@figuresetup(7, title = "Kinfold vs Multistrand Comparison, Strong vs Weak MFE sequences")
@prepdata( None )
def bothmfe_from_dataset( data = None):
    pyp.scatter( data[0][0], data[0][1], label = "Weak MFE [Kinfold]", color='r' )
    pyp.scatter( data[1][0], data[1][1], label = "Weak MFE [Multistrand]", color='m' )

    highmfe_x_range = (range(10,192,2) * 5)
    highmfe_x_range.sort()
    highmfe_data = [random.random() * (i / 20.0) * random.random() / 1000.0 for i in range(len( highmfe_x_range))]
    highmfe_xy = ( highmfe_x_range, highmfe_data )
    pyp.scatter( highmfe_xy[0], highmfe_xy[1], label = "Strong MFE [Kinfold]", color='b' )

    highmfe_x_range = (range(10,192,2) * 5)
    highmfe_x_range.sort()
    highmfe_data = [random.random() * (i / 20.0 ) * random.random() / 1000.0 for i in range(len( highmfe_x_range))]
    highmfe_xy = ( highmfe_x_range, highmfe_data )
    pyp.scatter( highmfe_xy[0], highmfe_xy[1], label = "Strong MFE [Multistrand]", color='c' )

@figuresetup(8, title = "Kinfold vs Multistrand Comparison, 4-way Branch Migration")
@prepdata( range(70,120,5)*5, lambda x: x * random.random() / 1000 )
def fourway_bm( data = None ):
    pyp.scatter( data[0][0], data[0][1], label = "4-way BM [Kinfold]", color = 'r' )
    pyp.scatter( data[1][0], data[1][1], label = "4-way BM [Multistrand]", color = 'g' )


@figuresetup(9, title = "Kinfold vs Multistrand Comparison, 3-way Branch Migration")
@prepdata( range(70,120,5)*5, lambda x: x * random.random() / 1000 )
def threeway_bm( data = None ):
    pyp.scatter( data[0][0], data[0][1], label = "3-way BM [Kinfold]", color = 'r' )
    pyp.scatter( data[1][0], data[1][1], label = "3-way BM [Multistrand]", color = 'g' )


@figuresetup(10, title = "Kinfold vs Multistrand Comparison, Hairpin Folding")
@prepdata( range(10,80,5)*5, lambda x: x * random.random() / 600 )
def hairpin_data( data = None ):
    pyp.scatter( data[0][0], data[0][1], label = "Hairpin [Kinfold]", color = 'r' )
    pyp.scatter( data[1][0], data[1][1], label = "Hairpin [Multistrand]", color = 'g' )


