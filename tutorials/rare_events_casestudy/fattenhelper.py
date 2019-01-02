import sys
import cPickle as pickle
import gc
import time
from multistrand.builder import Builder
from test_string_method import associationNoInit
print  "importing associationNoInit  from test_string_method.py in fattenhelper.py "


start = sys.argv[1]
end = sys.argv[2]
pathbuilder = sys.argv[3]
pathspace= sys.argv[4 ]
pathsequences = sys.argv[5]
pathoptions= sys.argv[6]

print "In fattenhelper.py, which uses  associationNoInit from test_string_method to initialize a builder. Then it calls fattenStateSpace_batch, etc  to fatten states  "
#def do(start, end, pathbuilder, pathspace , pathsequences, pathoptions) :
mytime = open ("times.txt", "a")
mytime.write( pathbuilder + "   start " + str( start) + "end "  + str(end) + "\n" )
st = time.time()

with open(pathoptions  , "rb" ) as p:
	optionsArg  = pickle.load( p)

myBuilder=   Builder(associationNoInit,   optionsArg )

with open(pathspace  , "rb" ) as p:
	myBuilder.protoSpacebackup  = pickle.load( p)

with open( pathsequences  , "rb" ) as p :
	myBuilder.protoSequences  = pickle.load( p)

mytime.write( "load time " + str( time.time()  - st )+"\n")

st  = time.time ( )
myBuilder.fattenStateSpace_batch(start = int(start)  , end= int( end))

mytime.write( "fatten time time " + str( time.time()   - st ) +"\n")

st = time.time()

with open( pathbuilder + "pt" + str(start)+"-"  +str(end) , "wb" ) as p:
	pickle.dump(myBuilder.protoTransitions,  p)


mytime.write( "save time " + str( time.time()  - st ) +"\n")
mytime.close()

del myBuilder
gc.collect()
