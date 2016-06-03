import os
import subprocess
import time


#Sample usage:
#
# c = ClusterInfo()
# c.cluster_status( 'schaeffer' )
#
# c.machine_info('zig-zag')
# c.machine_info('nucleotide1')
# c.machine_info()
#
# It will not try to poll the servers for data until you do a
# cluster_status command (or the underlying helper function).
#  A lot of data is kept around, try inspecting various internals.

username = 'someone'

machines = ['nucleotide1','nucleotide2','nucleotide3','nucleotide4','nucleotide5','nucleotide6','nucleotide7','nucleotide8','nucleotide9','nucleotide10','nucleotide11','nucleotide12','nucleotide13','nucleotide14','nucleotide15','nucleotide16','a-form','b-form','z-form','zig-zag','twisted','melted']

class ClusterInfo(object):
    def __init__(self):
        self.res = []
        self.procdata = {}
        self.lastupdate = ('None',-500)


    def check_machines(self, username ):
        self.lastupdate = (time.asctime(),time.time())
        self.res = []
        for line in machines:
            try:
                address = line + '.dna.caltech.edu'
                p1 = subprocess.Popen(["ssh", username + '@' + address],stdout = subprocess.PIPE, stderr = subprocess.PIPE, stdin=subprocess.PIPE)

                p1.stdin.write("cat /proc/loadavg\n")

                loads_line = p1.stdout.readline()
                loadavg = loads_line.split()[0:3]

                p1.stdin.write('cat /proc/cpuinfo | grep processor | wc -l\n')
                processor = p1.stdout.readline().strip()


                p1.stdin.write("cat /proc/cpuinfo | grep 'core id' | wc -l\n")
                cores = p1.stdout.readline().strip()

                if cores == '': cores = processor
                if processor == '': processor = cores
                if processor == '': processor = '-1000'
                

                p1.stdin.write('top -n 1 -b\n')

                (data_stdout,data_stderr) = p1.communicate()

                data_lines = data_stdout.split('\n')
                for l in data_lines:
                    if len(l) == 0:
                        rest_data = data_lines[data_lines.index(l)+1:]
                        break

                if line in self.procdata:
                    del self.procdata[line]
                self.procdata[line] = []
                # data starts from rest_data[1], shoudl be sorted by cpu.
                for l in rest_data[1:]:
                    if len(l) == 0:
                        break
                    pid,user,pr,ni,virt,resident,shr,stat,cpu,mem,cputime,cmd = l.split(None,11)
                    self.procdata[line].append( (user, stat, cpu, mem, cmd ) )

                self.res.append((address, loadavg, processor, cores, line))
                print("{0}: Cores: {cpucount} Load: {load[0]}/{load[1]}/{load[2]}".format( line, load=loadavg, cpucount=max(cores,processor)))
                time.sleep(.2)
            except IOError:
                print("{0} is not responding, or there were other issues.".format( address ))

    def machine_info(self, name = None ):
        """ display info about the machine from our last request, or all machines if no param. """
        if name == None:
            for nm in machines:
                if nm != None:
                    self.machine_info( nm )
            return
        if not name in self.procdata:
            return
        
        def machine_load( machine_name ):
            hits = [i for i in self.res if i[4] == machine_name ]
            if len(hits) == 0:
                raise IndexError("No machine with that name found in our results list.")
            else:
                return hits[-1]  # the most recent result for that processor.

        loadinfo = machine_load( name )
        print("Machine status: {0}\nCPU Processors/Cores: {1[2]}/{1[3]} Load: {1[1][0]}/{1[1][1]}/{1[1][2]} (Average # of CPUs used in last 1m/5m/15m)".format( name, loadinfo ))
        print("{0[0]:<12}| {0[1]:4} | {0[2]:>5} {0[3]:>5} | {0[4]:<40}".format(['user','stat','cpu','mem','cmd']))
        procs = self.procdata[ name ]
        procs.sort(key = lambda x: float(x[2]), reverse=True)
        for p in procs[:10]:
            print("{0[0]:<12}| {0[1]:4} | {0[2]:>5} {0[3]:>5} | {0[4]:<40}".format(p))
        if len(procs) > 10:
            print("... [[{0} more processes omitted]] ... ".format( len(procs) - 10 ) )
        print("")
        
    def cluster_status( self, uname = None, timeoverride=False ):
        if not timeoverride and ( time.time() - self.lastupdate[1]) <= 60:
            print("Warning: you recently ran this command, possibly you should wait before running it again?")
            return
        
        def reducer_cores(x,y):
            return [x[0] + float(y[1][0]),
                    x[1] + float(y[1][1]),
                    x[2] + float(y[1][2]),
                    x[3] + int(y[2]),
                    x[4] + int(y[3])]

        def corecount(d):
            return reduce( reducer_cores, d, [0.0,0.0,0.0,0,0])

        try:
            self.check_machines( uname or username)
        except NameError:
            self.check_machines('dummy')
        corevals = corecount( self.res )
        print("DNA Cluster Usage, {ctime}\nLoads: [{load[0]}/{cpucount}] - {usepercent:.2}%  \
        5/15 Load: {load[1]}/{load[2]}".format(load=corevals[0:3],cpucount=corevals[4], usepercent =100.0 * corevals[0]/corevals[4], ctime=time.asctime()))

        try:
            self.data_timeseries[time.asctime()] = corevals
        except:
            self.data_timeseries = {}
            self.data_timeseries[time.asctime()] = corevals

    
    
