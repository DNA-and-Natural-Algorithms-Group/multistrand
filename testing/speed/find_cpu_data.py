import wx,os
import subprocess

machinedata = {}
  # dictionary to contain all the data for each cluster node.

f = open('/Users/zifnab/compute_machines.txt','r')
for line in f:
    if not line.startswith('#'):
        address = 'schaeffer@' + line.rstrip(' \n')+ '.dna.caltech.edu'
        p1 = subprocess.Popen(["ssh", address, "cat /proc/loadavg | awk '{print $1}'; cat /proc/cpuinfo | grep processor | wc -l; cat /proc/cpuinfo | grep 'core id' | wc -l"], stdout = subprocess.PIPE, stderr = subprocess.PIPE )

        (data_stdout,data_stderr) = p1.communicate()
        data = data_stdout.split()

        if len(data) == 3:
            loadavg = data[0]
            processor = data[1]
            cores = data[2]
            print address, loadavg, processor, cores
        else:
            print address, "not reachable."
