# Mrinank Sharma, ms2314@cam.ac.uk August 2017
#
# Here we investigate the use of using 'supersampling', a technique
# in which we reuse each boltzmann sample a number of times.
# Here we will borrow most methods from case2.py

# python default packages
import time

# matlab style plotting packages
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.ticker import ScalarFormatter

# multistrand
from multistrand.concurrent import myMultistrand
from case2 import genOptions, Experiment, setupSimulationOptions, runExperiment, CL_LONG_GATE_A_SEQ, CL_LONG_GATE_B_SEQ
from multistrand.experiment import ClampedSeesawGate


# For a greater understanding of how this specific code actually runs, please see case2.py
def runSimulations():
    # Here we say we are going to use 2 threads, storing only succesful data.
    # We demand at 100 successful trials
    setupSimulationOptions(8, True, 40, 2.5e6, 1000)

    # Here we create two clamped seesaw gates, according to the defined interface.
    gateA = ClampedSeesawGate(*CL_LONG_GATE_A_SEQ)
    gateB = ClampedSeesawGate(*CL_LONG_GATE_B_SEQ)

    trialsPerThreadIncremenet = 100
    trialsIncrement = trialsPerThreadIncremenet * myMultistrand.numOfThreads

    xValues = []
    times = []

    for ssample in range(1, 76, 12):
        start = time.time()
        runExperiment(trialsIncrement, gateA,
                      Experiment.GATE_OUTPUT_PRODUCTION, None, ssample)
        end = time.time()
        elapsedTime = end-start
        timePerTrial = elapsedTime / myMultistrand.results.nTotal
        xValues.append(ssample)
        times.append(timePerTrial)
    
    return xValues, times

def plotFigure():
    plt.rcdefaults()
    results = runSimulations()
    xValues = results[0]
    times = results[1] 
    plt.plot(xValues, times, linewidth=2.0, marker='.', markersize=10.0)
    plt.xlim(xmin=1)
    plt.title("Time vs Supersampling Number")
    plt.xlabel('Supersample')
    plt.ylabel('Time per Trial(s)')
    plt.gca().set_ylim(bottom=0)
    plt.show()

if __name__ == "__main__":
    plotFigure()
