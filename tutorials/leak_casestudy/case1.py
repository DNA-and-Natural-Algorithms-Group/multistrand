# Frits Dannenberg, Caltech, 2016. 
# fdann@caltech.edu 

# Example of leak rates, AT ends vs CG ends.
# Compatible with Multistrand 2.1 or higher
# Visit www.multistrand.org

from msUtil import *
from msPlot import *

SCRIPT_DIR = "Hybridization_F4"
ATIME_OUT = 10.0

myMultistrand.setNumOfThreads(14) 

def first_step_simulation(strand_seq, trials, T=25, material="DNA"):

    print ("Running first step mode simulations for %s (with Boltzmann sampling)..." % (strand_seq))
    
    # Using domain representation makes it easier to write secondary structures.
    onedomain = Domain(name="itall", sequence=strand_seq)
    gdomain = Domain(name="Ge", sequence="GCGC")
    
    top = Strand(name="top", domains=[onedomain])
    bot = top.C
    dangle = Strand(name="Gee", domains=[onedomain, gdomain])
    
    duplex_complex = Complex(strands=[top, bot], structure="(+)")
    invader_complex = Complex(strands=[dangle], structure="..")
    duplex_invaded = Complex(strands=[dangle, bot], structure="(.+)")
    
    # Declare the simulation complete if the strands become a perfect duplex.
    success_stop_condition = StopCondition("SUCCESS", [(duplex_invaded, Exact_Macrostate, 0)])
    failed_stop_condition = StopCondition(STR_FAILURE, [(duplex_complex, Exact_Macrostate, 0)])
    
    for x in [duplex_complex, invader_complex]:
        
        x.boltzmann_count = trials
        x.boltzmann_sample = True
         
    
    def getOptions(material, trials, duplex_complex, dangle, success_stop_condition, failed_stop_condition):
        o = Options(simulation_mode="First Step", parameter_type="Nupack", substrate_type=material,
                    rate_method="Metropolis", num_simulations=trials, simulation_time=ATIME_OUT,
                    dangles="Some", temperature=T, rate_scaling="Calibrated", verbosity=0)
        o.start_state = [duplex_complex, dangle]
        o.stop_conditions = [success_stop_condition, failed_stop_condition]
        return o
    
    myMultistrand.setOptionsFactory6(getOptions, material, trials, duplex_complex, invader_complex, success_stop_condition, failed_stop_condition)
    myMultistrand.run()
    dataset = myMultistrand.results

    # Now determine the reaction model parameters from the simulation results.  (Simplified from hybridization_first_step_mode.py.)
    was_success = np.array([1 if i.tag == "SUCCESS" else 0 for i in dataset])
    was_failure = np.array([0 if i.tag == "SUCCESS" else 1 for i in dataset])

    print("Was success:", np.sum(was_success))
    print("Was failure:", np.sum(was_failure))



def doFirstStepMode(seq, T=25, material="DNA", numOfRuns=500):
    
    first_step_simulation(seq, numOfRuns, T=T, material=material) 
   

def makePlots():


    seqs = list()

    seqs.append('GTCGATGC')
    seqs.append('TCGAGTGA')

   
    for seq in seqs:
        doFirstStepMode(seq, numOfRuns=500000)
     


# The actual main method
if __name__ == '__main__':
        
        makePlots()















