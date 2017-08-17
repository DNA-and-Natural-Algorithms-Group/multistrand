
from multistrand.objects import Strand, Complex, Domain, StopCondition
from multistrand.options import Options
from multistrand.system import SimSystem
from multistrand.experiment import makeComplex

def createOptions(start_complex, stop_complex, simMode):

    full_sc = StopCondition("CLOSED", [(stop_complex, Options.dissocMacrostate, 2)])    
    
    o1 = Options(simulation_mode=simMode,  # "First Passage Time", 
                 parameter_type="Nupack", substrate_type="DNA", temperature=273.15 + 25.0,
                num_simulations=10,
                simulation_time=0.00001,
                rate_scaling='Calibrated',
                verbosity=0,
                join_concentration=1e-9,
                rate_method="Metropolis",
                start_state=[start_complex],
                stop_conditions=[full_sc])
   

    
    return o1

        

def create_test0():


    toehold_seq = "CCCC"
    domain_seq = "CATTAAC"

    # build complexes with domain-level information    
    toehold = Domain(name="toehold", sequence=toehold_seq, length=4)
    branch_migration = Domain(name="branch_migration", sequence=domain_seq, seq_length=7)        
        
    incoming = branch_migration.C + toehold.C  
    substrate = toehold + branch_migration
    
    start_complex = Complex(strands=[incoming, substrate], structure="((+))")
    stop_complex = Complex(strands=[incoming, substrate], structure="..+..") 

    
    return createOptions(start_complex, stop_complex, "First Passage Time")





def create_test0B():


    toehold_seq = "CC"
    domain_seq = "CAAC"

    # build complexes with domain-level information    
    toehold = Domain(name="toehold", sequence=toehold_seq, length=2)
    branch_migration = Domain(name="branch_migration", sequence=domain_seq, seq_length=4)
        
    incoming = branch_migration.C + toehold.C  
    substrate = toehold + branch_migration
        
#     start_complex = Complex(strands=[incoming, substrate], structure="..+..")
#     stop_complex = Complex(strands=[incoming, substrate], structure="((+))") 


    start_complex = Complex(strands=[incoming, substrate], structure="((+))")
    stop_complex = Complex(strands=[incoming, substrate], structure="..+..") 

    
    
    return createOptions(start_complex, stop_complex, "First Passage Time")



def create_test0C():


    domain_seq = "AGT"
    domain_seq2 = "GTA"


    left = Domain(name="branch_migration", sequence=domain_seq, seq_length=3)
    right =     Domain(name="branch_migration2", sequence=domain_seq2, seq_length=3)
        
    incoming = left + right
    substrate =  incoming.C
        
    start_complex1 = Complex(strands=[incoming], structure="..") 
    start_complex2 = Complex(strands=[substrate], structure="..")
    stop_complex = Complex(strands=[incoming, substrate], structure="((+))")

    full_sc = StopCondition("CLOSED", [(stop_complex, Options.dissocMacrostate, 2)])  
    
    
    o1 = Options(simulation_mode=Options.firstPassageTime,  # "First Passage Time", 
                temperature=273.15 + 25.0,
                num_simulations=10,
                simulation_time=0.00001,
                rate_scaling='Calibrated',
                verbosity=0,
                join_concentration=1e-9,
                rate_method="Metropolis",
                start_state=[start_complex1, start_complex2],
                stop_conditions=[full_sc])
   
    return o1




def create_test1():
 
 
    top0 = "ACT"
    top1 = "GAC"
    toehold = "TG"
    bottom = "TG"
 
    # build complexes with domain-level information    
    right_d = Domain(name="toehold0", sequence=top0, length=3)
    left_d = Domain(name="toehold1", sequence=top1, length=3)
    branch0 = Domain(name="branch_migration", sequence=bottom, seq_length=2)
    toehold = Domain(name="oToehold", sequence=toehold, length=2)
            
    substrate = toehold.C + branch0.C + toehold.C
    left = toehold + left_d
    right = right_d + toehold  
     
 
    # Note that "+" is used to indicate strand breaks.  
    # So the initial structures represent the incoming strand bound by its toehold,
    # and we'll see that either it completes strand displacement, or it dissociates.
    start_complex = Complex(strands=[left, right, substrate], structure="(.+.(+).)")
    stop_complex = Complex(strands=[left, right, substrate], structure="..+..+...") 
     
    return createOptions(start_complex, stop_complex, "First Passage Time")
    
    
def create_test1B():


    top0 = "ACT"
    top1 = "GAC"
    toehold = "TG"
    bottom = ""

    # build complexes with domain-level information    
    right_d = Domain(name="toehold0", sequence=top0, length=3)
    left_d = Domain(name="toehold1", sequence=top1, length=3)
    branch0 = Domain(name="branch_migration", sequence=bottom, seq_length=2)
    toehold = Domain(name="oToehold", sequence=toehold, length=2)
           
    substrate = toehold.C + branch0.C + toehold.C
    left = toehold + left_d
    right = right_d + toehold  
    

    # Note that "+" is used to indicate strand breaks.  
    # So the initial structures represent the incoming strand bound by its toehold,
    # and we'll see that either it completes strand displacement, or it dissociates.
    start_complex = Complex(strands=[left, right, substrate], structure="(.+.(+).)")
    stop_complex = Complex(strands=[left, right, substrate], structure="..+..+...") 
    
    return createOptions(start_complex, stop_complex, "First Passage Time")
    





def create_test2():


    top0 = "T"
    top1 = "T"
    top2 = "G"

    bottom0 = "C"
    bottom1 = "A"
    bottom2 = "T"


    # build complexes with domain-level information    
    strand0 = Domain(name="toehold0", sequence=top0, length=1)
    strand1 = Domain(name="toehold1", sequence=top1, length=1)
    strand2 = Domain(name="toehold2", sequence=top2, length=1)
    
    strand3 = Domain(name="toehold3", sequence=bottom0, length=1)
    strand4 = Domain(name="toehold4", sequence=bottom1, length=1)
    strand5 = Domain(name="toehold5", sequence=bottom2, length=1)
   
    substrate = strand3 + strand4 + strand5
    invading = strand0 + strand1 + strand2
    

    start_complex = Complex(strands=[substrate, invading], structure=".(.+.).")
    stop_complex = Complex(strands=[substrate, invading], structure="...+...") 
    
    return createOptions(start_complex, stop_complex, "First Passage Time")
    
    



def create_test3():

    
    strand_seq = "CTGA"
    num_traj = 10

   # Essentially, creates the options object and prepares to simulate the hybridization of the strand and its complement.
    onedomain = Domain(name="itall", sequence=strand_seq)
    top = Strand(name="top", domains=[onedomain])
    bot = top.C

    # Note that the structure is specified to be single stranded, but this will be over-ridden when Boltzmann sampling is turned on.
    start_complex_top = Complex(strands=[top], structure=".")
    start_complex_bot = Complex(strands=[bot], structure=".")
    start_complex_top.boltzmann_count = num_traj
    start_complex_bot.boltzmann_count = num_traj
    start_complex_top.boltzmann_sample = True
    start_complex_bot.boltzmann_sample = True
    # Turns Boltzmann sampling on for this complex and also does sampling more efficiently by sampling 'num_traj' states.

    # Stop when the exact full duplex is achieved. (No breathing!)
    success_complex = Complex(strands=[top, bot], structure="(+)")
    success_stop_condition = StopCondition("SUCCESS", [(success_complex, msUtil.Exact_Macrostate, 0)])

    # Declare the simulation unproductive if the strands become single-stranded again.
    failed_complex = Complex(strands=[top], structure=".")
    failed_stop_condition = StopCondition("FAILURE", [(failed_complex, msUtil.Dissoc_Macrostate, 0)])

    o = Options(simulation_mode="First Step",
                parameter_type="Nupack",
                substrate_type="DNA",
                rate_method="Metropolis",
                num_simulations=num_traj,
                simulation_time=1.0,
                dangles="Some",
                temperature=273.15 + 25.0,
                rate_scaling="Calibrated",
                useArrRates=True,
                verbosity=0)

    o.start_state = [start_complex_top, start_complex_bot]
    o.stop_conditions = [success_stop_condition, failed_stop_condition]
    return o







colors = ['blue', 'red', 'cyan', 'magenta', 'green', 'k', 'darkblue', 'darkred', 'darkcyan', 'darkmagenta', 'darkgreen']

def create_test4():

    toehold_t = "CTGC"
    toehold_dd = "CATATC"
    domain_R = "CATTAAC"

    # build complexes with domain-level information    
    toehold = Domain(name="toehold", sequence=toehold_t, length=6)
    toehold_2 = Domain(name="toehold", sequence=toehold_dd, length=6)
    branch_migration = Domain(name="branch_migration", sequence=domain_R, seq_length=7)
        
    incoming = branch_migration.C + toehold.C
    substrate = toehold + branch_migration

    # Note that "+" is used to indicate strand breaks.  
    # So the initial structures represent the incoming strand bound by its toehold,
    # and we'll see that either it completes strand displacement, or it dissociates.
    start_complex = Complex(strands=[incoming, substrate], structure=".(+).")
    stop_complex = Complex(strands=[incoming, substrate], structure="..+..") 
    
    full_sc = StopCondition("CLOSED", [(stop_complex, msUtil.Dissoc_Macrostate, 2)]) 
    
    return createOptions(start_complex, stop_complex, "First Passage Time")
    





def create_test5():


    
    top0 = "ACT"
    top1 = "GAC"
    toehold = "TG"
    bottom = "TG"
    connected = "ATA";
    

    # build complexes with domain-level information    
    right_d = Domain(name="toehold0", sequence=top0, length=3)
    left_d = Domain(name="toehold1", sequence=top1, length=3)
    branch0 = Domain(name="branch_migration", sequence=bottom, seq_length=2)
    toehold = Domain(name="oToehold", sequence=toehold, length=2)
    connect = Domain(name="cToehold", sequence=connected, length=3)
           
    substrate = toehold.C + branch0.C + toehold.C
    left = toehold + left_d + connect
    right = connect.C + right_d + toehold 
    

    # Note that "+" is used to indicate strand breaks.  
    # So the initial structures represent the incoming strand bound by its toehold,
    # and we'll see that either it completes strand displacement, or it dissociates.
    start_complex = Complex(strands=[left, right, substrate], structure="(.(+).(+).)")
    stop_complex = Complex(strands=[left, right, substrate], structure="...+...+...") 
    
    return createOptions(start_complex, stop_complex, "First Passage Time")    


def create_test6():

    toehold_seq = "CTGC"
    toehold_seq2 = "CAT"
    domain_seq = "CATGCTAAC"

    # build complexes with domain-level information    
    toehold = Domain(name="toehold", sequence=toehold_seq, length=4)
    toehold2 = Domain(name="toehold", sequence=toehold_seq2, length=3)
    branch_migration = Domain(name="branch_migration", sequence=domain_seq, seq_length=9)

        
    incoming = toehold2.C + branch_migration.C + toehold.C  
    substrate = toehold + branch_migration + toehold2
        
    start_complex = Complex(strands=[incoming, substrate], structure="(.(+).)")
    stop_complex = Complex(strands=[incoming, substrate], structure="...+...") 
    
    return createOptions(start_complex, stop_complex, "First Passage Time")


def create_test6B():

    toehold_seq = "CCC"
    toehold_seq2 = "TTT"
    domain_seq = "AA"

    # build complexes with domain-level information    
    toehold = Domain(name="toehold", sequence=toehold_seq, length=3)
    toehold2 = Domain(name="toehold", sequence=toehold_seq2, length=3)
    branch_migration = Domain(name="branch_migration", sequence=domain_seq, seq_length=2)

        
    incoming = toehold2.C + branch_migration.C + toehold.C  
    substrate = toehold + branch_migration + toehold2
        
    start_complex = Complex(strands=[incoming, substrate], structure="(.(+).)")
    stop_complex = Complex(strands=[incoming, substrate], structure="...+...") 
    
    return createOptions(start_complex, stop_complex, "First Passage Time")




def create_test7():

    toehold_seq = "CTGC"
    toehold_seq2 = "CAT"
    domain_seq = "CATGCTAAC"

    # build complexes with domain-level information    
    toehold = Domain(name="toehold", sequence=toehold_seq, length=6)
    toehold2 = Domain(name="toehold", sequence=toehold_seq2, length=6)
    branch_migration = Domain(name="branch_migration", sequence=domain_seq, seq_length=25)
    dangle = Domain(name="branch_migration", sequence="T", seq_length=1) 
        
    incoming = toehold2.C + branch_migration.C + toehold.C  
    substrate = toehold + toehold2
        
    start_complex = Complex(strands=[incoming, substrate], structure="(.(+))")
    stop_complex = Complex(strands=[incoming, substrate], structure="...+..") 
    
    return createOptions(start_complex, stop_complex, "First Passage Time")

def create_test8():

    toehold_seq = "CTGC"
    domain_seq = "CATGCTACAG"

    # build complexes with domain-level information    
    toehold = Domain(name="toehold", sequence=toehold_seq, length=6)
    branch_migration = Domain(name="branch_migration", sequence=domain_seq, seq_length=25)
    dangle = Domain(name="branch_migration", sequence="T", seq_length=1) 
        
    incoming = toehold + branch_migration.C + toehold.C  

        
    start_complex = Complex(strands=[incoming], structure="(.)")
    stop_complex = Complex(strands=[incoming], structure="...") 
    
    return createOptions(start_complex, stop_complex, "First Passage Time")


def create_test9():


    seq0 = "GTGT"
    seq1 = "T"
    
    

    # build complexes with domain-level information    
    branch = Domain(name="toehold0", sequence=seq0, length=3)
    toehold = Domain(name="toehold1", sequence=seq1, length=1)
           
    ghost = Domain(name="toeholdG", sequence="T", length=1)
   
    
    substrate = toehold + branch
    left = toehold.C + ghost
    right = branch.C + ghost
    

    # Note that "+" is used to indicate strand breaks.  
    # So the initial structures represent the incoming strand bound by its toehold,
    # and we'll see that either it completes strand displacement, or it dissociates.
    start_complex = Complex(strands=[right, left, substrate], structure="(.+(.+))")
    stop_complex = Complex(strands=[left, right, substrate], structure="..+..+..") 
    
    return createOptions(start_complex, stop_complex, "First Passage Time")


def create_test10():


    seq0 = "GTCACTGCTTTT"
    seq1 = "GCAGTGAC"
    dotparen1 = "..((((((....+)))))).."
    
    
    seq2 = "GTCACTGC"
    dotparen2 = "........"
    
    complex1 = makeComplex([seq0, seq1], dotparen1)
    complex2 = makeComplex([seq2], dotparen2 )
        
    return createOptions(complex1, complex2, "First Passage Time")




def main():
    
        
#     o1 = create_test0()      # just a fully hybridized strand.
#     o1 = create_test0B()      # just a fully hybridized strand.
#     o1 = create_test0C()      # just a fully hybridized strand.
#     o1 = create_test1()      # testing open-loop 
#     o1 = create_test1B()      # testing open-loop with the initialiation penalty 
#     o1 = create_test2()      # a very simple test
#     o1 = create_test3()  # this is the  bi-molecular test
#     o1 = create_test4()  # this is the lightbulb
#     o1 = create_test5()      # testing multi-loop code
#    o1 = create_test6()      # interior loop.
#     o1 = create_test6B()      # small interior loop.
#     o1 = create_test7()      # Bulge loop.
#     o1 = create_test8()      # Hairpin loop.
#     o1 = create_test9()      # Small open-loop code
    o1 = create_test10()      # half open duplex



#     setArrheniusConstantsNM1(o1)
    
    s = SimSystem(o1)
    # s.start()
    s.initialInfo()
    # s.calculate_rate
#     print "Testing loop internals"

    
    # s.start()
    # s.InitializeSystem()


main()








