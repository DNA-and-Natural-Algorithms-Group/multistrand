Multistrand basics
=========================


Multistrand sample input
----------------------------------

Our example sets up threeway branch migration and follows the threeway branch migration demo. The following sequences are used:
:: 
  
    toehold_seq = "TCTCCATGTCACTTC"	
    bm_design = "CCCTCATTCAATACCCTACG" 

To start the simulation, Multistrand should know what the initial state of the simulation is. In this example, we also specify the final state, so that Multistrand knows when to stop
the simulation. The initial state is specified as a list of complexes. Complexes, together with Strands and Domains, are multistrand-provided datastructures. We now creat the domain objects:
::
    toehold = Domain(name="toehold",sequence=toehold_seq[0:toehold_length])
    branch_migration = Domain(name="bm", sequence=bm_design)
    
Note that toehold_seq[0:toehold_length] specifies only a prefix of the toehold. 
After creating the two domains, we now create Strand objects by combining existing domains.  
The modifier .C specifies the reverse complement of the domain.
::
    incoming = branch_migration + toehold
    incoming.name = "incoming"
    
    toehold_extra = Domain(name="toehold_extra",sequence=toehold_seq[toehold_length:]) 
    substrate = toehold_extra.C + toehold.C + branch_migration.C

When a strand consists of just one domain, we have to explicitly call the Domain constructor:
::
    
    incumbent = Strand(name="incumbent",domains=[branch_migration])

Using the incoming strand, the incumbent strand and the substrate strand, we now construct the two complexes that form the initial state:
::
    # The start complex  
    start_complex_substrate_incumbent = Complex(strands=[substrate, incumbent],structure="..(+)")

    # The incoming complex. 
    start_complex_incoming = Complex(strands=[incoming],structure="..") 

The stop conditions are specified as follows:
::

    # creates a "succcessful displacement" complex. This is the incumbent strand forming a complex of its own which means it has been displaced.
    complete_complex_success = Complex(strands=[incumbent],structure=".")
    success_stop_condition = StopCondition("SUCCESS",[(complete_complex_success,Dissoc_Macrostate,0)])

    # A failed displacement stop condition; incumbent falls off.   
    failed_complex = Complex(strands = [incoming], structure="..")  
    failed_stop_condition = StopCondition("FAILURE",[(failed_complex,Dissoc_Macrostate,0)]) 

The simulation itself is started as follows:
::    
    o = Options(simulation_mode="First Step", parameter_type="Nupack", substrate_type="DNA", 
                rate_method = rate_method_k_or_m, num_simulations = num_traj, simulation_time=10.0,  
                dangles = "Some", temperature = 25 + 273.15, rate_scaling = "Calibrated", verbosity = 0)

    o.start_state = [start_complex_incoming, start_complex_substrate_incumbent]
    o.stop_conditions = [success_stop_condition,failed_stop_condition]
    return o


Modifications
-------------

To set up Boltzmann sampling for a complex <start>, we write: 
::

  start[0].boltzmann_sample = True

The boltzmann sampling will be more efficient if we also specify how many samples are needed: 
::

  start[0].boltzmann_count = 100
