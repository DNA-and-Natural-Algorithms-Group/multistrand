Input File Format Samples
=========================

Multistrand 1.0 style input file
--------------------------------

:: 
 
  #Strands
  S,AGTACGGACACTAGCTGGATCTGAGGATTAGT
  Q,ACTAATCCTCAGATCCAGCTAGTGTC
  T6,ACTAATCCTCAGATCCAGCTAGTGTCCGTACT
  #StartStructure
  S,Q
  ......((((((((((((((((((((((((((_))))))))))))))))))))))))))
  T6
  ................................
  #SimTime=2000000
  #NumSims=1
  #OutputInterval=0
  #StopOptions=2
  #Energymodel=NUPACK_DNA_2_3
  #Temperature = 20.0
  #Concentration = .001
  #StopStructures
  Q
  ..........................
  TAG: finalstop
  ##
  ## end old file

Multistrand 2.0 style 'input' file
----------------------------------

::

   start = [Complex("Starting Gate", ["S","Q"], ["AGTACGGACACTAGCTGGATCTGAGGATTAGT",
                                                 "ACTAATCCTCAGATCCAGCTAGTGTC"],
                                                ["......((((((((((((((((((((((((((",
                                                 "))))))))))))))))))))))))))"]),
            Complex("Input","T6","ACTAATCCTCAGATCCAGCTAGTGTCCGTACT")]
   stop = StopCondition("Output","Q")
   my_options = Options(sim_time=2e6,
                        num_sims=1,
                        parameter_type='Nupack',
                        substrate_type='DNA',
                        temperature=20.0,
                        concentration=1e-6,  # 1 microM, equivalent to the
                                             # old which was in units of mM
                        start_structure = start,
                        stop_conditions = [stop])  #or we can use `stop_condition=stop`

   run_system( my_options ) # run_system by default uses the
                           #  standard logging/interfacing settings.

Commentary
----------

Both of these styles are supported in 2.0, though we have yet to exactly define how a new style file would be run; it likely would be something like "multistrand newstylefile.py", which implies some things about what the 'multistrand' executable looks like that we haven't finalized in the new format.

Modifications
-------------

The new style is pretty easy to modify, for example we may want to do a first step simulation after our original normal mode sim, so we add this:

:: 

  forward = StopCondition("Forward Stop State","Q")
  reverse = StopCondition("Reverse Stop State",start, tolerance=2)

  firststep_options = my_options.copy()
  firststep_options.stop_conditions = [forward,reverse]
  firststep_options.simulation_mode = 'First Step'
  run_system( firststep_options )

Similarly, adding Boltzmann sampling to that is easy: ::

  start[0].boltzmann_sample = True

Possibly if we know we'd need a lot of samples, we should do this as well, to give the simulation a hint about how many samples may be needed: ::

  start[0].boltzmann_count = 100

Note that we don't have to update anything else if we want to use the same stop conditions in this case, they're already referring to the right object.

This format for setting up objects is pretty flexible, we could instead do: ::

  start = [Complex("Starting Gate",
                 [['S','AGTACGGACACTAGCTGGATCTGAGGATTAGT','......(((((((((((((((((((((((((('],
                  ['Q', 'ACTAATCCTCAGATCCAGCTAGTGTC', '))))))))))))))))))))))))))']]),
           Complex("Input","T6","ACTAATCCTCAGATCCAGCTAGTGTCCGTACT")]

Or even: ::

  strand_0 = Strand('S','AGTACGGACACTAGCTGGATCTGAGGATTAGT')
  strand_1 = Strand('Q', 'ACTAATCCTCAGATCCAGCTAGTGTC')
  start = [Complex("Starting Gate", [strand_0, strand_1],
                           '......((((((((((((((((((((((((((+))))))))))))))))))))))))))'),
           Complex("Input","T6","ACTAATCCTCAGATCCAGCTAGTGTCCGTACT")]

Amusingly enough, we can also create a input file straight from our interactive session, so if we happened to be messing around with stuff until we find the right set of options, we can store it for later use!

  >>> print( repr(my_options) )
  Options( join_concentration=1e-06,
           simulation_time=2000000.0,
           temperature=293.15,
           start_state=<state>,
           stop_conditions=<conditions>)
  >>> my_options.to_file('myoptions.py')

The keen eye among you may have noted that this temperature is not the
20 degrees C we had originally stated! Cross referencing the Options
object's :attr:`temperature <multistrand.options.Options.temperature>`, we see
that it is always stored in degrees K internally, and if it does any
conversion there's a message in the error log. Let's check that:

::

  >>> o = Options()
  >>> o.temperature = 20.0
  >>> print o.temperature
  300.15
  >>> print o.errorlog[-1]
  Warning: Temperature was set at the value [27]. We expected a value in Kelvin, or with appropriate units.
           Temperature was automatically converted to [300.15] degrees Kelvin.

See the above link to the parameter for details on how/when it
converts.

Notes
~~~~~

1. Some of this is in the <future> - right now there's no repr for states/conditions, but it's on my list.
#. Options that were equivalent to defaults were not shown in the 'repr' example above. So rather than the spam of every single possible option here, you get only the ones that were meaningful. There's currently on the order of 20-ish options you could have set, so printing out every one isn't usually very useful (though possible, if you want to see it).
#. None of this new stuff is really set in stone, and it's easy to manipulate; the old style file was a pain in the neck for adding new options. The new style it's easy, and the best part is that nearly all of the error checking is actually Python's problem, rather than my trying to catch every possible syntax error.
#. :func:`run_system <multistrand.system.run_system>` could be implicit in a python version input file, and we could easily define that style as having no variable definitions at all. Supporting that would be simple as we'd just add a dictionary off all strands, etc, created via input files and the user can just refer to them with whatever names they used.
#. Some of the names I used were edited slightly between the two examples - the TAG, additional state names, etc. The old style had some issues with spaces in identifiers, but it's no longer a worry as they must be distinct strings anyways.
#. Similarly, in some cases I didn't use a parameter - e.g. the stop state which is completely open; 'reasonable' defaults are used in cases like this, so when structure isn't specified for a complex, it is always completely unpaired. There's an error if you try to do so when the complex has more than one strand.
