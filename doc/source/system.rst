Multistrand Core
================

.. module:: multistrand.system
   :synopsis: Multistrand's core simulation components

.. class:: multistrand.system.SimSystem( self, *args )
      
   Python Wrapper for Multistrand's C++ SimulationSystem object.

   Provides a very simple interface to the SimulationSystem's 
   StartSimulation method. 

   .. method:: __init__(self, options )

      Initializes the SimSystem object with the given options.
      
      :param options: Defines the simulation options to be used.
      :type options: :class:`Options <multistrand.options.Options>`

   .. method:: start(self)

      Starts the simulation, only returns once the simulation has
      finished.

   .. attribute:: options

      The :class:`Options <multistrand.options.Options>` object that
      initialized this simulation.


.. function:: energy( start_state [, options = None, energy_type = 0])

   Compute the energy of the input state using the energy model
   settings, either from the options argument or the default.

   If no options object is given, uses the energymodel already
   initialized via :func:`initialize_energy_model
   <multistrand.system.initialize_energy_model>` or via a SimSystem call.

   :param start_state: A list of :class:`Complex <multistrand.objects.Complex>` or :class:`RestingState <multistrand.objects.RestingState>` objects.
   :param options: - With the default :obj:`None`, assume that the energy model
                     has already been initialized and raise an error
                     if there was energy model found.

                   - Otherwise, this should be an :class:`Options <multistrand.options.Options>`
                     object. 

                   *IMPORTANT NOTE*: If an energy model has
                   already been used, either by creating a
                   :class:`~multistrand.system.SimSystem` object or
                   via :func:`initialize_energy_model <multistrand.system.initialize_energy_model>`,
                   this will not overwrite that model; if you wish
                   to do so, you should first do
                   :func:`initialize_energy_model(None) <multistrand.system.initialize_energy_model>`
                   to explicitly clear out the model. See that
                   function for more details.

   :param energy_type: 0. do not include dG\ :sub:`volume` or dG\ :sub:`assoc` [Vienna default]
                       1. include dG\ :sub:`volume`
                       2. include dG\ :sub:`assoc` [NUPACK default]
                       3. include dG\ :sub:`volume` + dG\ :sub:`assoc`

                       Note that these are all equivalent for a single
                       stranded system. Which default is used depends
                       on the energy model chosen; The Vienna
                       parameter set does not include a dG_assoc term.

.. function:: initialize_energy_model( options = None )
   
   Initialize the Multistrand module's energy model using the options
   object given. If a model already exists, this will remove the old
   model and create a new one.

   **WARNING**\ :This changes out the energy model any existing
   SimSystem objects are using, so if you use it in an interactive
   session you should be very careful that you mean to do so and have
   no SimSystem objects running.

   This function is NOT required to use other parts of the module;
   they will never cause the energy model to be changed out if one
   already exists. This odd behaviour of energy models is a design
   defect and will likely be fixed in the future, so that multiple
   energy models can co-exist as necessary. 
   
   :param options: When present, replaces the existing energy model
                   with the one defined by the given 
                   :class:`~multistrand.options.Options`
                   object. If :obj:`None`, this removes the old energy
                   model and does not create a new one.

.. function:: run_system( options )
   
   Run a simulation.

   This is a shortcut for creating a SimSystem object and
   then calling the .start() method; 

   :param options: The simulation to run.
   :type options: :class:`Options <multistrand.options.Options>`
   :rtype: :obj:`None`
