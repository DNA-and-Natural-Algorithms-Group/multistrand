// Multistrand nucleic acid kinetic simulator
// Copyright (c) 2010-2017 California Institute of Technology. All rights reserved.
// The Multistrand Team (help@multistrand.org)

This directory contains a number of commented, runnable Multistrand scripts that
illustrate how to use a variety of simulation modes.  These examples are for
people who want to understand the full power and generality of the Multistrand
simulation environment.  (And even so, they leave much out.)

Although you can just run the python code and see it work, and you can look
at the code and read the comments about how it works, the BEST way to learn
from these examples is to try them out line-by-line in the python interpreter,
taking the time to examine the objects being created -- using 'print' and 
'help' and playing around to see what's there.

It is suggested to examine and run the tutorial scripts in the following order.

tutorial-00-the-NUPACK-interface.py [changed name and edited]
--- Examples of how to call NUPACK directly to get thermodynamic information.

tutorial-01-hairpin-energies.py              
--- See how to make a complex and calculate its energy using Multistrand.
--- Plot a simple energy landscape for hairpin folding.

tutorial-02-hairpin-trajectories.py
--- See how to set up a simulation and how to recover the states visited.

tutorial-03-hairpin-transitions.py 
--- See how to define macrostates, i.e. sets of related secondary structure microstates,
--- so that the simulation doesn't return every elementary step, 
--- but just returns when a macrostate is entered or exited.
--- Tabulate transition frequencies between macrostates, for coarse-grained analysis.
--- Examine the consequences of using different coarse-grainings, e.g. "loose" macrostates.

tutorial-04-hairpin-first-passage-times.py  [changed name and edited]
--- Don't look at trajectories at all, just tabulate how long they took, and compute statistics.
--- See the effects of kinetic traps in folding, by comparing sequences.
--- Use the random number seed to reproduce trajectories of interest.

tutorial-05-threeway-energies.py
--- See how to define complexes with multiple strands, specifically
--- ones that undergo toehold-mediated three-way branch migration and strand displacement.
--- Consider the energy of a test tube, as opposed to the energy of a complex.
--- Plot test tube energy for a hypothetical pathway for toehold-mediated three-way strand displacement.

tutorial-06-threeway-trajectories.py       
--- Run simulations in which strands might dissociate.
--- Look at how to extract information from the trajectories.

tutorial-07-threeway_transitions.py [notebook needs to be created: Frits' job]
--- Use transition mode for multistranded systems.
--- Again, compare exact versus loose macrostates.

tutorial-08-threeway-first-passage-times.py [new, not yet complete]
--- Compare two sequence designs based on 
--- (a) their duration statistics, and (b) their frequencies of reaching either "success" or "failure" states.  

tutorial-09-threeway-first-step-mode.py [needs to be created from hybridization-first-step-mode: Erik's job]
--- This is the preferred method for deducing bimolecular rate constants:
--- it allows more efficient simulation that skips "equilibration" between interactions.
--- Deduce first-order and second-order rate constants from probability of success and time to completion.
--- And compute error bars! (Precision, not accuracy.)

tutorial-10-threeway-multiprocessing.py [needs to be renamed from threeway-first-step-mode: Erik's job]
--- Automate the running of a month's worth of first-step-mode simulations,
--- thus examining the dependence of toehold-mediated strand displacement rate constants on toehold length.

tutorial-11-hybridization-first-step-mode.py [needs to be simplified since material transfers to threeway-hybridization-first-step-mode: Erik's job]
--- Examine the rate with which two single strands hybridize.

tutorial-12-hybridization-comparison.py [notebook needs to be created: Erik's job]
--- Look at two strands interacting in a tiny box.
--- Compare sequences that are unstructured vs have hairpins in various locations.
--- Compute association and dissociation rate constants.
--- Compare using first passage time, transition mode, and first step mode.
--- See how the time to reach a full duplex scales with concentration.

tutorial-13-hybridization-scatterplot.py [notebook needs to be created: Erik's job]
--- Compare the association and dissociation rates for many random sequences.
--- Look for trends related to secondary structure.

tutorial-14-hybridization-transitions.py [notebook needs to be created: Erik's job]
--- Examine coarse-graining methods for hybridization of structured hairpins.
--- Look at "typical" pathways.
--- Attempt to extract dissociation rates, and understand the modeling challenges.

Contributions: The scripts and notebooks were developed by Joseph Schaeffer, Niranjan Srinivas, Xander Rudelis, 
Joseph Berleant, Chris Thachuk, Frits Dannenberg, and Erik Winfree
