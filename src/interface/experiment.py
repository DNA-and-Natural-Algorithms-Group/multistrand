
# Often recurring experimental setups
from multistrand.objects import Complex, Domain, Strand, StopCondition
from multistrand.options import Options


def setBoltzmann(complexIn, trials, supersample=1):
    complexIn.boltzmann_supersample = supersample
    complexIn.boltzmann_count = trials
    complexIn.boltzmann_sample = True


# easy handle for options creation
def standardOptions(simMode=Options.firstStep, tempIn=25.0, trials=10, timeOut=0.1):

    output = Options(simulation_mode=simMode,
                     num_simulations=trials,
                     simulation_time=timeOut,
                     temperature=tempIn
                     )

    output.DNA23Metropolis()
    output.rate_method = Options.metropolis

    return output


def makeComplex(seq, dotparen):

    strandList = []

    for seq in seq:

        onedomain = Domain(
            name="domain" + str(makeComplex.counter), sequence=seq)
        makeComplex.counter += 1

        onestrand = Strand(domains=[onedomain])
        strandList.append(onestrand)

    return Complex(strands=strandList, structure=dotparen)


makeComplex.counter = 0



def leakInvasion(options, mySeq, myTrials=0, doFirstPassage=False):
    
    onedomain = Domain(name="d1", sequence=mySeq)
    top = Strand(name="top", domains=[onedomain])
    bot = top.C

    initial_complex = Complex(strands=[top, bot], structure="(+)")
    invader_top = Complex(strands=[top], structure=".")

    # Turns Boltzmann sampling on for this complex and also does sampling more efficiently by sampling 'trials' states.
    if(myTrials > 0):
        setBoltzmann(initial_complex, myTrials)
        setBoltzmann(invader_top, myTrials)

    # Stop when the exact full duplex is achieved.
    stopSuccess = StopCondition(Options.STR_SUCCESS, [(
        success_complex, Options.exactMacrostate, 0)])

    # Declare the simulation unproductive if the strands become single-stranded again.
    failed_complex = Complex(strands=[top], structure=".")
    stopFailed = StopCondition(Options.STR_FAILURE, [(
        failed_complex, Options.dissocMacrostate, 0)])

    options.start_state = [startTop, startBot]

    # Point the options to the right objects
    if not doFirstPassage:

        options.stop_conditions = [stopSuccess, stopFailed]

    else:

        options.stop_conditions = [stopSuccess]
    
    



def hybridization(options, mySeq, myTrials=0, doFirstPassage=False):

    # Using domain representation makes it easier to write secondary structures.
    onedomain = Domain(name="itall", sequence=mySeq)
    top = Strand(name="top", domains=[onedomain])
    bot = top.C

    # Note that the structure is specified to be single stranded, but this will be over-ridden when Boltzmann sampling is turned on.
    startTop = Complex(strands=[top], structure=".")
    startBot = Complex(strands=[bot], structure=".")

    # Turns Boltzmann sampling on for this complex and also does sampling more efficiently by sampling 'trials' states.
    if(myTrials > 0):
        setBoltzmann(startTop, myTrials)
        setBoltzmann(startBot, myTrials)

    # Stop when the exact full duplex is achieved.
    success_complex = Complex(strands=[top, bot], structure="(+)")
    stopSuccess = StopCondition(Options.STR_SUCCESS, [(
        success_complex, Options.exactMacrostate, 0)])

    # Declare the simulation unproductive if the strands become single-stranded again.
    failed_complex = Complex(strands=[top], structure=".")
    stopFailed = StopCondition(Options.STR_FAILURE, [(
        failed_complex, Options.dissocMacrostate, 0)])

    options.start_state = [startTop, startBot]

    # Point the options to the right objects
    if not doFirstPassage:

        options.stop_conditions = [stopSuccess, stopFailed]

    else:

        options.stop_conditions = [stopSuccess]

 # Does not inherit NormalSeesawGate - basically have to rewrite everything

# simulates until a dissociation event occurs.
# Note: for long sequences, this could result in a timeout.


def dissociation(options, mySeq, myTrials=0):

    # Using domain representation makes it easier to write secondary structures.
    onedomain = Domain(name="itall", sequence=mySeq)
    top = Strand(name="top", domains=[onedomain])
    bot = top.C

    # Note that the structure is specified to be single stranded, but this will be over-ridden when Boltzmann sampling is turned on.
    duplex = Complex(strands=[top, bot], structure="(+)")

    # Turns Boltzmann sampling on for this complex and also does sampling more efficiently by sampling 'trials' states.
    if(myTrials > 0):
        setBoltzmann(duplex, myTrials)

    # Stop when the strands fall apart.
    successComplex = Complex(strands=[top], structure=".")
    stopSuccess = StopCondition(Options.STR_SUCCESS, [(
        successComplex, Options.dissocMacrostate, 0)])

    options.start_state = [duplex]
    options.stop_conditions = [stopSuccess]


#
def threewayDisplacement(options, toeholdSeq, domainSeq, doFirstPassage=False, myTrials=0, mySuperSample=1):
    
    toeholdDomain = Domain(name="toehold", sequence=toeholdSeq)
    disDomain = Domain(name="displacement", sequence=domainSeq)
        
    topS = Strand(name="top", domains=[disDomain])
    invaderS = Strand(name="invader", domains=[disDomain, toeholdDomain])
    botS = invaderS.C
    
    startComplex = Complex(strands=[topS, botS], structure="(+.)")
    invaderComplex = Complex(strands=[invaderS], structure="..")
    successComplex = Complex(strands=[botS, invaderS], structure="((+))")
    
    if(myTrials > 0):
        setBoltzmann(startComplex, myTrials, mySuperSample)
        setBoltzmann(invaderComplex, myTrials, mySuperSample)
        
    # stop when the invasion is complete, or when the invader dissociates
    stopSuccess = StopCondition(Options.STR_SUCCESS, [(successComplex, Options.dissocMacrostate, 0)])
    
    # Declare the simulation unproductive if the invader becomes single-stranded again.
    stopFailed = StopCondition(Options.STR_FAILURE, [(invaderComplex, Options.dissocMacrostate, 0)])
    
    # set the starting and stopping conditions
    options.start_state = [startComplex, invaderComplex]

    if doFirstPassage:
        options.stop_conditions = [stopSuccess]
    else:
        options.stop_conditions = [stopFailed, stopSuccess]

# # Strandlist -- stem domain is first argument
# # Strandlist -- loop domain is second argument
# # Try:   ['ACGTGGA', 'TTTTT']

# def hairpinclosing(options,strands_list, myTrials=0):
#     return hairpinclosing(options, strands_list[0], strands_list[1], myTrials=0);
    

def hairpinclosing(options, stemSeq, loopSeq, myTrials=0):

    # Using domain representation makes it easier to write secondary structures.
    stemdomain1 = Domain(name="stemdomain1", sequence=stemSeq)
    loopdomain = Domain(name="loopdomain", sequence=loopSeq)
    stemdomain2 = stemdomain1.C
    
    strand = Strand(name="top", domains=[stemdomain1, loopdomain, stemdomain2])
    start_complex = Complex(strands=[strand], structure="...")
    success_complex = Complex(strands=[strand], structure="(.)")

    # N.B.: myTrials input signature is considered "default", 
    # but in no circumstance will we enable Boltzmann sampling

    # Stop when the exact full duplex is achieved.
    stopSuccess = StopCondition(Options.STR_SUCCESS, [(
        success_complex, Options.exactMacrostate, 0)])

    options.start_state = [start_complex]
    options.stop_conditions = [stopSuccess]


def hairpinopening(options, stemSeq, loopSeq, myTrials=0):
    
    # Using domain representation makes it easier to write secondary structures.
    stemdomain1 = Domain(name="stemdomain1", sequence=stemSeq)
    loopdomain = Domain(name="loopdomain", sequence=loopSeq)
    stemdomain2 = stemdomain1.C
    
    strand = Strand(name="top", domains=[stemdomain1, loopdomain, stemdomain2])
    start_complex = Complex(strands=[strand], structure="(.)")
    success_complex = Complex(strands=[strand], structure="...")

    # N.B.: myTrials input signature is considered "default", 
    # but in no circumstance will we enable Boltzmann sampling

    # Stop when the exact full duplex is achieved.
    stopSuccess = StopCondition(Options.STR_SUCCESS, [(
        success_complex, Options.exactMacrostate, 0)])

    options.start_state = [start_complex]
    options.stop_conditions = [stopSuccess]


# Aug 2017: this is Mrinanks' implementation of the clamped seesaw gate.
# this is how much the domain overlaps
# see Figure 1, Winfree Qian 2010 -- A simple DNA gate motif for synthesizing large-scale circuits
SEESAW_DELTA = 5


# Domain list as a list of strings, defined as follows input_sequence, base_sequence, output_sequence, fuel_sequence,
# toehold_sequence, clamp_sequence which defaults to clamp_sequence
# This method takes a gate output complex and its input and then calculate the rate of output production
def seesaw_gate_output_production(options, gate, trials, supersample=25, doFirstPassage=False):
    two_input(options, gate.gate_output_complex, gate.input_complex,
              gate.output_complex, trials, supersample, doFirstPassage=doFirstPassage)

# Domain list as a list of strings, defined as follows input_sequence, base_sequence, output_sequence, fuel_sequence,
# toehold_sequence, clamp_sequence which defaults to clamp_sequence
# This method takes a gate input complex and its fuels and then calculates the rate of input regeneration


def seesaw_gate_fuel_catalysis(options, gate, trials, supersample=25, doFirstPassage=False):
    two_input(options, gate.gate_input_complex, gate.fuel_complex,
              gate.input_complex, trials, supersample, doFirstPassage=doFirstPassage)


# Domain list as a list of strings, defined as follows input_sequence, base_sequence, output_sequence, fuel_sequence,
# toehold_sequence, clamp_sequence which defaults to clamp_sequence
# This method takes a gate output complex with its fuel complex and calculates the ***leak*** rate at which
# the fuel displaces the output i.e. the rate of leak output production.
def seesaw_gate_fuel_leak(options, gate, trials, supersample=25, doFirstPassage=False):
    two_input(options, gate.gate_output_complex, gate.fuel_complex,
              gate.output_complex, trials, supersample, doFirstPassage=doFirstPassage)


# Takes two gate objects and calcualtes their leak!
def seesaw_gate_gate_leak(options, gateA, gateB, trials, supersample=25, doFirstPassage=False):
    two_input_two_success(options, gateA.gate_output_complex, gateB.gate_output_complex,
                          gateA.output_complex, gateB.output_complex, trials, supersample, doFirstPassage)


def two_input(options, input_complex_A, input_complex_B, output_complex, trials=0, supersample=25, doFirstPassage=False):

    if(trials > 0):
        for x in [input_complex_A, input_complex_B]:
            setBoltzmann(x, trials, supersample)

    successful_stop_condition = StopCondition(
        Options.STR_SUCCESS, [(output_complex, Options.dissocMacrostate, 0)])
    failure_stop_condition = StopCondition(
        Options.STR_FAILURE, [(input_complex_B, Options.dissocMacrostate, 0)])

    options.start_state = [input_complex_A, input_complex_B]

    # Point the options to the right objects
    if not doFirstPassage:
        options.stop_conditions = [
            successful_stop_condition, failure_stop_condition]
    else:
        options.stop_conditions = [successful_stop_condition]


def two_input_two_success(trials, options, input_complex_A, input_complex_B, output_complex_A, output_complex_B, supersample=25, doFirstPassage=False):

    if(trials > 0):
        for x in [input_complex_A, input_complex_B]:
            setBoltzmann(x, trials, supersample)

    successful_stop_condition = StopCondition(
        Options.STR_SUCCESS, [(output_complex_A, Options.dissocMacrostate, 0)])
    alt_successful_stop_condition = StopCondition(
        Options.STR_ALT_SUCCESS, [(output_complex_B, Options.dissocMacrostate, 0)])
    failure_stop_condition = StopCondition(
        Options.STR_FAILURE, [(input_complex_B, Options.dissocMacrostate, 0)])

    options.start_state = [input_complex_A, input_complex_B]

    # Point the options to the right objects
    if not doFirstPassage:
        options.stop_conditions = [successful_stop_condition,
                                   alt_successful_stop_condition, failure_stop_condition]
    else:
        options.stop_conditions = [
            successful_stop_condition, alt_successful_stop_condition]


class ClampedSeesawGate(object):

    Gate_Count = 1

    def __init__(self, input_sequence, base_sequence, output_sequence, fuel_sequence,
                 toehold_sequence, clamp_sequence="CG", sameID=False):

        count_str = str(ClampedSeesawGate.Gate_Count) + '_Cl '
        self.input_domain = Domain(
            name="input_domain_" + count_str, sequence=input_sequence)
        self.base_domain = Domain(
            name="base_domain_" + count_str, sequence=base_sequence)
        self.output_domain = Domain(
            name="output_domain_" + count_str, sequence=output_sequence)
        self.fuel_domain = Domain(
            name="fuel_domain_" + count_str, sequence=fuel_sequence)
        self.toehold_domain = Domain(
            name="toehold_domain_" + count_str, sequence=toehold_sequence)
        self.clamp_domain = Domain(
            name="clamp_domain_" + count_str, sequence=clamp_sequence)

        # Use the convention of always adding 5' to 3'
        # Setup stuff for this type of gate

        # Clamp domain setup - add clamp domains either side of each recognition domain
        self.input_strand = self.clamp_domain + self.base_domain + self.clamp_domain + \
            self.toehold_domain + self.clamp_domain + self.input_domain + self.clamp_domain

        self.fuel_strand = self.clamp_domain + self.fuel_domain + self.clamp_domain + \
            self.toehold_domain + self.clamp_domain + self.base_domain + self.clamp_domain

        self.base_strand = self.clamp_domain.C + self.toehold_domain.C + self.clamp_domain.C + \
            self.base_domain.C + self.clamp_domain.C + \
            self.toehold_domain.C + self.clamp_domain.C

        self.output_strand = self.clamp_domain + self.output_domain + self.clamp_domain + \
            self.toehold_domain + self.clamp_domain + self.base_domain + self.clamp_domain

        self.input_partial = Domain(name="partial",
                                    sequence=self.input_domain.sequence[:SEESAW_DELTA])

        self.threshold_base = self.input_partial.C + self.clamp_domain.C + \
            self.toehold_domain.C + self.clamp_domain + \
            self.base_domain.C + self.clamp_domain
        self.base_dom_strand = self.clamp_domain + self.base_domain + self.clamp_domain
        
        self.input_strand.name = "input"
        self.fuel_strand.name = "fuel"
        self.base_strand.name = "base" 
        self.output_strand.name = "output" 
        self.input_partial.name = "inputP" 
        self.threshold_base.name = "thres"
        self.base_dom_strand.name = "baseD"

        self.gate_output_complex = Complex(strands=[self.base_strand,
                                                    self.output_strand],
                                           structure="..(((((+..)))))")
        self.gate_fuel_complex = Complex(strands=[self.base_strand,
                                                  self.fuel_strand],
                                         structure="..(((((+..)))))")
        self.gate_input_complex = Complex(strands=[self.base_strand,
                                                   self.input_strand],
                                          structure="(((((..+)))))..")
        self.threshold_complex = Complex(strands=[self.threshold_base,
                                                  self.base_dom_strand],
                                         structure="...(((+)))")
        self.input_complex = Complex(strands=[self.input_strand],
                                     structure='.' * 
                                     len(self.input_strand.sequence))
        self.fuel_complex = Complex(strands=[self.fuel_strand],
                                    structure='.' * 
                                    len(self.fuel_strand.sequence))
        self.output_complex = Complex(strands=[self.output_strand],
                                      structure='.' * 
                                      len(self.output_strand.sequence))

        if not sameID :
            ClampedSeesawGate.Gate_Count += 1
