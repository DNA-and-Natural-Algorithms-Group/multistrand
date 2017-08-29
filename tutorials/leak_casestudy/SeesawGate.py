from multistrand.objects import Complex, Domain, Strand
from multistrand.concurrent import Bootstrap, FirstStepRate

MISMATCH_TYPE = 'C'

SEESAW_DELTA = 5


class NormalSeesawGate(object):
    Gate_Count = 1

    # S1, S2, S5, S7, T
    def __init__(self, input_sequence, base_sequence, output_sequence, fuel_sequence,
                 toehold_sequence):

        count_str = str(NormalSeesawGate.Gate_Count)
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

        # Use the convention of always adding 5' to 3'
        # Setup stuff for this type of gate
        self.input_strand = self.base_domain + self.toehold_domain + self.input_domain
        self.fuel_strand = self.fuel_domain + self.toehold_domain + self.base_domain
        self.base_strand = self.toehold_domain.C + \
            self.base_domain.C + self.toehold_domain.C
        self.output_strand = self.output_domain + \
            self.toehold_domain + self.base_domain
        self.input_partial = Domain(name="partial",
                                    sequence=self.input_domain.sequence[:SEESAW_DELTA])
        self.threshold_base = self.input_partial.C + self.toehold_domain.C + \
            self.base_domain.C
        self.base_dom_strand = Strand(
            name="base strand", domains=[self.base_domain])
        self.threshold_free_waste_complex = Complex(
            strands=[self.base_dom_strand], structure='.' * len(self.base_dom_strand.sequence))

        self.gate_output_complex = Complex(strands=[self.base_strand,
                                                    self.output_strand],
                                           structure=".((+.))")
        self.gate_fuel_complex = Complex(strands=[self.base_strand,
                                                  self.fuel_strand],
                                         structure=".((+.))")
        self.gate_input_complex = Complex(strands=[self.base_strand,
                                                   self.input_strand],
                                          structure="((.+)).")
        self.threshold_complex = Complex(strands=[self.threshold_base,
                                                  self.base_dom_strand],
                                         structure="..(+)")
        self.input_complex = Complex(strands=[self.input_strand],
                                     structure='.' *
                                     len(self.input_strand.sequence))
        self.fuel_complex = Complex(strands=[self.fuel_strand],
                                    structure='.' *
                                    len(self.fuel_strand.sequence))
        self.output_complex = Complex(strands=[self.output_strand],
                                      structure='.' *
                                      len(self.output_strand.sequence))
        self.threshold_complex_output_occluded = Complex(strands=[self.threshold_base,
                                                  self.base_dom_strand, self.output_strand],
                                         structure=".((+)+.).")
        self.gate_input_complex_input_occluded = Complex(strands=[self.base_strand,
                                                   self.input_strand, self.input_strand],
                                          structure="(((+.).+)).")     
        self.gate_output_complex_output_occluded = Complex(strands=[self.base_strand,
                                                    self.output_strand, self.output_strand],
                                           structure="(((+.))+.).")                                                              
        NormalSeesawGate.Gate_Count += 1



class ClampedSeesawGate(object):

    Gate_Count = 1

    def __init__(self, input_sequence, base_sequence, output_sequence, fuel_sequence,
                 toehold_sequence, clamp_sequence="CG"):

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
        self.threshold_complex_output_occluded = Complex(strands=[self.threshold_base,
                                                  self.base_dom_strand, self.output_strand],
                                         structure=".(((((+)))+...))..")
        self.gate_input_complex_input_occluded = Complex(strands=[self.base_strand,
                                                   self.input_strand, self.input_strand],
                                           structure="(((((((+..))...+)))))..")
        self.gate_output_complex_output_occluded = Complex(strands=[self.base_strand,
                                                    self.output_strand, self.output_strand],
                                           structure="(((((((+..)))))+...))..")      
        ClampedSeesawGate.Gate_Count += 1 
    
# MS: Please note that placing mismatches in double stranded regions is
#     currently unsupported
# 
#  Please note that the occulsion reactions are currently unsuported with this
#  type of SeesawGate.
# 
class MismatchedSeesawGate(NormalSeesawGate):

    # Utility function for placing mismatches
    def placeMismatchInDomain(self, position, sequence):
        pos = position
        length = len(sequence)
        s = list(sequence)
        if position > (length - 2) or position < 1:
            print "Invalid mismatch position - placing mismatch at position 1"
            pos = 1
        self.mismatch_type = MISMATCH_TYPE
        mis = Domain(name="mismatch", sequence=self.mismatch_type)
        premis = Domain(name="premis", sequence=sequence[:pos])
        postmis = Domain(name="postmis", sequence=sequence[pos + 1:])
        strand = premis + mis + postmis
        strand_prime = postmis.C + mis + premis.C
        return strand, strand_prime, strand.C

    def placeMismatchInOutput(self, position):
        mismatched_strands = self.placeMismatchInDomain(
            position, self.output_domain.sequence)
        self.output_strand = mismatched_strands[0] + \
            self.toehold_domain + self.base_domain
        self.gate_output_complex = Complex(strands=[self.base_strand,
                                                    self.output_strand],
                                           structure=".((+...))")
        self.output_complex = Complex(strands=[self.output_strand],
                                      structure='.' *
                                      len(self.output_strand.sequence))

    def placeMismatchInInputWire(self, position):
        mismatched_strands = self.placeMismatchInDomain(
            position, self.base_domain.sequence)
        # So now the input strand has a 'C' in the designated position
        recog_len = len(self.base_domain.sequence)
        # Place a mismatch in the input wire
        self.input_strand = mismatched_strands[0] + \
            self.toehold_domain + self.input_domain
        # make the new base strand have a 'C-C' mismatch
        self.base_strand = self.toehold_domain.C + \
            mismatched_strands[1] + self.toehold_domain.C
        # we want to keep our the gate
        self.base_domain = mismatched_strands[1].C

        # Update the relevant strands and complexes to take the mismatch into account
        # Use the convention of always adding 5' to 3'
        # Setup stuff for this type of gate
        self.fuel_strand = self.fuel_domain + self.toehold_domain + self.base_domain
        self.output_strand = self.output_domain + \
            self.toehold_domain + self.base_domain

        # We do want the threshold to still act well!!
        self.threshold_base = self.input_partial.C + self.toehold_domain.C + \
            mismatched_strands[0].C
        self.base_dom_strand = Strand(
            name="base strand", domains=[mismatched_strands[0]])

        self._redefineMismatchedComplexes()
        self.gate_input_complex = Complex(strands=[self.base_strand,
                                                   self.input_strand],
                                          structure="((.(.+).)).")

    def placeMismatchInFuelWireBase(self, position):
        # alter the base domain, as per usual
        mismatched_strands = self.placeMismatchInDomain(
            position, self.base_domain.sequence)
        # So now the input strand has a 'C' in the designated position
        # Get the standard recognition domain length
        recog_len = len(self.base_domain.sequence)

        self.fuel_strand = self.fuel_domain + \
            self.toehold_domain + mismatched_strands[0]
        # make the new base strand have a 'C-C' mismatch with the fuel
        self.base_strand = self.toehold_domain.C + \
            mismatched_strands[1] + self.toehold_domain.C
        # we want to keep our the gate - as above, the complement for the prime in most places!
        self.base_domain = mismatched_strands[1].C

        # Update the relevant strands and complexes to take the mismatch into account
        # Use the convention of always adding 5' to 3'
        # Setup stuff for this type of gate
        self.input_strand = self.base_domain + self.toehold_domain + self.input_domain
        self.output_strand = self.output_domain + \
            self.toehold_domain + self.base_domain

        # We do want the threshold to still act well!!
        self.threshold_base = self.input_partial.C + self.toehold_domain.C + \
            mismatched_strands[0].C
        self.base_dom_strand = Strand(
            name="base strand", domains=[mismatched_strands[0]])

        self._redefineMismatchedComplexes()
        self.gate_fuel_complex = Complex(strands=[self.base_strand,
                                                  self.fuel_strand],
                                         structure=".(.((+.)).)")

    # This method redefines all complexes using a base domain comprised of
    # three subdomains rather than just one normally. This method assumes that
    # all complexes are meant to be fully bound, which is clearly not the case given
    # that there is a mismatch in at least one of them!!
    def _redefineMismatchedComplexes(self):
        # Redefine complexes
        # print "Warning: one of these complexes will be invalid! Must replace"
        self.threshold_free_waste_complex = Complex(
            strands=[self.base_dom_strand], structure='.' * len(self.base_dom_strand.sequence))
        self.gate_output_complex = Complex(strands=[self.base_strand,
                                                    self.output_strand],
                                           structure=".((((+.))))")
        self.gate_fuel_complex = Complex(strands=[self.base_strand,
                                                  self.fuel_strand],
                                         structure=".((((+.))))")
        self.gate_input_complex = Complex(strands=[self.base_strand,
                                                   self.input_strand],
                                          structure="((((.+)))).")
        self.threshold_complex = Complex(strands=[self.threshold_base,
                                                  self.base_dom_strand],
                                         structure="..(((+)))")

        self.input_complex = Complex(strands=[self.input_strand],
                                     structure='.' *
                                     len(self.input_strand.sequence))
        self.fuel_complex = Complex(strands=[self.fuel_strand],
                                    structure='.' *
                                    len(self.fuel_strand.sequence))
        self.output_complex = Complex(strands=[self.output_strand],
                                      structure='.' *
                                      len(self.output_strand.sequence))



