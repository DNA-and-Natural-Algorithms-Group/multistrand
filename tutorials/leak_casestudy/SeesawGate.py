from multistrand.objects import Complex, Domain, Strand


class NormalSeesawGate(object):
    # S1, S2, S5, S7, T
    def __init__(self, input_domain, base_domain, output_domain, fuel_domain,
                 toehold_domain):
        # Use the convention of always adding 5' to 3'
        # Setup stuff for this type of gate
        self.input_strand = base_domain + toehold_domain + input_domain
        self.fuel_strand = fuel_domain + toehold_domain + base_domain
        self.base_strand = toehold_domain.C + base_domain.C + toehold_domain.C
        self.output_strand = output_domain + toehold_domain + base_domain
        self.input_partial = Domain(name="partial",
                                    sequence=input_domain.sequence[:3])
        self.threshold_base = self.input_partial.C + toehold_domain.C + \
            base_domain.C
        base_dom_strand = Strand(name="base strand", domains=[base_domain])
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
                                                  base_dom_strand],
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


# MS: Please note that placing mismatches in double stranded regions is
#     currently unsupported
class MismatchedSeesawGate(NormalSeesawGate):
    input_mismatch = 0
    base_mismatch = 1
    output_mismatch = 2
    fuel_mismatch = 3
    toehold_mismatch = 3

    def __init__(self, input_domain, base_domain, output_domain, fuel_domain,
                 toehold_domain, mismatch_domain=None, mismatch_position=None,
                 mismatch_type='C'):
        # Setup intialially as a normal seesaw gate
        super(MismatchedSeesawGate, self).__init__(input_domain, base_domain,
                                                   output_domain, fuel_domain,
                                                   toehold_domain)

        # Then modify for the mismatch itself
        if (mismatch_domain != None):

            if(mismatch_position == MismatchedSeesawGate.input_mismatch):
                print("Placing a mismatch within the input domain")
                self.placeMismatchInInput(
                    mismatch_position, input_domain, base_domain, toehold_domain)
            elif(mismatch_domain == MismatchedSeesawGate.base_mismatch):
                print("It is currently unsuported to place a mismatch within")
                print("a double stranded region. Seesaw will act as normal")
            elif(mismatch_domain == MismatchedSeesawGate.output_mismatch):
                print("Placing a mismatch within the output domain")
                self.placeMismatchInOutput(
                    mismatch_position, output_domain, base_domain, toehold_domain)
            elif(mismatch_domain == MismatchedSeesawGate.fuel_mismatch):
                print("Placing a mismatch within the input domain")
            elif(mismatch_domain == MismatchedSeesawGate.toehold_mismatch):
                print("It is currently unsupported to place a mismatch \n \
                    within a double stranded region, seesaw gate will \n \
                    behave as a normal gate")

        else:
            print "Normal seesaw gate created, but with the option to place mismatches"

     # Utility function for placing mismatches
    def placeMismatchInStrand(self, position, sequence):
        length=len(sequence)
        s=list(sequence)
        if position > (length - 1):
            print "Invalid mismatch position"
        if not s[position] == 'C':
            self.mismatch_type='C'
        else:
            # Don't want it to turn out that 'C' doesn't result in a mismatch!
            # Choose 'G' since most likely to be comparable
            self.mismatch_type='G'
        mis=Domain(name="mismatch", sequence=self.mismatch_type)
        premis=Domain(name="premis", sequence=sequence[:position])
        postmis=Domain(name="postmis", sequence=sequence[position + 1:])
        dom=premis + mis + postmis
        compDom=postmis.C + mis + premis.C
        return dom, compDom

    def placeMismatchInOutput(self, position, output_domain, base_domain, toehold_domain):
        sequence=output_domain.sequence
        mismatched_domain=self.placeMismatchInStrand(position,
                                                       sequence)[0]
        self.output_strand=mismatched_domain + toehold_domain + base_domain
        self.gate_output_complex=Complex(strands=[self.base_strand,
                                                    self.output_strand],
                                           structure=".((+...))")
        self.output_complex=Complex(strands=[self.output_strand],
                                      structure='.' *
                                      len(self.output_strand.sequence))


    def placeMismatchInInput(self, position, input_domain, base_domain, toehold_domain):
        # Exact same as output domains replacement
        sequence=input_domain.sequence
        mismatched_domain=self.placeMismatchInStrand(position,
                                                       sequence)[0]
        self.input_strand=base_domain + toehold_domain + mismatched_domain
        self.input_partial=Domain(name="partial",
                                    sequence=mismatched_domain.sequence[:3])
        self.threshold_base=self.input_partial.C + toehold_domain.C + \
            base_domain.C
        self.gate_input_complex=Complex(strands=[self.base_strand,
                                                   self.input_strand],
                                          structure="((.+)).")
        self.input_complex=Complex(strands=[self.input_strand],
                                     structure='.' *
                                     len(self.input_strand.sequence))
