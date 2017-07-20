from multistrand.objects import Complex, Domain, Strand


class NormalSeesawGate:
    # S1, S2, S5, S7, T
    def __init__(self, input_domain, base_domain, output_domain, fuel_domain,
                 toehold_domain, type="normal"):
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
        self.output_complex = Complex(strands=[self.output_strand],
                                      structure='.' *
                                      len(self.output_strand.sequence))


class MismatchedSeesawGate:
    def __init__(self, input_domain, base_domain, output_domain, fuel_domain,
                 toehold_domain, ):
        print ("test")