# nupack_calls.py
#
# Multistrand is equipped with Python wrapper functions for many of Nupack's executables, including
# pfunc, energy, defect, prob, sample, mfe, and pairs. Here we provide usage examples for all
# provided wrapper functions. See the Nupack User Manual for more information about the arguments to
# each executable. Many of these examples are adapted from examples provided with Nupack.
# Note that these functions assume you have Nupack version 3.0.2 or higher. Certain functions, such
# as sample(), will not work with earlier versions of Nupack.
#
# In general, the following arguments are required:
#   sequences  - a list of DNA/RNA sequences
#   structure  - a dot-paren structure (required by energy, prob, and defect)
#   samples    - an integer specifying how many samples to draw (required by sample)
#
# The following optional arguments may also be included in all function calls:
#   ordering   - a list of integer indices indicating the ordering of the given strand sequences
#                in the complex, or None if each strand is distinguishable, occurs exactly once,
#                and is ordered as in the sequences argument. Indices may be repeated to indicate
#                identical, indistinguishable strands. This argument is ignored if multi == False.
#                (default: None)
#   material   - a string representing the material parameters to use (e.g. 'rna' or 'dna')
#                The Nupack User Manual has a complete list of possible values.
#                (default: 'rna')
#   multi      - whether or not to allow multistrand structures with identical strands
#                only works with DNA (default: True)
#   pseudo     - whether or not to allow pseudoknotted structures
#                only works with RNA (default: False)
#   dangles    - 'none', 'all', or 'some' to specify what types of dangle corrections to use (default: 'some')
#   T          - temperature of the system in degrees C (default: 37)
#   sodium     - concentration of Na+ in molar (default: 1.0)
#   magnesium  - concentration of Mg++ in molar (default: 0.0)
#
# The following optional arguments are allowed by only some functions:
#   cutoff     - (for use with pairs) probabilities less than cutoff are excluded (default: 0.001)
#   degenerate - (for use with mfe) degenerate structures are all returned (default: False)
#   mfe        - (for use with defect) returns mfe defect rather than ensemble defect (default: False)


## Sequences used throughout this file:
rna_seq = ['GGGCUGUUUUUCUCGCUGACUUUCAGCCCCAAACAAAAAAUGUCAGCA'];
dna_seqs = ['AGTCTAGGATTCGGCGTGGGTTAA',
            'TTAACCCACGCCGAATCCTAGACTCAAAGTAGTCTAGGATTCGGCGTG',
            'AGTCTAGGATTCGGCGTGGGTTAACACGCCGAATCCTAGACTACTTTG']
dna_struct = '((((((((((((((((((((((((+))))))))))))))))))))))))((((((((((((.(((((((((((+........................))))))))))).))))))))))))'


## Examples with pfunc()
# Find the partition function of this set of three DNA strands.
pfunc(dna_seqs, material = 'dna')
# Output: -62.66453454292297

# Find the partition function of this single RNA strand, including pseudoknotted structures.
pfunc(rna_seq, material = 'rna', multi = False, pseudo = True)
# Output: -17.2786076


## Example with energy()
# Find the energy (in kcal/mol) of this DNA structure.
energy(dna_seqs, dna_struct, material = 'dna')
# Output: -61.79270283815892

## Example with defect()
# Find the ensemble defect of this DNA structure.
defect(dna_seqs, dna_struct, material = 'dna')
# Output: 8.297

# Find the mfe defect of this DNA structure.
defect(dna_seqs, dna_struct, material = 'dna', mfe = True)
# Output: 2.0

## Examples with prob()
# Find the probability of this DNA structure.
prob(dna_seqs, dna_struct, material = 'dna')
# Output: 7.992e-05

## Examples with sample()
# Boltzmann sample 3 structures with these DNA strands.
sample(dna_seqs, 3, material = 'dna')
# Example output:
#  ['.(((((((((((((((((((((((+))))))))))))))))))))))).((((((((((((((((((((((((+..............(...).....))))))))))))))))))))))))',
#   '(((((((((((((((((((((((.+.))))))))))))))))))))))).(((((.(((((((((((((((((+........................))))))))))))))))).))))).',
#   '((((((((((((((((((((((..+..))))))))))))))))))))))..((((((((((((((((((((((+........................))))))))))))))))))))))..']

## Examples with mfe()
# Find the MFE structure of this ordering of DNA sequences (the energy is also provided).
mfe(dna_seqs, material = 'dna')
# Output:
#  [('((((((((((((((((((((((((+))))))))))))))))))))))))((((((((((((((((((((((((+........................))))))))))))))))))))))))',
#   '-65.913')]

## Examples with pairs()
# Find the pair probabilities of the following unpseudoknotted DNA strands.
# The output is a list of tuples of the form (i, j, p) where p is the
# probability that the ith and jth bases are bound.
pairs(['ACGT','GCTT'], material = 'dna')
# Output:
#  [('1', '7', '3.1271e-02'),
#   ('1', '8', '4.9495e-02'),
#   ('2', '5', '2.1630e-01'),
#   ('3', '6', '7.3677e-01'),
#   ('1', '9', '9.1923e-01'),
#   ('2', '9', '7.8370e-01'),
#   ('3', '9', '2.6323e-01'),
#   ('4', '9', '1.0000e+00'),
#   ('5', '9', '7.8370e-01'),
#   ('6', '9', '2.6323e-01'),
#   ('7', '9', '9.6873e-01'),
#   ('8', '9', '9.5051e-01')]
