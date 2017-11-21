# Erik Winfree and Frits Dannenberg
# May 2017

"""
Call this using arguments 

generate iso-random 15 10 ir15-10
generate iso-random 25 100 ir25-100
generate iso-structured-random 15 5 5 10 sr15-5-5-10
generate iso-structured-random 25 5 4 17 sr25-5-4-17
generate structured-random 25 3 4 10 sr25-10
"""
# from msUtil import myMultistrand


from multistrand.experiment import hybridization, standardOptions
from multistrand.concurrent import FirstStepRate,FirstStepLeakRate,  myMultistrand
from constantsgao import goa2006_P0, goa2006_P3, goa2006_P4, setSaltGao2006

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm


from datetime import datetime
import random, string
import numpy as np
import nupack
import sys, os


SCRIPT_DIR = 'case2_scatterplots/'
myMultistrand.setNumOfThreads(8)

GLOBAL_TEMPERATURE = 20.0
ATIME_OUT = 1.0
TIME_MEASUREMENT = 0.0;

OFFSET = 0;
RANGE1 = range(1);
RANGE2 = range(1);
RANGE3 = range(1);

MARKER_SIZE = 40

TOTAL_NUM_OF_SEQ = 0;



def create_setup(num_traj, strand_seq):
    
    options = standardOptions("First Step", GLOBAL_TEMPERATURE, num_traj, ATIME_OUT)
    hybridization(options, strand_seq, num_traj)    
    options.uniformRates()
    setSaltGao2006(options)    
    
    return options

def first_step_simulation(strand_seq, num_traj, rate_method_k_or_m="Metropolis", concentration=50e-9):

    print "Running first step mode simulations for %s (with Boltzmann sampling)..." % (strand_seq)
    
    myMultistrand.setOptionsFactory2(create_setup, num_traj, strand_seq)
    myMultistrand.setLeakMode()
    myMultistrand.run()
    
    return myMultistrand.results, int(myMultistrand.nForward.value), int(myMultistrand.nReverse.value)
    

def WC(seq):
    """Computes the Watson-Crick complement for a DNA sequence."""
    return seq.translate(string.maketrans('ACTG', 'TGAC'))[::-1]

def randomseq(length, bases='ACTG'):
    """Chooses a random DNA sequence."""
    return ''.join(random.choice(bases) for i in range(length))

def stemsize(seq, T=GLOBAL_TEMPERATURE, material='dna'):
    result = nupack.mfe([seq], T=T, material=material)
    struct = result[0][0]
    n = len(struct)
    try: 
        stemstart5 = struct.index('(')
        stemlen5 = struct[stemstart5:].index('.')
        revstruct = struct[::-1]
        stemstart3 = revstruct.index(')')
        stemlen3 = revstruct[stemstart3:].index('.')
        stemlen = min(stemlen5, stemlen3)
    except ValueError:
        stemlen = 0
    return stemlen

def toeholds(seq, T=GLOBAL_TEMPERATURE, material="dna"):
    result = nupack.mfe([seq], T=T, material=material)
    struct = result[0][0]
    n = len(struct)
    try: 
        toe5len = struct.index('(')
        toe3len = struct[::-1].index(')')
    except ValueError:
        toe5len = n
        toe3len = n
    # print "    %s = %s , toes %d and %d" % (seq, struct, toe5len, toe3len)
    return (toe5len, toe3len)
    
def binding_dG(seq, T=GLOBAL_TEMPERATURE, material='dna'):
    dG_top = nupack.pfunc([seq], T=T, material=material)
    dG_bot = nupack.pfunc([ WC(seq) ], T=T, material=material)
    dG_duplex = nupack.pfunc([ seq, WC(seq) ], T=T, material=material)
    return (dG_duplex - dG_top - dG_bot)

def duplex_dG(seq, T=GLOBAL_TEMPERATURE, material='dna'):
    return nupack.pfunc([ seq, WC(seq) ], T=T, material=material)


def strand_dG(seq, T=GLOBAL_TEMPERATURE, material='dna'):
    return nupack.pfunc([ seq ], T=T, material=material)


def reverse_rate(seq, kf, T=GLOBAL_TEMPERATURE, material='dna'):
    dG = binding_dG(seq, T=GLOBAL_TEMPERATURE, material='dna')
    RT = 0.001987 * (273.15 + GLOBAL_TEMPERATURE)  # kcal/mol at 25 C
    kr = kf * np.exp(dG / RT)
    return kr

def trim(prior, topN, mean_dG=None):
    
    candidates = list(set(prior))  # get rid of duplicates
    dGs = [ duplex_dG(seq) for seq in candidates ]
    
    if(mean_dG == None):    
        mean_dG = np.mean(dGs)
    dist_dG = [ abs(dG - mean_dG) for dG in dGs ]
    ranked_candidates = sorted(zip(candidates, dist_dG), key=lambda z: z[1])
    candidates = [ s for (s, d) in ranked_candidates[0:topN] ]
    
    return candidates, mean_dG


# this just looks at struct of sense strand.  it would be better to look at both, even if they have different structures. TO DO.
def toebinding(seq, T=GLOBAL_TEMPERATURE, material='dna'):
    
    n = len(seq)
    (toe5len, toe3len) = toeholds(seq, T=T, material=material)   
    
    if toe5len < n:
        toe5seq = seq[:toe5len] if toe5len > 0 else ''
        toe3seq = seq[-toe3len:] if toe3len > 0 else ''
        base = seq[toe5len]
        fake_seq1 = toe5seq + base + 'GGGGTTTTCCCC' + WC(base) + toe3seq
        fake_struct1 = '.' * toe5len + '(((((....)))))' + '.' * toe3len
        fake_seq2 = WC(fake_seq1)
        fake_struct2 = '.' * toe3len + '(((((....)))))' + '.' * toe5len
        fake_struct_toe5 = '(' * toe5len + '(((((....)))))' + '.' * toe3len + '+' + '.' * toe3len + '(((((....)))))' + ')' * toe5len
        fake_struct_toe3 = '.' * toe5len + '(((((....)))))' + '(' * toe3len + '+' + ')' * toe3len + '(((((....)))))' + '.' * toe5len
        fake_struct_toes = '(' * toe5len + '(((((....)))))' + '(' * toe3len + '+' + ')' * toe3len + '(((((....)))))' + ')' * toe5len
        # print "Sequence %s:  %s %s and %s %s" % (seq,fake_seq1,fake_struct1,fake_seq2,fake_struct2)
        dG1 = nupack.energy([fake_seq1], fake_struct1, T=T, material=material)
        dG2 = nupack.energy([fake_seq2], fake_struct2, T=T, material=material)
        if toe5len > 0:
            dG_toe5 = nupack.energy([fake_seq1, fake_seq2], fake_struct_toe5, T=T, material=material) - dG1 - dG2
        else:
            dG_toe5 = 0
        if toe3len > 0:
            dG_toe3 = nupack.energy([fake_seq1, fake_seq2], fake_struct_toe3, T=T, material=material) - dG1 - dG2
        else:
            dG_toe3 = 0
        if toe5len > 0 and toe3len > 0:
            dG_toes = nupack.energy([fake_seq1, fake_seq2], fake_struct_toes, T=T, material=material) - dG1 - dG2
        else:
            dG_toes = 0
        return max(-dG_toe5, -dG_toe3, -dG_toes)
    else:
        # print "Sequence %s: no structure!" % seq
        return -duplex_dG(seq, T=T, material=material)

if __name__ == '__main__':
    
    
    TIME_MEASUREMENT = datetime.now()



    if len(sys.argv) < 2:
        print """Usage:
              python -i hybridization_scatterplot generate random <len> <num seqs> <data file name>
              python -i hybridization_scatterplot generate iso-random <len> <num seqs> <data file name>
              python -i hybridization_scatterplot generate structured <len> <max toe> <max stem> <num seqs> <data file name>
              python -i hybridization_scatterplot generate iso-structured <len> <max toe> <max stem> <num seqs> <data file name>
              python -i hybridization_scatterplot generate fixed-stem <len> <max toe> <stem size> <num seqs> <data file name>
              python -i hybridization_scatterplot generate iso-fixed-stem <len> <max toe> <stem size> <num seqs> <data file name>
              python -i hybridization_scatterplot plot <data file names>
              """
        sys.exit()

    if sys.argv[1] == 'generate':
        if sys.argv[2] in ['random', 'iso-random', 'structured-random'] and len(sys.argv) == 6:
            L = int(sys.argv[3])
            N = int(sys.argv[4])
            filename = sys.argv[5]
        elif sys.argv[2] in ['structured', 'structured-random', 'iso-structured', 'iso-structured-random', 'fixed-stem', 'iso-fixed-stem'] and len(sys.argv) == 8:
            L = int(sys.argv[3])
            toemax = int(sys.argv[4])
            stemmax = int(sys.argv[5])
            N = int(sys.argv[6])
            filename = sys.argv[7]
        else:   
            print """Usage:
              python -i hybridization_scatterplot generate random <len> <num seqs> <data file name>
              python -i hybridization_scatterplot generate iso-random <len> <num seqs> <data file name>
              python -i hybridization_scatterplot generate structured <len> <max toe> <max stem> <num seqs> <data file name>
              python -i hybridization_scatterplot generate iso-structured <len> <max toe> <max stem> <num seqs> <data file name>
              python -i hybridization_scatterplot generate fixed-stem <len> <max toe> <stem size> <num seqs> <data file name>
              python -i hybridization_scatterplot generate iso-fixed-stem <len> <max toe> <stem size> <num seqs> <data file name>
              python -i hybridization_scatterplot plot <data file names>
              """
            sys.exit()
        if not set(filename) <= set('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-'):  # don't allow crazy symbols
            print "Data file postfix names must be alphanumeric only."
            sys.exit()

    if sys.argv[1] == 'plot':
        filenames = []
        for fn in sys.argv[2:]:
            if not set(fn) <= set('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-'):  # don't allow crazy symbols
                print "Data file postfix names must be alphanumeric only."
                sys.exit()
            filenames.append(fn)

    print
    if sys.argv[1] == 'generate' and (sys.argv[2] in ['random', 'iso-random', 'structured-random', 'structured', 'iso-structured', 'iso-structured-random', 'fixed-stem', 'iso-fixed-stem']):   

        
        candidates = []
        
        if(L == 25):  # # 25 mers        

            # # The first 3 sequences will always be P0, P3, P4. 
            candidates.append(goa2006_P0)
            candidates.append(goa2006_P3)
            candidates.append(goa2006_P4)
            
            # adding some very slow sequences. 
            candidates.append('TGAGGGGTTAGACTGAGAACCCCTG')
            candidates.append('ATCCGGTATCTTGGGTTCTCCGGAT')
            candidates.append('GAACCTGTTGTACACGGAGGGTTCT')
            candidates.append('ACCCGCTAACCATAGCGGTTCTTAT')
            candidates.append('ACCATGTGCACTTGCTATCACATGG')
            
            candidates.append('CACACTCGTAGGGGTCTTTGAGTGA')
            candidates.append('CACGCAGTCAATTGACTTCGTGAAC')
            candidates.append('TGCAATCCCTTGTCTCAAGGGATCC')
            candidates.append('TACGCCGTATGTACCTTACGGGTAA')
            candidates.append('AGGGCATAGTTAGGATACGGACCCC')
            
            candidates.append('TCCCGCACATAGATTCAAAGGAGGG')
            candidates.append('TTGCATGATGAGCAGCTCATGCATA')
            candidates.append('TCGATCACTCACTAGGAGTGAACGG')
            candidates.append('AGGCCAATTTTTACCCAGGGCCTAC')
            candidates.append('TTCGGATAAAGTCAACCAGTCCGCT')
            
            candidates.append('AACGTAGAGCAAGGCTAGAAGACGT')
            candidates.append('ATGTTTGTGAGCCCGTCTGAAAACC')
            
            
            OFFSET = 20
            RANGE3 = range(OFFSET)
            
        M = 0
        
        if 'structured' in sys.argv[2] or 'fixed-stem' in sys.argv[2]:   
            M = int(0.5 * N)
        
        
        TOTAL_NUM_OF_SEQ = N + M + OFFSET
        
        RANGE1 = range(OFFSET, N + OFFSET)
        RANGE2 = range(N + OFFSET, TOTAL_NUM_OF_SEQ)  
        
        
        
        randomCandidates = []
        structuredCandidates = []
        
        if 'iso' in sys.argv[2]:
            Ngenerate = 80 * N  # we'll generate many sequences and choose the ones with the closest-to-average duplex energy
            Mgenerate = 80 * M
            print "Generating %d sequences each of length %d, with similar duplex energies..." % (N, L)
        else:
            Ngenerate = N
            Mgenerate = M
            
            
        # first choose the sequences
        if 'random' in sys.argv[2]:
            for seq in [ randomseq(L, 'ACTG') for i in range(Ngenerate) ]:
                randomCandidates.append(seq)
                
        
        if 'structured' in sys.argv[2] or 'fixed-stem' in sys.argv[2]:    
            # generate sequences according to domain-level specifications with a range of hairpin stem and toehold lengths
            # (L,am,bm,cm) = (30,6,8,6)   # up to size 6+2*8+6 = 28 nt in the toeholds and stems.
            # (L,am,bm,cm) = (21,5,6,3)    # up to size 5+2*6+3 = 20 nt in the toeholds and stems.
            (am, bm, cm) = (toemax, stemmax, toemax)
            for i in range(Mgenerate):
                d = 0
                while d <= 0:
                    a = random.randint(3, am) if sys.argv[2] == 'structured' else random.randint(min(2, am), am)
                    b = random.randint(4, bm) if sys.argv[2] == 'structured' else bm
                    c = random.randint(3, cm) if sys.argv[2] == 'structured' else random.randint(min(2, cm), cm)
                    d = L - a - 2 * b - c
                stem = randomseq(b)
                # seq  = randomseq(a,'ATC')+stem+randomseq(d,'T')+WC(stem)+randomseq(c,'ATC')   # tends to have different duplex_dG, which confounds analysis
                seq = randomseq(a, 'ATCG') + stem + randomseq(d, 'ACTG') + WC(stem) + randomseq(c, 'ATCG')
                structuredCandidates.append(seq)
                
         
        if 'iso' in sys.argv[2]:  
            randomCandidates, mean_dG = trim(randomCandidates, N)
            structuredCandidates, mean_dG = trim(structuredCandidates, M, mean_dG)

                        
        for seq in structuredCandidates:
            candidates.append(seq)
        for seq in randomCandidates:
            candidates.append(seq)
        
        
        for seq in candidates:
            print("The sequence is ", seq)

        
        results = []
        timeOut = list(range(len(candidates)))
        
        for i in range(len(candidates)):
            
            timeOut[i] = False
            
#             print
            seq = candidates[i]
            print "%d: Simulating %s DNA strand at 20 C:  %s with toeholds %d and %d, stem %d, ss dG = %g and %g, ds dG = %g " % \
                    (i, sys.argv[2], seq, toeholds(seq)[0], toeholds(seq)[1], stemsize(seq), strand_dG(seq), strand_dG(WC(seq)), duplex_dG(seq))

            # Accumulate statistics from all the runs, until enough succesful simulations are collected.
            trials = 60
            CONCENTRATION = 50e-9
            totalRates = FirstStepLeakRate()

            nForward = 0;
            nReverse = 0;
            
            while ((nForward < 25 or (nReverse == 0)) and not timeOut[i]) :  # this is a bit wasteful, but "only" by a factor of about 1.5.  Aims for 20% error bars on k1.

                print 
                batchRates, nF, nR = first_step_simulation(seq, trials, concentration=CONCENTRATION) 
                

                if(totalRates == None):
                    totalRates = batchRates
                else:
                    totalRates.merge(batchRates,deepCopy=True) 

                print "Inferred rate constants with analytical error bars (all runs so far):"
                print totalRates

                nForward = nForward + nF
                nReverse = nReverse + nR


                if totalRates.nReverse > 4500000:  # abort this by extrapolating the rate when 0.5 million trajectories have been run.                    
                    timeOut[i] = True

                if trials < 50000:  # memory problems running huge numbers of trials.  maybe we need to do better:  compress dataset incrementally...?
                    trials *= 2
                                      
                    
            print "Final statistics are "
            print totalRates

#             results.append([seq, nForward, nReverse, totalRates.kEff(concentration=CONCENTRATION)])
            results.append([seq, nForward, nReverse, totalRates.k1()])
        

    else:
        print "Didn't understand args."
        sys.exit()
        


    def writeResults(inResults, name):
    
        # first make sure /scatterdata is available
        if not os.path.exists(SCRIPT_DIR):
            os.makedirs(SCRIPT_DIR)
    
        # now also save the sequences 
        f = open(SCRIPT_DIR + name + "-details" + ".txt", 'w')
        
        f.write("SIMULATION TIME WAS: " + str(datetime.now() - TIME_MEASUREMENT) + "          \n \n");
        f.write("OFFSET, N, M =  " + str(OFFSET) + "  " + str(N) + "  " + str(M) + "  \n");
        
        f.close()
        
        # now also save the sequences 
        f = open(SCRIPT_DIR + name + "-seqs" + ".txt", 'w')
        
        
        for i in range(len(inResults)):
        
            f.write(str(i))
            f.write(" " + str(inResults[i][0]))
            f.write(" " + str(timeOut[i]))
            f.write(" " + str(inResults[i][1]))
            f.write(" " + str(inResults[i][2]))
            f.write(" " + str(np.log10(inResults[i][3])))
            f.write(" " + str(np.log10(reverse_rate(inResults[i][0], inResults[i][3]))) + "\n")
        
        f.close()

    
    writeResults(results, filename)
    
    
    RANGE1 = [i for i in RANGE1  if not timeOut[i]]
    RANGE2 = [i for i in RANGE2  if not timeOut[i]]
    RANGE3 = [i for i in RANGE3  if not timeOut[i]]    
    
    resultsFiltered = [results[i]  for i in range(TOTAL_NUM_OF_SEQ) if timeOut[i] ]
   

    # some old data files contain only (seq,kf); newer ones contains (seq, Nf, Nr, k1, k1p, k2, k2p); but here we only need kf == k1
    results = [ ((data[0], data[3]) if len(data) == 4 else (data[0], data[1])) for data in results ]  
    resultsFiltered = [ ((data[0], data[3]) if len(data) == 4 else (data[0], data[1])) for data in resultsFiltered  ]  
    
    
#     # NOW : Write to file! 
#     
#     
#     f = open(SCRIPT_DIR + filename + "-storage" + ".txt", 'w')
#     
#     f.write("SIMULATION TIME WAS: " + str(datetime.now() - TIME_MEASUREMENT) + "          \n \n");
#     f.write("OFFSET, N, M =  " + str(OFFSET) + "  " + str(N) + "  " + str(M) + "  \n");
#     
       
    
    
    # now scrap all results that did not finish
    




    print
    print "Calculating stats on database of simulated strands..."
    toe_adj = 2
    toes = [ min(max(max(toeholds(seq)), sum(toeholds(seq)) - toe_adj), 17) for (seq, kf) in results ]  # effective length of toeholds
    stems = [ stemsize(seq) for (seq, kf) in results ]  # length of duplex stem
    toe_dGs = [ toebinding(seq) for (seq, kf) in results ]  # (positive number, 1 or 2) for best energy of binding by 1 or 2 toeholds
    stem_dGs = [ max(-strand_dG(seq), -strand_dG(WC(seq))) for (seq, kf) in results ]  # positive number, whichever is stronger
    binding_dGs = [ binding_dG(seq) for (seq, kf) in results ]
    duplex_dGs = [ duplex_dG(seq) for (seq, kf) in results ]
    strand_dGs = [ strand_dG(seq) + strand_dG(WC(seq)) for (seq, kf) in results ]   
    kfs = [ kf for (seq, kf) in results ]
    krs = [ reverse_rate(seq, kf) for (seq, kf) in results ]
    log_kfs = [ np.log10(kf) for (seq, kf) in results]  # forward rates
    log_krs = [ np.log10(kr) for kr in krs]  # reverse rates


    extremes = ''
    i = np.argmax(kfs)
    extremes += "\nFastest association: %s at %g /M/s" % (results[i][0], kfs[i])
    i = np.argmin(kfs)
    extremes += "\nSlowest association: %s at %g /M/s" % (results[i][0], kfs[i])
    i = np.argmax(krs)
    extremes += "\nFastest dissociation: %s at %g /s" % (results[i][0], krs[i])
    i = np.argmin(krs)
    extremes += "\nSlowest dissociation: %s at %g /s" % (results[i][0], krs[i])


    with PdfPages(SCRIPT_DIR + filename + ".pdf") as pdf:
        print "Drawing plots, saving to '" + SCRIPT_DIR + "'" + filename + ".pdf'..." 

        # Two subplots, the axes array is 1-d
        plt.figure(1)
        plt.subplots_adjust(hspace=1.0)
        plt.subplot(211)
        n, bins, patches = plt.hist(log_kfs, 20, normed=1, facecolor='green', alpha=0.75)
        plt.title("Association and dissociation rate distributions")
        plt.ylabel("frequency of\nassociation rates", fontsize='large')
        plt.yticks(fontsize='larger', va='bottom')
        plt.xlabel("Log10 rate constant kf (/M/s)" + extremes, fontsize='large')
        plt.xticks(fontsize='larger')
        plt.subplot(212)
        n, bins, patches = plt.hist(log_krs, 20, normed=1, facecolor='green', alpha=0.75)
        plt.ylabel("frequency of\ndissociation rates", fontsize='large')
        plt.yticks(fontsize='larger', va='bottom')
        plt.xlabel("Log10 rate constant kr (/s)", fontsize='large')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        # Do the rates depend upon how much secondary structure the single strands have?
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(strand_dGs, log_kfs, s=[10 * s for s in stems], alpha=0.5, cmap=cm.jet)
        plt.title("Association and dissociation rates")
        plt.ylabel("Log10 rate constant kf (/M/s) \n", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.subplot(212)
        plt.scatter(strand_dGs, log_krs, s=[10 * s for s in stems], alpha=0.5, cmap=cm.jet)
        plt.ylabel("Log10 rate constant kr (/s)", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.xlabel("single-strand secondary structure (kcal/mol)", fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        # Do the rates depend upon the overall dG of the reaction
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(binding_dGs, log_kfs, s=[10 * s for s in stems], alpha=0.5, cmap=cm.jet)
        plt.title("Association and dissociation rates")
        plt.ylabel("Log10 rate constant kf (/M/s) \n", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.subplot(212)
        plt.scatter(binding_dGs, log_krs, s=[10 * s for s in stems], alpha=0.5, cmap=cm.jet)
        plt.ylabel("Log10 rate constant kr (/s)", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.xlabel("net binding strength (kcal/mol)", fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        # Do the rates depend upon the strength of the duplex state?
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(duplex_dGs, log_kfs, s=[10 * s for s in stems], alpha=0.5, cmap=cm.jet)
        plt.title("Association and dissociation rates")
        plt.ylabel("Log10 rate constant kf (/M/s) \n", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.subplot(212)
        plt.scatter(duplex_dGs, log_krs, s=[10 * s for s in stems], alpha=0.5, cmap=cm.jet)
        plt.ylabel("Log10 rate constant kr (/s)", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.xlabel("duplex strength (kcal/mol)", fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        singletoes = [ max(toeholds(seq)) for (seq, kf) in results ]
        toecolors = ['red' if s == t else 'blue' for (s, t) in zip(singletoes, toes)]

        # Do the rates depend upon toehold lengths?
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(toes, log_kfs, s=[10 * s for s in stems], color=toecolors, alpha=0.5, cmap=cm.jet)
        plt.title("Association and dissociation rates (double-toeholds are blue)")
        plt.ylabel("Log10 rate constant kf (/M/s) \n", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.subplot(212)
        plt.scatter(toes, log_krs, s=[10 * s for s in stems], color=toecolors, alpha=0.5, cmap=cm.jet)
        plt.ylabel("Log10 rate constant kr (/s)", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.xlabel("effective toehold length (nt)", fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        toes = np.array(toes)
        stems = np.array(stems)
        toe_dGs = np.array(toe_dGs)  # max of the MFE difference for toehold binding configurations (i.e. positive number)
        stem_dGs = np.array(stem_dGs)  # max of the two (negated) ss pfuncs (i.e. positive number)
        binding_dGs = np.array(binding_dGs)
        duplex_dGs = np.array(duplex_dGs)
        strand_dGs = np.array(strand_dGs)  # sum of both ss pfuncs  (negative number)
        kfs = np.array(kfs)
        krs = np.array(krs)
        log_kfs = np.array(log_kfs)
        log_krs = np.array(log_krs)
        long = np.array([ s > 3 for s in stems ])  # numpy, when indexed by Boolean array, selects that subset of items
        short = np.array([ s < 3 for s in stems ])

        # Do the rates depend upon toehold lengths, for long-stemmed strands?
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(toes[long], log_kfs[long], s=[10 * s for s in stems[long]], color=toecolors, alpha=0.5, cmap=cm.jet)
        plt.title("Association and dissociation rates for long-stem strands")
        plt.ylabel("Log10 rate constant kf (/M/s) \n", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.subplot(212)
        plt.scatter(toes[long], log_krs[long], s=[10 * s for s in stems[long]], color=toecolors, alpha=0.5, cmap=cm.jet)
        plt.ylabel("Log10 rate constant kr (/s)", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.xlabel("effective toehold length (nt)", fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        # Do the rates depend upon toehold lengths, for short-stemmed strands?
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(toes[short], log_kfs[short], s=[10 * s for s in stems[short]], color=toecolors, alpha=0.5, cmap=cm.jet)
        plt.title("Association and dissociation rates for short-stem strands")
        plt.ylabel("Log10 rate constant kf (/M/s) \n", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.subplot(212)
        plt.scatter(toes[short], log_krs[short], s=[10 * s for s in stems[short]], color=toecolors, alpha=0.5, cmap=cm.jet)
        plt.ylabel("Log10 rate constant kr (/s)", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.xlabel("effective toehold length (nt)", fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        # Do the rates depend upon toehold strengths?
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(toe_dGs, log_kfs, s=[10 * s for s in stems], color=toecolors, alpha=0.5, cmap=cm.jet)
        plt.title("Association and dissociation rates (double-toeholds are blue)")
        plt.ylabel("Log10 rate constant kf (/M/s) \n", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.subplot(212)
        plt.scatter(toe_dGs, log_krs, s=[10 * s for s in stems], color=toecolors, alpha=0.5, cmap=cm.jet)
        plt.ylabel("Log10 rate constant kr (/s)", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.xlabel("toehold binding strength (kcal/mol)", fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()


        r1 = np.array([ (random.random() - 0.5) / 2.0 for i in range(len(toes)) ])  # jitter to avoid super-positioning
        r2 = np.array([ (random.random() - 0.5) / 2.0 for i in range(len(toes)) ])  # jitter to avoid super-positioning

        # Not quite a contour plot, but can we illustrate kf as a function of toehold length and stem length
        plt.figure(1)
        plt.subplot(111)
        plt.scatter(toes + r1, stems + r2, s=50, c=log_kfs, alpha=0.5, cmap=cm.jet)
        # plt.gray()
        cb = plt.colorbar()
        cb.solids.set_edgecolor("face")
        cb.set_label('log10 of kf')
        plt.title("Association rates in toehold-length vs stem-length space\n")
        plt.ylabel("stem length", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.xlabel("effective toehold length", fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()

        plt.figure(1)
        plt.subplot(111)
        plt.scatter(toes + r1, stems + r2, s=50, c=log_krs, alpha=0.5, cmap=cm.jet)
        # plt.gray()
        cb = plt.colorbar()
        cb.set_label('log10 of kr')
        cb.solids.set_edgecolor("face")
        plt.title("Disassociation rates in toehold-length vs stem-length space\n")
        plt.ylabel("stem length", fontsize='larger')
        plt.yticks(fontsize='larger', va='bottom')
        plt.xlabel("effective toehold length", fontsize='larger')
        plt.xticks(fontsize='larger')
        pdf.savefig()
        plt.close()
        
        
        def filter(array, range):
        
            return array[timeOut[range] == False]
        
        def genericScatter(toe, stem, rates):
            
            plt.figure(1)
#             ax = plt.subplot(111)  
            
            plt.scatter(toe_dGs[RANGE1], stem_dGs[RANGE1], s=MARKER_SIZE, c=log_kfs[RANGE1], alpha=0.5, cmap=cm.jet)

            if(len(toe_dGs) > (N + OFFSET)):
                plt.scatter(toe_dGs[RANGE2], stem_dGs[RANGE2], s=MARKER_SIZE, c=log_kfs[RANGE2], alpha=0.5, marker='D', cmap=cm.jet)
                          
            if L == 25:               
                plt.scatter(toe_dGs[RANGE3], stem_dGs[RANGE3], s= (2* MARKER_SIZE), c=log_kfs[RANGE3], alpha=1.0, marker='*', cmap=cm.jet)

            cb = plt.colorbar()
            cb.solids.set_edgecolor("face")
            

        def strengthPlotA(toe_dGs, stem_dGs, log_kfs):

            # Not quite a contour plot, but can we illustrate kf as a function of toehold strength and stem strength
            
            genericScatter(toe_dGs, stem_dGs, log_kfs)
                
            cb.set_label('log10 of kf')
            cb.solids.set_edgecolor("face")
            plt.title("Association rates in toehold-strength vs stem-strength space\n")
            plt.ylabel("stem strength", fontsize='larger')
            plt.yticks(fontsize='larger', va='bottom')
            plt.xlabel("toehold strength", fontsize='larger')
            plt.xticks(fontsize='larger')
        
        
        def strengthPlotB(toe_dGs, stem_dGs, log_krs):         
            
            genericScatter(toe_dGs, stem_dGs, log_krs)

            cb.set_label('log10 of kr')
            cb.solids.set_edgecolor("face")
            plt.title("Disassociation rates in toehold-strength vs stem-strength space\n")
            plt.ylabel("stem strength", fontsize='larger')
            plt.yticks(fontsize='larger', va='bottom')
            plt.xlabel("toehold strength", fontsize='larger')
            plt.xticks(fontsize='larger')
            
        def annotateThat(toe_dGs, stem_dGs):
            for label, x, y in zip(range(len(toe_dGs)), toe_dGs, stem_dGs):
                plt.annotate(label, xy=[x, y], size=2)
        
        strengthPlotA(toe_dGs, stem_dGs, log_kfs)
        pdf.savefig()
        plt.close()

                
        strengthPlotA(toe_dGs, stem_dGs, log_kfs)
        annotateThat(toe_dGs, stem_dGs)
        pdf.savefig()
        plt.close()

         
        strengthPlotB(toe_dGs, stem_dGs, log_krs)   
        pdf.savefig()
        plt.close()


        strengthPlotB(toe_dGs, stem_dGs, log_krs)   
        annotateThat(toe_dGs, stem_dGs)
        pdf.savefig()
        plt.close()



        number_distinct_lengths = len(set([len(seq) for (seq, kf) in results]))
        
        duplex_DG_std = np.std(duplex_dGs)
        
        if(L == 25):
            duplex_DG_std = np.std(duplex_dGs[3:(N + 3)])
        
        print "duplex dG standard deviation = %g" % duplex_DG_std
#         if number_distinct_lengths == 1 and duplex_DG_std < .2:
        if True:
            log_fastest_assoc = np.max(log_kfs)
            log_slowest_dissoc = np.min(log_krs)
            dissoc_speed_ups = log_krs - log_slowest_dissoc
            assoc_slow_downs = log_fastest_assoc - log_kfs
            updown_range = round(max(max(dissoc_speed_ups), max(assoc_slow_downs)))


            def finalPlots(dissoc_speed_ups, assoc_slow_downs, colorsIn):
                                
                plt.figure(1)
                plt.subplot(111)
                
#                 iconSize = 20
                
                plt.scatter(dissoc_speed_ups[RANGE1], assoc_slow_downs[RANGE1], c=colorsIn[RANGE1], s=MARKER_SIZE, alpha=0.5, cmap=cm.jet)

                if(len(dissoc_speed_ups)) > (N + OFFSET):
                    plt.scatter(dissoc_speed_ups[RANGE2], assoc_slow_downs[RANGE2], c=colorsIn[RANGE2], s=MARKER_SIZE, alpha=0.5, marker='D', cmap=cm.jet)
                    
                                        
                plt.plot([0, updown_range], [0, updown_range], 'r--')
                plt.title("Secondary structure slows down association \n or speeds up dissociation?")
                plt.ylabel("Log10 slow-down of k+ (relative to fastest)", fontsize='larger')
                plt.yticks(fontsize='larger', va='bottom')
                plt.xlabel("Log10 speed-up of k- (relative to slowest)", fontsize='larger')
                plt.xticks(fontsize='larger')
                
                if L == 25:
                    plt.scatter(dissoc_speed_ups[RANGE3], assoc_slow_downs[RANGE3], c=colorsIn[RANGE3], s=(2 * MARKER_SIZE), alpha=1.0, marker='*', cmap=cm.jet)


            # Does secondary structure primarily speed up dissociation rather than slow down association?
                
            blueColors = np.array(['b' for x in range(TOTAL_NUM_OF_SEQ)])
            finalPlots(dissoc_speed_ups, assoc_slow_downs, blueColors)
            pdf.savefig()
            plt.close()

            finalPlots(dissoc_speed_ups, assoc_slow_downs, stem_dGs)
            cb = plt.colorbar()
            cb.set_label('stem strength (kcal/mol)')
            cb.solids.set_edgecolor("face")
            pdf.savefig()
            plt.close()
            
            finalPlots(dissoc_speed_ups, assoc_slow_downs, toe_dGs)
            cb = plt.colorbar()
            cb.solids.set_edgecolor("face")
            pdf.savefig()
            plt.close()


            finalPlots(dissoc_speed_ups, assoc_slow_downs, stem_dGs - toe_dGs)           
            cb = plt.colorbar() 
            cb.set_label('stem strength - toehold strength (kcal/mol)')
            cb.solids.set_edgecolor("face")
            pdf.savefig()
            plt.close()
            
            
            finalPlots(dissoc_speed_ups, assoc_slow_downs, stem_dGs - toe_dGs)           
            cb = plt.colorbar() 
            cb.set_label('stem strength - toehold strength (kcal/mol)')
            cb.solids.set_edgecolor("face")
            annotateThat(dissoc_speed_ups, assoc_slow_downs)
            pdf.savefig()
            plt.close()
            
