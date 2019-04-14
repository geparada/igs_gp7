import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

Genome = {}
def Genomictabulator(fasta):
	


	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq

	f.close()


def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0	

def PWM_to_dict(file):
	
	reader = csv.DictReader(file, delimiter = '\t')
	
	A_frec = []
	C_frec = []
	G_frec = []
	T_frec = []
	N_freq = []
	
	matrix = {}
	
	for row in reader:
			
		A_frec.append(float(row["A"]))
		C_frec.append(float(row["C"]))
		G_frec.append(float(row["G"]))
		T_frec.append(float(row["T"]))
		N_freq.append(0)
			
	matrix["A"] = A_frec
	matrix["C"] = C_frec
	matrix["G"] = G_frec
	matrix["T"] = T_frec
	matrix["N"] = N_freq
	
	return matrix


def U2Score(ME_chr, ME_strand, ME_start, ME_end, U2_GTAG_3, U2_GTAG_5,  max_U2_score_3, max_U2_score_5):

    ME5 = str(Genome[ME_chr][ME_start-14:ME_start+3]).upper()
    ME3 = str(Genome[ME_chr][ME_end-3:ME_end+10]).upper()


    if ME_strand == "-":

        ME5 = str(Genome[ME_chr][ME_end-3:ME_end+14].reverse_complement()).upper()
        ME3 = str(Genome[ME_chr][ME_start-10:ME_start+3].reverse_complement()).upper()



    U2_ME5 = 0
    U2_ME3 = 0

    i = 0

    for N in ME5:
        U2_ME5 += U2_GTAG_3[N][i]
        i += 1

    i = 0

    for N in ME3:
        U2_ME3 += U2_GTAG_5[N][i]
        i += 1

    U2_ME5 = percent(U2_ME5, max_U2_score_3)
    U2_ME3 = percent(U2_ME3, max_U2_score_5)

    return(U2_ME5, U2_ME3)


  
def main(exon_bed, GT_AG_U2_5, GT_AG_U2_3):
	
	
	U2_GTAG_5 = PWM_to_dict(U2_GTAG_5_file)
	U2_GTAG_3 = PWM_to_dict(U2_GTAG_3_file)	
	
	with open(GT_AG_U2_5) as U2_GTAG_5_file, open(GT_AG_U2_3) as U2_GTAG_3_file:


		U2_GTAG_5 = PWM_to_dict(U2_GTAG_5_file)
		U2_GTAG_3 = PWM_to_dict(U2_GTAG_3_file)

		U2_GTAG_5_max_score = 0
		U2_GTAG_3_max_score = 0

		for index in range(13):
			U2_GTAG_5_max_score += max(U2_GTAG_5['A'][index], U2_GTAG_5['C'][index], U2_GTAG_5['T'][index], U2_GTAG_5['G'][index])

		for index in range(17):
			U2_GTAG_3_max_score += max(U2_GTAG_3['A'][index], U2_GTAG_3['C'][index], U2_GTAG_3['T'][index], U2_GTAG_3['G'][index])

	
	for row in csv.reader(open(exon_bed), delimiter="\t"):
		
		chrom, start, end, name, cero, strand = row
		
		estart, eend = map(int, name.split("_"))
		
		U2_E5, U2_E3 = U2Score(chrom, strand, estart, eend, U2_GTAG_3, U2_GTAG_5, U2_GTAG_3_max_score, U2_GTAG_5_max_score)

		print name, U2_E5, U2_E3

	
	
  
if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2], sys.argv[3], sys.argv[4])
  
  
