import sys
import csv

##  q1 + b1 ] ------[q2     q2+b2] ----- [q3 

def main(bed12, VASTDB_PSI):

	exon_introns = dict()
	
	for row in csv.reader(open(bed12), delimiter = '\t'):


		qstarts = list(map(int, row[11].strip(",").split(",")))
		blocksizes = list(map(int, row[10].strip(",").split(",")))

		start = int(row[1])
		strand = row[5]
		bn = int(row[9])
		chrom = row[0]
		qstart = 0

		for q1, q2, q3, b1, b2 in zip(qstarts, qstarts[1:], qstarts[2:], blocksizes, blocksizes[1:]):

			istart_up = start + q1 + b1
			iend_up = start + q2

			estart = start + q2
			eend = start + q2 + b2
			
			istart_down = start + q2 + b2
			iend_down = start + q3
			
			exon = (chrom , estart, eend)
			intron_up = (chrom,  istart_up, iend_up)
			intron_down = (chrom,  istart_down, iend_down)
			
			exon_introns[exon] = (intron_up, intron_down)
			
	for row in csv.reader(open(VASTDB_PSI), delimiter = '\t'):
		
		exon_coords = row[1]
		chrom = exon_coords.split(":")[0]
		estart, eend =  exon_coords.split(":")[1].split("-")
		
		try:
		
			estart = int(estart)
			eend = int(eend)
			
			exon = (chrom , estart, eend)

			if exon in exon_introns:
				intron_up, intron_down = exon_introns[exon]
				intron_up_len = intron_up[-1] - intron_up[-2]
				intron_down_len = intron_down[-1] - intron_down[-2]

				print( row, intron_up, intron_down, str(intron_up_len), str(intron_down_len), sep="\t" )
		
		except ValueError:
			pass
		

if __name__ == '__main__':
	main (sys.argv[1], sys.argv[2])
