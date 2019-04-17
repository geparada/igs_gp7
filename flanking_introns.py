import sys
import csv


def main(bed12, ME_len):

	n = 100

	transcript_intron_info = defaultdict(list)

	min_intron_lenght = 80

	for row in csv.reader(open(bed12), delimiter = '\t'):

		try:


			qName = row[3]
			seq = Transcriptome[qName]

			qstarts = map (int, row[11].strip(",").split(","))
			blocksizes = map(int, row[10].strip(",").split(","))

			start = int(row[1])
			strand = row[5]
			bn = int(row[9])
			chr = row[0]
			qstart = 0

			for q1, q2, b, b2 in zip(qstarts, qstarts[1:], blocksizes, blocksizes[1:]):
			
      
      

if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	Transcriptometabulator(sys.argv[2])
	main (sys.argv[3], int(sys.argv[4]))
