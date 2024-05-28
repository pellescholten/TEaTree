import subprocess
import sys

from collections import defaultdict
nested_dict = lambda: defaultdict(nested_dict)

def removefrags(gffp,threshold,outfile):

	gff = str(gffp)

	# print track
	stdout = sys.stdout
	sys.stdout = open(outfile, 'w')

	# Check number of row of header
	print("##gff-version 3")
	print("##repeatcraft")

	with open(gff, "r") as f:
		for line in f:
			if line.startswith("#"):
				continue
			col = line.rstrip().split("\t")

		if int(col[4]) - int(col[3]) >= threshold:
			print(*col,sep="\t")

	sys.stdout.close()
	sys.stderr.write("\n")
	sys.stdout = stdout