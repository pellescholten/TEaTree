import sys
import re
import subprocess
from collections import defaultdict
nested_dict = lambda: defaultdict(nested_dict)

def filterlength(gffp,outfile,lengthfile):


	gff = gffp

	out = str(outfile)

	# print track
	stdout = sys.stdout
	sys.stdout = open(out, 'w')

	sys.stderr.write("\rfiltering alignments based on concensus length (for alignment)...\n")

	lengths = {}
	with open(lengthfile, 'r') as file:
		for line in file:
			parts = line.strip().split('\t')
			if len(parts) == 2:
				fam_info, fam_length = parts
				fam_name = fam_info.split('#')[0]
				fam_threshold = int(int(fam_length) * 0.3)
				lengths[fam_name] = fam_threshold

	
	# Check number of row of header
	print("##gff-version 3")
	print("##repeatcraft")

	with open(gff, "r") as f:
		for line in f:
			if line.startswith("#"):
				continue

			col = line.rstrip().split("\t")
			classification = col[2]
			tempfam = col[8].split(";")[-1].split("=")[1].split("|")
			
			if len(tempfam) > 1:
				fam = tempfam[1]
			else:
				fam = tempfam[0]
			LTR = classification.split('/')[0] == "LTR"
			size = int(col[4]) - int(col[3])

			longLTR = LTR and size > 250

			if fam[0] == "(":
				continue

			minlength = lengths[fam]

			if size > minlength or longLTR:
				print(line, sep="\t",end="")
				
	sys.stdout.close()
	sys.stdout = stdout
