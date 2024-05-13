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
	#sys.stdout = open(out, 'w')

	sys.stderr.write("\rMerging fragmented elements (for alignment)...\n")

	lengths = {}
	with open(lengthfile, 'r') as file:
		for line in file:
			parts = line.strip().split('\t')
			if len(parts) == 2:
				fam_info, fam_length = parts
				fam_name = fam_info.split('#')[0].split("=")[1]
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

			order = col[2]
			fam = col[8].split(";")[-1]
			LTR = fam.split('/')[0] == "LTR"
			size = int(col[4]) - int(col[3])

			minlength = lengths[fam]
			breakpoint()
			if size > minlength and not LTR:
				print(line)
				
	sys.stdout.close()
	sys.stdout = stdout
