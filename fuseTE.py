import subprocess
import sys

from collections import defaultdict
nested_dict = lambda: defaultdict(nested_dict)


def update_dictionary(d, col, cattrD, grouped):
	d[col[0]][cattrD["ID"]]["lastcol"] = col
	d[col[0]][cattrD["ID"]]["lastend"] = int(col[4])
	d[col[0]][cattrD["ID"]]["Tstart"] = cattrD["Tstart"]
	d[col[0]][cattrD["ID"]]["Tend"] = cattrD["Tend"]
	d[col[0]][cattrD["ID"]]["grouped"] = grouped
	return d

def truefusete(gffp,gapsize,outfile, mergemode, allowed_consensus_overlap):

	gff = str(gffp)

	sort = str(outfile) + "to_be_sorted"

	# print track
	stdout = sys.stdout
	sys.stdout = open(sort, 'w')

	# quick check number of line of the file
	sh = subprocess.run(['wc', '-l',gff], stdout=subprocess.PIPE)
	totalline = int(sh.stdout.split()[0])

	d = nested_dict()
	lastchrom = "" 

	#check  allowed consensus overlap
	if allowed_consensus_overlap.endswith("p"):
		consensus_proportion = float(allowed_consensus_overlap.rstrip('pP')) / 100
	else:
		allowed_consensus_overlap = int(allowed_consensus_overlap)
		consensus_proportion = 0

	# Progress
	dcnt = 0
	lcnt = 0
	sys.stderr.write("\n\nLabelling fragmented elements to be merged (for alignment)...\n\n")


	print("##gff-version 3")
	print("##TEatree labelled")
	print("##These annotations are not merged yet ")

	with open(gff, "r") as f:
		for line in f:
			if line.startswith("#"):
				continue

			# Progress
			dcnt += 1
			if dcnt % 10000 == 0 or (dcnt == totalline and totalline > 10000):
				sys.stderr.write("\rProgress: " + str(dcnt) + "/"+ str(totalline)+ "\n")

			col = line.rstrip().split("\t")

			# Extract attribute
			cattrD = {}
			cattr = col[8].split(";")
			for i in cattr:
				k, v = i.split("=")
				cattrD[k] = v
			cattrD["Tstart"] = int(cattrD["Tstart"])
			cattrD["Tend"] = int(cattrD["Tend"])

			# if changing to the last column, need to print the what have not printed  (lastcol for all families in last chrom)
			if  col[0] != lastchrom:
				if lastchrom != "":
					for family in d[lastchrom]:
						print(*d[lastchrom][family]["lastcol"],sep="\t")

			lastchrom = col[0] # Update lastchrom

			if not d[col[0]][cattrD["ID"]]:
				# first of family on this chrom, cannot merge, just add to dictionary
				d = update_dictionary(d, col, cattrD, False)
				continue

			# check 
			if mergemode == "threshold":
				match = (int(col[3]) - d[col[0]][cattrD["ID"]]["lastend"]) <= gapsize
			elif mergemode == "ID":
				match = d[col[0]][cattrD["ID"]]["lastcol"][9] == col[9]
			else:
				match = d[col[0]][cattrD["ID"]]["lastcol"][9] == col[9] or (int(col[3]) - d[col[0]][cattrD["ID"]]["lastend"]) <= gapsize
	
			if not match:
				# no need to group the two records
				# print the last record of the last family group without adding new label, and update dictionary
				col2print = d[col[0]][cattrD["ID"]]["lastcol"]
				print(*col2print, sep = "\t")
				d = update_dictionary(d, col, cattrD, False)
				continue

			# Is the previous annotation already grouped?
			if not d[col[0]][cattrD["ID"]]["grouped"]: 
				# the lastcol is the first element is this group, just need to check if last and current copies overlap
				
				#check whether strandedness is readable
				if col[6] != "-" and col[6] != "C" and col[6] != "+":
					sys.stderr.write("\r Error, unknown strandedness while labeling annotations to be merged (\'+\', \'-\' or \'C\' expected)\n\r"+col[6]+"\n")
					sys.exit(1)

				# check whether strandness is the same of the annotations a.k.a. Is the orientation of the annotation the same
				if d[col[0]][cattrD["ID"]]["lastcol"][6] != col[6]:
					print(*d[col[0]][cattrD["ID"]]["lastcol"], sep="\t")
					d = update_dictionary(d, col, cattrD, False)
					continue
				
				#calculate annotation on consensus overlap
				consensus_overlap = min(int(cattrD["Tend"]), int(d[col[0]][cattrD["ID"]]["Tend"])) - max(int(cattrD["Tstart"]), int(d[col[0]][cattrD["ID"]]["Tstart"]) + 1)
				
				#if consensus overlap allowed is specified as a proportion of the consensus familyt length, calculate this proportion
				if consensus_proportion:
					allowed_consensus_overlap = int(col[8].split(";")[2].split("=")[1]) * consensus_proportion

				#if consensus family overlap is more than is allowed, default = 25, then do not merge
				if consensus_overlap >= allowed_consensus_overlap:
					print(*d[col[0]][cattrD["ID"]]["lastcol"], sep="\t")
					d = update_dictionary(d, col, cattrD, False)
					continue

				# passed all tests
				lcnt += 2
				# label for merger, first merger of group
				# can open a new group now and update the attr of the last and current copies

				# Make a new label for lastcol and current col
				d[col[0]][cattrD["ID"]]["groupnumber"] = d[col[0]][cattrD["ID"]].get("groupnumber", 0) + 1

				# Mark down the consensus coverage
				# fix this, inefficient, do ranges instead of individual positions
				groupnumber = d[col[0]][cattrD["ID"]]["groupnumber"]

				#write down consensus coverage
				#make a list of ranges that, one range for each of the two grouped annotations
				d[col[0]][cattrD["ID"]][groupnumber]["consensus_coverage"] = [(d[col[0]][cattrD["ID"]]["Tstart"],d[col[0]][cattrD["ID"]]["Tend"]), (cattrD["Tstart"],cattrD["Tend"])]

				# Print lastcol and update dictionary
				lastcol2print = d[col[0]][cattrD["ID"]]["lastcol"]
				attr = ";TEgroup=" + col[0] + "|" + cattrD["ID"] + "|" + str(d[col[0]][cattrD["ID"]]["groupnumber"])
				lastcol2print[8] = lastcol2print[8] + attr
				print(*lastcol2print,sep="\t")
				col[8] = col[8] + attr # Update lastcol

				d = update_dictionary(d, col, cattrD, True)

			else: # previous annotation of this family is already grouped
				# check consensus information (all in last group)
				groupnumber = d[col[0]][cattrD["ID"]]["groupnumber"]
				

				#FIX
				#do not check all positions, inefficient
				#check list of ranges instead, similar to check later
				#also, idem to later, do I need to be that strict?
				#do i also need to check strandedness etc here? i think so?

				#check whether strandedness is readable
				if col[6] != "-" and col[6] != "C" and col[6] != "+":
					sys.stderr.write("\r Error, unknown strandedness while labeling annotations to be merged (\'+\', \'-\' or \'C\' expected)\n\r"+col[6]+"\n")
					sys.exit(1)

				# check whether strandness is the same of the annotations a.k.a. Is the orientation of the annotation the same
				if str(d[col[0]][cattrD["ID"]]["lastcol"][6]) != str(col[6]):
					print(*d[col[0]][cattrD["ID"]]["lastcol"], sep="\t")
					d = update_dictionary(d, col, cattrD, False)
					continue

				#if consensus overlap allowed is specified as a proportion of the consensus familyt length, calculate this proportion
				if consensus_proportion:
					allowed_consensus_overlap = int(col[8].split(";")[2].split("=")[1]) * consensus_proportion

				#check whether annotation on consensus overlaps 
				consensus_overlap = False
				for r in d[col[0]][cattrD["ID"]][groupnumber]["consensus_coverage"]:
					c_start, c_end = r
					consensus_overlap = min(int(cattrD["Tend"]), c_end) - max(int(cattrD["Tstart"]), c_start + 1)

					if consensus_overlap >= allowed_consensus_overlap:
						consensus_overlap = True	
						break

				if consensus_overlap: 
					print(*d[col[0]][cattrD["ID"]]["lastcol"], sep="\t")
					d = update_dictionary(d, col, cattrD, False)
					continue

				# passed all tests, and group already present
				lcnt += 1 
				# add annotation to existing group
				print(*d[col[0]][cattrD["ID"]]["lastcol"], sep="\t")
				
				# Update consensus coverage of grsys;oup
				d[col[0]][cattrD["ID"]][groupnumber]["consensus_coverage"].append((cattrD["Tstart"],cattrD["Tend"]))
				
				# Update last col using label from last label
				attr = ";TEgroup=" + col[0] + "|" + cattrD["ID"] + "|" + str(d[col[0]][cattrD["ID"]]["groupnumber"])
				col[8] = col[8] + attr
				d = update_dictionary(d, col, cattrD, True)
				

		# print the last record for all families from the last chrom
		for family in d[lastchrom]:
			print(*d[lastchrom][family]["lastcol"],sep="\t")

	sys.stdout.close()

	#resort the outfile based on position in shell
	c = "sort -k1,1 -k4,4n -k5,5n " + sort + " >" + outfile + "&& rm " + sort
	subprocess.run(c, shell=True)
	sys.stderr.write(str(lcnt)+ " fragments labelled for merger\n\n")
	sys.stdout = stdout

	return lcnt