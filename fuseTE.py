import subprocess
import sys

from collections import defaultdict
nested_dict = lambda: defaultdict(nested_dict)

def truefusete(gffp,gapsize,outfile, mergemode):

	gff = str(gffp)

	sort = str(outfile) + "to_be_sorted"

	# print track
	stdout = sys.stdout
	sys.stdout = open(sort, 'w')

	# quick check number of line of the file
	sh = subprocess.run(['wc', '-l',gff], stdout=subprocess.PIPE)
	totalline = str(int(sh.stdout.split()[0]))

	cnt = 0
	d = nested_dict()
	lastchrom = ""

	# Progress
	dcnt = 0

	sys.stderr.write("\rLabelling fragmented elements (for alignment)...\n")

	print("##gff-version 3")
	print("##repeatcraft")
	with open(gff, "r") as f:
		for line in f:
			if line.startswith("#"):
				continue
			# Progress
			dcnt += 1
			sys.stderr.write("\rProgress:" + str(dcnt) + "/"+ totalline+ "...")

			col = line.rstrip().split("\t")

			# Extract attribute
			cattrD = {}
			cattr = col[8].split(";")
			for i in cattr:
				k, v = i.split("=")
				cattrD[k] = v
			cattrD["Tstart"] = int(cattrD["Tstart"])
			cattrD["Tend"] = int(cattrD["Tend"])

			# if changing to the last column, need to print the what havn't print out (lastcol for all families in last chrom)
			if  col[0] != lastchrom:
				if lastchrom != "":
					#print("new chrom, print remaining lastcol") # debug
					for family in d[lastchrom]:
						print(*d[lastchrom][family]["lastcol"],sep="\t")

			lastchrom = col[0] # Update lastchrom

		

			if d[col[0]][cattrD["ID"]] or d[col[0]][cattrD["ID"]]["lastcol"][9] == col[9]:  # not the first of this family on this chrom
				#print("not first family") # debug

				#if col[0] != d[col[0]][cattrD["ID"]]["lastcol"][0] # new chromosome, maybe unecessary
	
				if mergemode == "threshold":
					match = (int(col[3]) - d[col[0]][cattrD["ID"]]["lastend"]) <= gapsize
				elif mergemode == "ID":
					match = d[col[0]][cattrD["ID"]]["lastcol"][9] == col[9]
				else:
					threshold = (int(col[3]) - d[col[0]][cattrD["ID"]]["lastend"]) <= gapsize
					idmatch = d[col[0]][cattrD["ID"]]["lastcol"][9] == col[9]
					match = threshold or idmatch
				
				if not match: # or col[0] != d[col[0]][cattrD["ID"]]["lastcol"][0]:
					#print("larger than 150") # debug
					# don't need to group the two records
					# print the lastest record of the lastest family group without adding new label
					col2print = d[col[0]][cattrD["ID"]]["lastcol"]
					print(*col2print, sep = "\t")
					# update the dictionary
					d[col[0]][cattrD["ID"]]["lastcol"] = col
					d[col[0]][cattrD["ID"]]["lastend"] = int(col[4])
					d[col[0]][cattrD["ID"]]["Tstart"] = cattrD["Tstart"]
					d[col[0]][cattrD["ID"]]["Tend"] = cattrD["Tend"]
					d[col[0]][cattrD["ID"]]["grouped"] = False
				else:
					#print("less than 150") # debug
					# Is the lastcol carrying a TEgroup label?
					if d[col[0]][cattrD["ID"]]["grouped"]:
						#print("last one have label") # debug
						# check consensus information (all in last group)
						groupnumber = d[col[0]][cattrD["ID"]]["groupcnt"]
						o = False
						for i in range(cattrD["Tstart"],cattrD["Tend"]):
							if i in list(d[col[0]][cattrD["ID"]][groupnumber]):
								o = True
								break
						if o: # overlap with some copies in the group, break the grouping
							#print("consensus overlap") # debug
							print(*d[col[0]][cattrD["ID"]]["lastcol"], sep="\t")
							d[col[0]][cattrD["ID"]]["lastcol"] = col
							d[col[0]][cattrD["ID"]]["lastend"] = int(col[4])
							d[col[0]][cattrD["ID"]]["Tstart"] = cattrD["Tstart"]
							d[col[0]][cattrD["ID"]]["Tend"] = cattrD["Tend"]
							d[col[0]][cattrD["ID"]]["grouped"] = False #????
						else:
							#print("consensus pass") # debug
							# print the lastcol directly
							print(*d[col[0]][cattrD["ID"]]["lastcol"], sep="\t")
							# Update consensus coverage
							for i in range(cattrD["Tstart"],cattrD["Tend"]):
								d[col[0]][cattrD["ID"]][groupnumber][i] = 1

							# Update last col using label from last label
							attr = ";TEgroup=" + col[0] + "|" + cattrD["ID"] + "|" + str(d[col[0]][cattrD["ID"]]["groupcnt"])
							col[8] = col[8] + attr
							d[col[0]][cattrD["ID"]]["lastcol"] = col
							d[col[0]][cattrD["ID"]]["lastend"] = int(col[4])
							d[col[0]][cattrD["ID"]]["Tstart"] = cattrD["Tstart"]
							d[col[0]][cattrD["ID"]]["Tend"] = cattrD["Tend"]
							d[col[0]][cattrD["ID"]]["grouped"] = True


					else: # the lastcol is the first element is this group, just need to check if last and current copies overlap
						#print("last copy no label") # debug

						#if d[col[0]][cattrD["ID"]]["Tstart"] >= cattrD["Tstart"]:
						#	o = cattrD["Tend"] - d[col[0]][cattrD["ID"]]["Tstart"]
						#else:
						#	o = d[col[0]][cattrD["ID"]]["Tend"] - cattrD["Tstart"]
						#o = min(d[col[0]][cattrD["ID"]]["Tend"], cattrD["Tstart"]) - max(d[col[0]][cattrD["ID"]]["Tstart"], cattrD["Tend"]) #old

						strandedness_not_equal = d[col[0]][cattrD["ID"]]["lastcol"][6] != col[6]
						if col[6] == "+":
							if int(cattrD["Tstart"]) > int(d[col[0]][cattrD["ID"]]["Tend"]):
								overlapping_concensus = False
							else:
								overlapping_concensus = True
						elif col[6] == "-":
							if int(cattrD["Tend"]) < int(d[col[0]][cattrD["ID"]]["Tstart"]):
								overlapping_concensus = False
							else:
								overlapping_concensus = True	
						else:
							sys.stderr.write("\r error,unknown strandedness")
							sys.stderr.write("\r"+col[6])
							sys.exit(1)


						if overlapping_concensus or strandedness_not_equal:  # Consensus position overlap, don't need to group them just print
							#print("consensus overlap") # debug
							print(*d[col[0]][cattrD["ID"]]["lastcol"], sep="\t")
							d[col[0]][cattrD["ID"]]["lastcol"] = col
							d[col[0]][cattrD["ID"]]["lastend"] = int(col[4])
							d[col[0]][cattrD["ID"]]["Tstart"] = cattrD["Tstart"]
							d[col[0]][cattrD["ID"]]["Tend"] = cattrD["Tend"]
							d[col[0]][cattrD["ID"]]["grouped"] = False

						else: # can open a new group now and update the attr of the last and current copies
							#print("consensus pass") # debug
							# Make a new label for lastcol and current col
							if d[col[0]][cattrD["ID"]]["groupcnt"]: # Is there a previous family group in the same chrom?
								d[col[0]][cattrD["ID"]]["groupcnt"] = d[col[0]][cattrD["ID"]]["groupcnt"] + 1
							else:
								d[col[0]][cattrD["ID"]]["groupcnt"] = 1

							# Mark down the consensus coverage
							groupnumber = d[col[0]][cattrD["ID"]]["groupcnt"]
							for i in range(d[col[0]][cattrD["ID"]]["Tstart"],d[col[0]][cattrD["ID"]]["Tend"]):
								d[col[0]][cattrD["ID"]][groupnumber][i] = 1
							for i in range(cattrD["Tstart"],cattrD["Tend"]):
								d[col[0]][cattrD["ID"]][groupnumber][i] = 1

							# Print lastcol
							lastcol2print = d[col[0]][cattrD["ID"]]["lastcol"]
							attr = ";TEgroup=" + col[0] + "|" + cattrD["ID"] + "|" + str(d[col[0]][cattrD["ID"]]["groupcnt"])
							lastcol2print[8] = lastcol2print[8] + attr
							print(*lastcol2print,sep="\t")

							# Update lastcol
							col[8] = col[8] + attr
							d[col[0]][cattrD["ID"]]["lastcol"] = col
							d[col[0]][cattrD["ID"]]["lastend"] = int(col[4])
							d[col[0]][cattrD["ID"]]["Tstart"] = cattrD["Tstart"]
							d[col[0]][cattrD["ID"]]["Tend"] = cattrD["Tend"]
							d[col[0]][cattrD["ID"]]["grouped"] = True

			else: # first family on this chrom
				#print("first element on this chrom") # debug
				d[col[0]][cattrD["ID"]]["lastcol"] = col
				d[col[0]][cattrD["ID"]]["lastend"] = int(col[4])
				d[col[0]][cattrD["ID"]]["Tstart"] = cattrD["Tstart"]
				d[col[0]][cattrD["ID"]]["Tend"] = cattrD["Tend"]
				d[col[0]][cattrD["ID"]]["grouped"] = False

		# print the last record for all families from the last chrom
		#print("print remaining lastcol") # debug
		for family in d[lastchrom]:
			print(*d[lastchrom][family]["lastcol"],sep="\t")

	sys.stdout.close()

	c = "sort -k1,1 -k4,4n -k5,5n " + sort + " >" + outfile + "&& rm " + sort
	subprocess.run(c, shell=True)
	sys.stderr.write("\n")
	sys.stdout = stdout