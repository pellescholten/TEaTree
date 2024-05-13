import sys
import re

def freqalign(infile, labelfile, outfileclass, outfilefamily):

	# print track
	stdout = sys.stdout


	# Read rlabel
	# Stat variables

	famsin = {}
	classificationin = {}

	with open(infile,"r") as fi:
		for i in range(3):
			next(fi)
		for line in fi:
			col = line.rstrip().split()

			if col[10] == "Simple_repeat" or col[10] == "Low_complexity":
				continue
			if famsin.get(col[9]):
				famsin[col[9]] += 1
			else:
				famsin[col[9]] = 1

			if classificationin.get(col[10]):
				classificationin[col[10]] += 1
			else:
				classificationin[col[10]] = 1


	famslabel = {}
	classificationlabel = {}

	with open(labelfile,"r") as fm:
		for line in fm:
			if line.startswith("#"):
				continue
			col = line.rstrip().split("\t")

			fam = col[8].split(";")[2].split("=")[1]

			if famslabel.get(fam):
				famslabel[fam] += 1
			else:
				famslabel[fam] = 1

			if classificationlabel.get(col[2]):
				classificationlabel[col[2]] += 1
			else:
				classificationlabel[col[2]] = 1



	for k in list(classificationin.keys()):
		if not classificationlabel.get(k):
			classificationlabel[k] = 0

	for k in list(famsin.keys()):
		if not famslabel.get(k):
			famslabel[k] = 0


	classificationin["total"] = sum(classificationin.values())
	classificationlabel["total"] = sum(classificationlabel.values())

	sys.stdout = open(outfileclass, 'w')
	print("#Number of repeats by class before and after overlap resolving and filtering")
	print("#=============================================================")
	print(*["#repeat class","no. before overlap resolving","no. after"], sep="\t")
	for c in list(classificationin.keys()):
		print(*[c,classificationin[c],classificationlabel[c]],sep="\t")

	sys.stdout.close()

	sys.stdout = open(outfilefamily, 'w')
	print("#Number of repeats by family before and after overlap resolving and filtering")
	print("#=============================================================")
	print(*["#repeat family","no. before overlap resolving","no. before merge"], sep="\t")
	for c in list(famsin.keys()):
		print(*[c,famsin[c],famslabel[c]],sep="\t")

	sys.stdout.close()
	sys.stdout = stdout

def freqcontent(infile, bedfile, outfileclass, outfilefamily):

	# print track
	stdout = sys.stdout

	# Read rlabel
	# Stat variables

	famsin = {}
	classificationin = {}

	with open(infile,"r") as fi:
		for i in range(3):
			next(fi)
		for line in fi:
			col = line.rstrip().split()
			if famsin.get(col[9]):
				famsin[col[9]] += 1
			else:
				famsin[col[9]] = 1

			if classificationin.get(col[10]):
				classificationin[col[10]] += 1
			else:
				classificationin[col[10]] = 1

	famsbed = {}
	classificationbed = {}

	with open(bedfile,"r") as fm:
		for line in fm:
			if line.startswith("#"):
				continue
			col = line.rstrip().split("\t")

			info = col[3]
			info_fam = info.split(".")[1]
			info_class = info.split(".")[2]
			if famsbed.get(info_fam):
				famsbed[info_fam] += 1
			else:
				famsbed[info_fam] = 1

			if classificationbed.get(info_class):
				classificationbed[info_class] += 1
			else:
				classificationbed[info_class] = 1


	for k in list(classificationin.keys()):
		if not classificationbed.get(k):
			classificationbed[k] = 0

	for k in list(famsin.keys()):
		if not famsbed.get(k):
			famsbed[k] = 0

	classificationin["total"] = sum(classificationin.values())
	classificationbed["total"] = sum(classificationbed.values())

	sys.stdout = open(outfileclass, 'w')
	print("#Number of repeats by class before overlap and after resolving")
	print("#=============================================================")
	print(*["#repeat class","no. before overlap resolving","no. after overlap resolving"], sep="\t")
	for c in list(classificationin.keys()):
		print(*[c,classificationin[c],classificationbed[c]],sep="\t")
	
	sys.stdout.close()
	sys.stdout = open(outfilefamily, 'w')
	print("#Number of repeats by family before and after overlap resolving")
	print("#=============================================================")
	print(*["#repeat class","no. before overlap resolving","no. after overlap resolving"], sep="\t")
	for c in list(famsin.keys()):
		print(*[c,famsin[c],famsbed[c]],sep="\t")

	sys.stdout.close()
	sys.stdout = stdout

def bpcontent(infile, bedfile, outfile):

	# print track
	stdout = sys.stdout
	sys.stdout = open(outfile, 'w')

	# Read rlabel
	# Stat variables

	#famsin = {}
	classificationin = {}

	with open(infile,"r") as fi:
		for i in range(3):
			next(fi)
		for line in fi:
			col = line.rstrip().split()

			if classificationin.get(col[10]):
				classificationin[col[10]] += int(col[6]) - int(col[5]) + 1
			else:
				classificationin[col[10]] = int(col[6]) - int(col[5]) + 1

	#famsbed = {}
	classificationbed = {}

	with open(bedfile,"r") as fm:
		for line in fm:
			if line.startswith("#"):
				continue
			col = line.rstrip().split("\t")

			info = col[3]
			#info_fam = info.split(".")[1]
			info_class = info.split(".")[2]
			#if famsbed.get(info_fam):
		#		famsbed[info_fam] += 1
			#else:
			#	famsbed[info_fam] = 1

			if classificationbed.get(info_class):
				classificationbed[info_class] += int(col[2]) - int(col[1])
			else:
				classificationbed[info_class] = int(col[2]) - int(col[1])


	for k in list(classificationin.keys()):
		if not classificationbed.get(k):
			breakpoint()
			classificationbed[k] = 0

	classificationin["total"] = sum(classificationin.values())
	classificationbed["total"] = sum(classificationbed.values())


	#for k in list(famsin.keys()):
	#	if not famsbed.get(k):
	#		famsbed[k] = 0

	for c in list(classificationbed.keys()):
		print(*[c,classificationin[c],classificationbed[c]],sep="\t")

	sys.stdout.close()
	sys.stdout = stdout

