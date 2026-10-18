#!/bin/bash
# Compare TEaTree with RM2Bed.py and RepeatCraft (written for macOS)
# Usage: cd benchmark; bash benchmark.sh dm6   (or hg38)
# RM2Bed.py needs pandas
set -e
G=$1
rm -rf work_$G
mkdir work_$G
cd work_$G

# download the RepeatMasker annotation and the other tools
curl -s https://hgdownload.soe.ucsc.edu/goldenPath/$G/bigZips/$G.fa.out.gz | gunzip > $G.fa.out
curl -s -O https://raw.githubusercontent.com/rmhubley/RepeatMasker/master/util/RM2Bed.py
git clone -q https://github.com/niccw/repeatcraftp.git
git -C repeatcraftp checkout -q eb7a30d

# 1. resolve overlaps (higher score wins)
/usr/bin/time -l python3 ../../TEaTree.py -i $G.fa.out -o teatree -remove_short_fragments 0 > /dev/null 2> teatree.time
/usr/bin/time -l python3 RM2Bed.py -p rm2bed -o higher_score $G.fa.out 2> rm2bed.time

# 2. merge fragments only (same family, < 150 bp apart)
/usr/bin/time -l python3 ../../TEaTree.py -i $G.fa.out -o teatree_merge -mode alignment -mergemode threshold -gapsize 150 -leave_overlap -remove_short_fragments 0 > /dev/null 2> teatree_merge.time

# RepeatCraft needs the RepeatMasker gff, make it from the .out file
awk 'BEGIN{OFS="\t"} NR>3 {
  if ($9 == "C") {strand = "-"; start = $14} else {strand = "+"; start = $12}
  print $5, "RepeatMasker", "similarity", $6, $7, $2, strand, ".", "Target \"Motif:" $10 "\" " start " " $13}' $G.fa.out > rm.gff
cat > repeatcraft.cfg <<CFG
shortTE_size: 100
mapfile: None
ltr_finder_gff: None
max_LTR_size: 10000
LTR_flanking_size: 200
gap_size: 150
CFG
/usr/bin/time -l python3 repeatcraftp/repeatcraft.py -r rm.gff -u $G.fa.out -c repeatcraft.cfg -o repeatcraft 2> repeatcraft.time

# 3. count annotations and overlaps, and get time and memory
python3 ../count.py $G.fa.out teatree.bed rm2bed_rm.bed teatree_merge.merged.gff repeatcraft.rmerge.gff > ../results/$G.txt
for f in *.time; do
  echo "$f: $(awk '/ real/{print $1}' $f) s, $(awk '/maximum resident/{printf "%.2f", $1/1e9}' $f) GB" >> ../results/$G.txt
done
cat ../results/$G.txt
