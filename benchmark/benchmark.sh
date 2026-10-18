#!/usr/bin/env bash
# Benchmark TEaTree against RM2Bed.py and RepeatCraft on a RepeatMasker .out file from UCSC.
# Usage (macOS): bash benchmark.sh dm6 [python_with_pandas]    (or hg38)
set -e
G=$1; PY=${2:-python3}
cd "$(dirname "$0")"; rm -rf work_$G; mkdir work_$G; cd work_$G

# Download the input and the other tools
curl -s https://hgdownload.soe.ucsc.edu/goldenPath/$G/bigZips/$G.fa.out.gz | gunzip > $G.fa.out
curl -s -O https://raw.githubusercontent.com/rmhubley/RepeatMasker/master/util/RM2Bed.py
git clone -q https://github.com/niccw/repeatcraftp.git && git -C repeatcraftp checkout -q eb7a30d

# 1. Overlap resolution (higher score wins)
/usr/bin/time -l python3 ../../TEaTree.py -i $G.fa.out -o teatree -remove_short_fragments 0 > /dev/null 2> teatree.time
/usr/bin/time -l $PY RM2Bed.py -p rm2bed -o higher_score $G.fa.out 2> rm2bed.time

# 2. Defragmentation only (same family, < 150 bp apart, overlaps left in)
/usr/bin/time -l python3 ../../TEaTree.py -i $G.fa.out -o teatree_defrag -mode alignment -mergemode threshold \
    -gapsize 150 -leave_overlap -remove_short_fragments 0 > /dev/null 2> teatree_defrag.time
awk 'BEGIN{OFS="\t"} NR>3 {s=($9=="C")?"-":"+"; ts=($9=="C")?$14:$12;
     print $5,"RepeatMasker","similarity",$6,$7,$2,s,".","Target \"Motif:"$10"\" "ts" "$13}' $G.fa.out > rm.gff
printf 'shortTE_size: 100\nmapfile: None\nltr_finder_gff: None\nmax_LTR_size: 10000\nLTR_flanking_size: 200\ngap_size: 150\n' > rc.cfg
/usr/bin/time -l python3 repeatcraftp/repeatcraft.py -r rm.gff -u $G.fa.out -c rc.cfg -o repeatcraft 2> repeatcraft.time

# 3. Results
{ python3 ../count.py $G.fa.out teatree.bed rm2bed_rm.bed teatree_defrag.merged.gff repeatcraft.rmerge.gff
  echo; echo -e "run\twall_s\tpeak_RAM_GB"
  for f in *.time; do awk -v f=${f%.time} '/ real/{w=$1} /maximum resident/{r=$1/1e9} END{printf "%s\t%s\t%.2f\n", f, w, r}' $f; done
} | tee ../results/$G.txt
