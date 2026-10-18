#!/usr/bin/env bash
# Compare TEaTree with RM2Bed.py and RepeatCraft on a RepeatMasker .out file from UCSC.
# Usage: bash benchmark.sh dm6 [python_with_pandas]     (or hg38; RM2Bed.py needs pandas)
# Written for macOS (/usr/bin/time -l). Results go to results/<assembly>.tsv
# All tools get the same input file, run one after the other on one thread, and no tool filters on length.
set -e
G=$1; PY=${2:-python3}
HERE=$(cd "$(dirname "$0")" && pwd)
TEATREE="python3 $HERE/../TEaTree.py"
R="$HERE/results/$G.tsv"; mkdir -p "$HERE/results"; rm -f "$R"
mkdir -p "$HERE/work_$G" && cd "$HERE/work_$G"
count() { python3 "$HERE/compare.py" "$R" "$@"; }

# 1. Download the input and the other tools
[ -f $G.fa.out ] || curl -s https://hgdownload.soe.ucsc.edu/goldenPath/$G/bigZips/$G.fa.out.gz | gunzip > $G.fa.out
[ -f RM2Bed.py ] || curl -s -O https://raw.githubusercontent.com/rmhubley/RepeatMasker/master/util/RM2Bed.py
[ -d repeatcraftp ] || { git clone -q https://github.com/niccw/repeatcraftp.git; git -C repeatcraftp checkout -q eb7a30d; }
count RepeatMasker_out $G.fa.out

# 2. Overlap resolution: keep the higher-scoring annotation where annotations overlap
/usr/bin/time -l $TEATREE -i $G.fa.out -o tt_overlap -remove_short_fragments 0 > /dev/null 2> tt_overlap.time
count TEaTree_overlap tt_overlap.bed tt_overlap.time
/usr/bin/time -l $PY RM2Bed.py -d . -p rm2bed -o higher_score $G.fa.out 2> rm2bed.time
count RM2Bed_overlap rm2bed_rm.bed rm2bed.time

# 3. Defragmentation only (overlaps left in): merge fragments of the same family that are
#    < 150 bp apart and do not overlap on the consensus
/usr/bin/time -l $TEATREE -i $G.fa.out -o tt_defrag -mode alignment -mergemode threshold -gapsize 150 \
    -leave_overlap -remove_short_fragments 0 > /dev/null 2> tt_defrag.time
count TEaTree_defrag tt_defrag.merged.gff tt_defrag.time
rm -f tt_defrag.gff tt_defrag.label.gff
#    RepeatCraft needs the RepeatMasker GFF, so make it from the .out file
awk 'BEGIN{OFS="\t"} NR>3 {s=($9=="C")?"-":"+"; ts=($9=="C")?$14:$12;
     print $5,"RepeatMasker","similarity",$6,$7,$2,s,".","Target \"Motif:"$10"\" "ts" "$13}' $G.fa.out > rm.gff
printf 'shortTE_size: 100\nmapfile: None\nltr_finder_gff: None\nmax_LTR_size: 10000\nLTR_flanking_size: 200\ngap_size: 150\n' > repeatcraft.cfg
/usr/bin/time -l python3 repeatcraftp/repeatcraft.py -r rm.gff -u $G.fa.out -c repeatcraft.cfg -o repeatcraft -m loose 2> repeatcraft.time
count RepeatCraft_defrag repeatcraft.rmerge.gff repeatcraft.time
rm -f rm.gff repeatcraft.rclabel.gff

# 4. Not compared: TEaTree alignment mode with default settings (overlap resolution + defragmentation)
/usr/bin/time -l $TEATREE -i $G.fa.out -o tt_full -mode alignment > /dev/null 2> tt_full.time
count TEaTree_full_alignment tt_full.merged.gff tt_full.time
rm -f tt_full.gff tt_full.label.gff

column -t "$R"
