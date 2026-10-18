#!/usr/bin/env python3
"""Append one line of metrics for a repeat annotation file to a results table.
Usage: compare.py results.tsv TOOL_NAME FILE [TIME_FILE]
FILE can be a RepeatMasker .out, a TEaTree .bed, an RM2Bed .bed or a .gff (TEaTree/RepeatCraft).
TIME_FILE is the output of /usr/bin/time -l (macOS)."""
import os, re, sys
from collections import defaultdict

results, tool, path = sys.argv[1:4]
CLASSES = ["LINE", "SINE", "LTR", "DNA", "RC", "Satellite", "Simple_Low", "Other"]

# read every annotation as (chromosome, start, end, class), 0-based half-open
intervals = defaultdict(list)
class_bp = defaultdict(int)
for i, line in enumerate(open(path)):
    if line.startswith("#") or not line.strip():
        continue
    p = line.split() if path.endswith(".out") else line.rstrip("\n").split("\t")
    if path.endswith(".out"):        # RepeatMasker: skip 3 header lines
        if i < 3: continue
        chrom, start, end, cls = p[4], int(p[5]) - 1, int(p[6]), p[10]
    elif path.endswith(".gff"):      # TEaTree / RepeatCraft: class in column 3
        chrom, start, end, cls = p[0], int(p[3]) - 1, int(p[4]), p[2]
    elif "rm2bed" in path:           # RM2Bed: class in column 7
        chrom, start, end, cls = p[0], int(p[1]), int(p[2]), p[6]
    else:                            # TEaTree bed: name = RM_n.repeat.class.id
        chrom, start, end, cls = p[0], int(p[1]), int(p[2]), p[3].rsplit(".", 2)[1]
    intervals[chrom].append((start, end))
    top = cls.split("/")[0].rstrip("?")
    if top in ("Simple_repeat", "Low_complexity"): top = "Simple_Low"
    class_bp[top if top in CLASSES else "Other"] += end - start

# sweep along each chromosome: covered bp, and annotations overlapping another one
n = summed = covered = n_overlap = 0
for chrom, iv in intervals.items():
    iv.sort()
    cluster_start, cluster_end, cluster_size = iv[0][0], iv[0][1], 1
    for start, end in iv[1:] + [(float("inf"), float("inf"))]:
        if start < cluster_end:      # overlaps the current cluster
            cluster_end = max(cluster_end, end); cluster_size += 1
        else:                        # close the cluster
            covered += cluster_end - cluster_start
            n_overlap += cluster_size if cluster_size > 1 else 0
            cluster_start, cluster_end, cluster_size = start, end, 1
    n += len(iv); summed += sum(e - s for s, e in iv)

wall = ram = "NA"
if len(sys.argv) > 4:
    t = open(sys.argv[4]).read()
    wall = re.search(r"([\d.]+) real", t).group(1)
    ram = round(int(re.search(r"(\d+)\s+maximum resident set size", t).group(1)) / 1e9, 2)

header = ["tool", "annotations", "overlapping_another", "summed_bp", "covered_bp", "counted_twice_bp", "wall_s", "peak_ram_gb"] + \
         [c + "_bp" for c in CLASSES]
row = [tool, n, n_overlap, summed, covered, summed - covered, wall, ram] + [class_bp[c] for c in CLASSES]
new = not os.path.exists(results)
with open(results, "a") as out:
    if new: out.write("\t".join(header) + "\n")
    out.write("\t".join(map(str, row)) + "\n")
