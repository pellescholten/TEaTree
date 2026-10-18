"""Count annotations, overlaps and covered bp in repeat annotation files.
Usage: python3 count.py FILE [FILE ...]      (RepeatMasker .out, .bed or .gff)"""
import sys
from collections import defaultdict

print("file\tannotations\toverlapping_another\tsummed_Mb\tcovered_Mb")
for path in sys.argv[1:]:
    # read all annotations per chromosome as (start, end), 0-based
    intervals = defaultdict(list)
    for i, line in enumerate(open(path)):
        p = line.split()
        if line.startswith("#") or len(p) < 5 or (path.endswith(".out") and i < 3):
            continue
        if path.endswith(".out"):
            intervals[p[4]].append((int(p[5]) - 1, int(p[6])))
        elif path.endswith(".gff"):
            intervals[p[0]].append((int(p[3]) - 1, int(p[4])))
        else:
            intervals[p[0]].append((int(p[1]), int(p[2])))

    # walk along each chromosome and join overlapping annotations into clusters
    n = summed = covered = overlapping = 0
    for iv in intervals.values():
        iv.sort()
        start, end, size = iv[0][0], iv[0][1], 1
        for s, e in iv[1:] + [(float("inf"), 0)]:
            if s < end:                        # overlaps the current cluster
                end, size = max(end, e), size + 1
            else:                              # cluster finished
                covered += end - start
                overlapping += size if size > 1 else 0
                start, end, size = s, e, 1
        n += len(iv)
        summed += sum(e - s for s, e in iv)
    print(f"{path}\t{n}\t{overlapping}\t{summed / 1e6:.1f}\t{covered / 1e6:.1f}")
