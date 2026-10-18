# Benchmark

Here TEaTree is compared to RM2Bed.py and RepeatCraft, using the RepeatMasker annotation from UCSC as input.

## Run
```
cd benchmark
bash benchmark.sh dm6
bash benchmark.sh hg38
```
The script downloads the RepeatMasker `.out` file and the other tools, runs everything, and writes the results to `results/`. It is written for macOS. RM2Bed.py needs pandas. hg38 takes about 10 minutes and needs about 7 GB of disk.

## What is compared
1. **Resolving overlaps.** TEaTree (TEcontent mode) vs [RM2Bed.py](https://github.com/rmhubley/RepeatMasker/blob/master/util/RM2Bed.py) with `--ovlp_resolution higher_score`. Both keep the annotation with the higher score where annotations overlap.
2. **Merging fragments.** TEaTree (`-mode alignment -mergemode threshold -gapsize 150 -leave_overlap`) vs [RepeatCraft](https://github.com/niccw/repeatcraftp) (loose mode, `gap_size: 150`, without LTR_FINDER). Both merge fragments of the same family that are less than 150 bp apart, and do not resolve overlaps.

All tools get the same input and no tool removes short fragments. Everything ran on one thread on an Apple M3 Pro (18 GB RAM).

`count.py` counts for each output the number of annotations, the number of annotations that overlap another annotation, and the total and covered length (bp covered by at least one annotation).

## Results (hg38)
The hg38 file (downloaded 17 October 2026) has 5,520,118 annotations, of which 270,851 overlap another annotation.

| Resolving overlaps | Annotations | Overlapping another | Covered (Mb) | Time (s) | Memory (GB) |
|---|---:|---:|---:|---:|---:|
| RepeatMasker (input) | 5,520,118 | 270,851 | 1,586.3 | | |
| TEaTree | 5,492,039 | 0 | 1,585.6 | 31.8 | 0.25 |
| RM2Bed.py | 5,519,544 | 239,423 | 1,586.3 | 30.4 | 5.97 |

| Merging fragments | Annotations | Removed by merging | Time (s) | Memory (GB) |
|---|---:|---:|---:|---:|
| TEaTree | 5,301,120 | 218,998 | 96.0 | 4.24 |
| RepeatCraft | 5,202,386 | 317,732 | 402.0 | 7.44 |

TEaTree removes all overlaps, while RM2Bed.py leaves 239,423 overlapping annotations.

RepeatCraft merges more fragments than TEaTree because it checks less:
- it does not check if fragments are on the same strand;
- for the first two fragments of a group, it does not check if they overlap on the consensus sequence.

In dm6, 28% of the groups merged by RepeatCraft have fragments on opposite strands or fragments that overlap on the consensus, so they cannot come from the same TE insertion. TEaTree only merges fragments on the same strand, and allows at most 25 bp of overlap on the consensus (`-allowed_consensus_overlap`).

Results for dm6 are in `results/dm6.txt`.
