# Benchmark

`benchmark.sh` compares TEaTree with other tools on the same RepeatMasker `.out` file (from UCSC). Each comparison uses matched settings:

| Task | TEaTree | Compared with |
|---|---|---|
| Overlap resolution (higher score wins) | `-remove_short_fragments 0` (TEcontent mode) | [RM2Bed.py](https://github.com/rmhubley/RepeatMasker/blob/master/util/RM2Bed.py) `--ovlp_resolution higher_score` |
| Defragmentation only (same family, < 150 bp apart, no overlap on the consensus) | `-mode alignment -mergemode threshold -gapsize 150 -leave_overlap -remove_short_fragments 0` | [RepeatCraft](https://github.com/niccw/repeatcraftp) loose mode, `gap_size: 150`, no LTR_FINDER input |

In both comparisons, all tools:
- get the same input;
- run one after the other on the same machine, on one thread, with the same Python;
- apply no length filter;
- have their outputs counted by the same script, `compare.py`.

TEaTree alignment mode with default settings (overlap resolution plus defragmentation) is also run, but it is not compared with another tool because no other tool does both.

## Run
```bash
bash benchmark.sh dm6 /path/to/python_with_pandas    # about 1 minute
bash benchmark.sh hg38 /path/to/python_with_pandas   # about 25 minutes, needs about 6 GB of disk
```
The input, the other tools and all outputs go to `work_<assembly>/`. The results go to `results/<assembly>.tsv`.
Requirements:
- macOS (`/usr/bin/time -l`), curl and git;
- a Python with pandas, for RM2Bed.py.

TEaTree, RepeatCraft and `compare.py` only use the Python standard library.

## Metrics (`compare.py`)
- `annotations`: number of annotations.
- `overlapping_another`: annotations that overlap at least one other annotation.
- `summed_bp`: summed length of all annotations.
- `covered_bp`: genome bp in at least one annotation.
- `counted_twice_bp`: `summed_bp - covered_bp`, i.e. bp counted more than once when TE content is summed over annotations.
- `wall_s`, `peak_ram_gb`: from `/usr/bin/time -l`.
- `<class>_bp`: summed length per repeat class.

The results in `results/` were run on an Apple M3 Pro with 18 GB RAM and Python 3.12.
