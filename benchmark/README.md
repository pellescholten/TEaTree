# Benchmark

Compares TEaTree with [RM2Bed.py](https://github.com/rmhubley/RepeatMasker/blob/master/util/RM2Bed.py) (overlap resolution) and [RepeatCraft](https://github.com/niccw/repeatcraftp) (defragmentation) on the same RepeatMasker `.out` file from UCSC.

```bash
bash benchmark.sh dm6 /path/to/python_with_pandas    # about 1 minute
bash benchmark.sh hg38 /path/to/python_with_pandas   # about 10 minutes, needs about 7 GB of disk
```

- `benchmark.sh` downloads the input and the tools into `work_<assembly>/`, runs everything (macOS), and writes `results/<assembly>.txt`.
- `count.py` counts, for each output, the annotations, the annotations that overlap another one, and the summed and covered Mb.
- RM2Bed.py needs a Python with pandas.

The results in `results/` were run on an Apple M3 Pro with 18 GB RAM.
