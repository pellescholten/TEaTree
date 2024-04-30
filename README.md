# Welcome to collapse Tree

### About
The script `collapseTree.py` collapses TE annotations from RepeatMasker (i.e. genome.fa.out file). The TE annotations from RepeatMasker are frequently overlap each other. If so, when counting NGS reads mapping to TEs, reads mapping to two or more TE annotations may not be counted. In such case, it is preferable to use non-overlapping TE annotations. 
  
### Requirement
- Python 3.6 or later.
- Python modules (all should be built-in): os, sys, gzip, datetime, collections, argparse, errno
  
### How to use
If you do not have `genome.fa.out`, you need to generate it by RepeatMasker.
The command below will output `example/example_input.out` file, which contains repeat annotations.
```
# this is the example of GRCm38
RepeatMasker \
GRCm38.p6.genome.fa \
-species 10090 \
-s -no_is \
-dir ./RM_out \
-pa 8
```
  
The script `collapseTree.py` will remove overlapping repeat annotations based on bit score.
  
The two flags below are required. Please specify an input `.fa.out` file with the `-i` flag, and an output file basename with the `-o` flag.
  
- `-i [input.fa.out]`  
- `-o [output basename]`
  
```
python collapseTree.py \
-i example/example_input.out \
-o example/example_output
```
  
### Output files
It will generate a`.bed` file, which can be used to extract TE content.

- `.bed`  
The 4th column (name field) in the bed file will be `RM_n`.`repeat`.`repeat_class`.`original_id`, which is the same naming convention as the `gene_id` in `.gtf.gz` file.  
The 5th column (score filed) will have bit score (1st column of the input `.fa.out` file).  
Zero based base pair positions are used.


If `-alignment` is specifed, two 3 files will be created that can be used for alignments of the annotations to the concensus sequences.
- `.gff`  
Same information as the bed file but in different format and where overlapping chimers are not resolved.
- `.label.gff`  
The `.gff` but with labels for groups that will be merged.
- `.merged.gff`  
The `.gff` but defragmented.
  
### Options that affect output results
- `-remove_simple_repeat`  
By default, the script will keep "Simple_repeat" and "Low_complexity" in output files.
If you want to remove such annotations, please add the `-remove_simple_repeat` option.  
```
python collapseTree.py \
-i example/example_input.out \
-o example/example_output \
-remove_simple_repeat
```
  
- `-lvl [int]` (default = 80)  
 When looking at overlaps, the programs removes annotations that are very similar. Meaning that, if an element/annotation is contained for more lvl% in a higher scoring annotation, it will be removed. This similarity level can be set.
```
python collapseTree.py \
-i example/example_input.out \
-o example/example_output \
-lvl 90
```

- `-alignment` (default = False)  
When this option is specified, two extra files .gff will be created. The annotations in these files can be used for alignment of annotations to their concensus. The annotations in these files are defragmented by merging TEs together that are of the same family, close together, and correctly positioned on their concensus sequence.
```
python collapseTree.py \
-i example/example_input.out \
-o example/example_output \
-alignment
```

- `-mergemode` (default = both)  
This option is only relevant if the alignment option is specified. Annotations are defragmented. How this is done can be determined with this option. If 'ID' is specified, the program will only consider for merging those elements that have the same RepeatMasker id. If 'threshold is specified, the program will consider annotations that have a distance between them lower than a certain threshold. This threshold can be determined with the option -gapsize. If 'both' is specified, the program will use both options to consider potential merges. 
```
python collapseTree.py \
-i example/example_input.out \
-o example/example_output \
-mergemode ID
```

- `-gapsize` (default = 150)  
Specify the threshold used for considering annotations for merging when alignment files are being made and -mergemode is 'threshold' or 'both'.
```
python collapseTree.py \
-i example/example_input.out \
-o example/example_output \
-threshold 300
```
- `-remove` (default = False)  
Specify this option if you want to remove short fragments. Note that this will work differently for normal mode and alignment mode. In the normal mode for TE content, short fragments are removed at the very start so are also not considered in overlap. In the alignment mode, short fragments are removed at the very end so they can still be used for defragmentation. This may lead to slightly different overlap resolving. e.g.:
```
Removal after (for alignment)
200  ———        fix overlap     —>      ———     remove bottom fragment      —>      ——— remove short fragments  —> 0
100  ——————                     —>      ---———  that has become too short   —>      0                           —> 0
Removal before (for TE content)
200 ———         remove short fragments  —>  0        
100 ——————                              —>  ——————
```
```
python collapseTree.py \
-i example/example_input.out \
-o example/example_output \
-remove
```

- `-min` (default = 50)  
Specify the minimum size that fragments must have after being cut to resolve potential overlap, to prevent very short fragments from being formed. This is also the minimum framgent length of elements if the -remove option is specified.
```
python collapseTree.py \
-i example/example_input.out \
-o example/example_output \
-min 80
```
  
### Further options that does not affect results
- `-quiet`  
- `-version`  
- `-h`  
  
Please see further information by the `-h` option.  
```
python collapseTree.py -h
```


