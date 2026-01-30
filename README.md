# Welcome to TEaTree

TEaTree: essential oil for your TE annotation machine

## About

TEaTree is a tool for refining transposable element and repeat annotations from RepeatMasker. RepeatMasker often contains many overlapping and fragmented TE annotations. For some analyses, non-overlapping and/or defragmented annotations are preferred or required. TEatree first uses an highly efficient interval tree algorithm find, sort and resolve overlapping TE annotations based on their genomic positions and Smith-Waterman score. Furthermore, TEaTree can merge fragmented copies based on their RepeatMasker ID and/or family and position. 

It also has an optional alignment mode that produces additional GFF files, including a defragmented merged GFF, useful for aligning genomic copies to consensus sequences.
  
## Requirement
- Python 3.6 or later.
- Python modules (all should be built-in): os, sys, gzip, datetime, collections, argparse, errno
  
The script `collapseTree.py` will remove overlapping repeat annotations based on the Smith–Waterman score of their alignments (see video https://github.com/pellescholten/TEaTree/blob/main/TEaTree.mp4).

https://github.com/pellescholten/TEaTree/assets/126644559/c0fb8bea-5841-48d0-b62d-dafbdd05b08f

## Installation

Clone the repository and run the script directly.

```bash
git clone https://github.com/pellescholten/TEaTree.git
cd TEaTree
python TEaTree.py -h
```

## How to use
The minimal required flags to run TEaTree are `-i` for the input file (.out or .align) and `-o` for the output file basename.

```
python collapseTree.py \
-i example/example_input.out \
-o example/example_output
```
  
TEatree takes RepeatMasker output files as input, either in .out or .align format (i.e. genome.fa.out or genome.fa.align file).

## Output files
A normal TEAtree run will generate a`.bed` file with non-overlapping and defragmented TE annotations. This file can for example be used to extract TE content.

- `.bed`  
The 4th column (name field) in the bed file will be `RM_n`.`repeat`.`repeat_class`.`original_id`, which is the same naming convention as the `gene_id` in `.gtf.gz` file.  
The 5th column (score field) contains the SW-score (1st column of the input `.fa.out` file).  
The 6th column shows the strandedness of the TE annotation. 
Zero based base pair positions are used.

If `-alignment` is specifed, three files will be created. These can be used for example for the re-alignment of the new annotations to the concensus sequences.
- `.gff`  
Same information as the bed file but in different format and where all overlapping annotations are resolved except for chimers (resolving these chimers would lead to odd re-alignments). 
- `.label.gff`  
The `.gff` but with labels for groups that will be merged/defragmented.
- `.merged.gff`  
The final `.gff` where annotations defragmented/merged.

Notes:
- The BED output is the primary nonoverlapping annotation set intended for counting and summarisation.
- The GFF outputs in alignment mode are intended for alignment to consensus and defragmentation workflows. The plain GFF output does not resolve overlapping chimeras in the same way as the BED output, so prefer the BED for counting.
  
## Options that affect output results
- `-remove_simple_repeat`  
By default, the script will keep "Simple_repeat" and "Low_complexity" in output files.
If you want to remove such annotations, please add the `-remove_simple_repeat` option.  

  
- `-lvl [int]` (default = 80)  
 When looking at overlaps, the programs removes annotations that are very similar. Meaning that, if an element/annotation is contained for more lvl% in a higher scoring annotation, it will be removed. This similarity level can be set.


- `-alignment` (default = False)  
When this option is specified, two extra files .gff will be created. The annotations in these files can be used for alignment of annotations to their concensus. The annotations in these files are defragmented by merging TEs together that are of the same family, close together, and correctly positioned on their concensus sequence.

- `-mergemode` (default = both)  
This option is only relevant if the alignment option is specified. Annotations are defragmented. How this is done can be determined with this option. If 'ID' is specified, the program will only consider for merging those elements that have the same RepeatMasker id. If 'threshold is specified, the program will consider annotations that have a distance between them lower than a certain threshold. This threshold can be determined with the option -gapsize. If 'both' is specified, the program will use both options to consider potential merges. 


- `-gapsize` (default = 150)  
Specify the threshold used for considering annotations for merging when alignment files are being made and -mergemode is 'threshold' or 'both'.

- `-remove` (default = False)  
Specify this option if you want to remove short fragments. 

- `-min` (default = 50)  
Specify the minimum size that fragments must have after being cut to resolve potential overlap, to prevent very short fragments from being formed. This is also the minimum framgent length of elements if the -remove option is specified.


  

