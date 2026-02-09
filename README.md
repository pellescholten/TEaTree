# Welcome to TEaTree
<p align="center">
<img width="517" height="440" alt="Image" src="https://github.com/user-attachments/assets/1d79a312-480f-428d-9f2f-12a1be981dd2" />
</p>
## About

TEaTree is a tool for refining transposable element and repeat annotations from RepeatMasker. RepeatMasker often contains many overlapping and fragmented TE annotations. For some analyses, non-overlapping and/or defragmented annotations are preferred or required. TEatree first uses an highly efficient interval tree algorithm find, sort and resolve overlapping TE annotations based on their genomic positions and Smith-Waterman score. TEaTree appropriately deals with chimeras and multiple overlapping sequences (>2) while also reassigning newly formed fragments that are under a certain length. Furthermore, TEaTree can merge fragmented copies based on their consensus family ID and consensus family position, and their RepeatMasker ID and/or sequence position. 

It also has an optional alignment mode that produces additional GFF files, including a defragmented merged GFF, useful for aligning genomic copies to consensus sequences.

A short animation explains the overlapping resolving and defragmentation (https://github.com/pellescholten/TEaTree/blob/main/TEaTree.mp4).
https://github.com/pellescholten/TEaTree/assets/126644559/c0fb8bea-5841-48d0-b62d-dafbdd05b08f
  
## Requirement
- Python 3.6 or later.
- Python modules (all should be built-in): os, sys, gzip, datetime, collections, argparse, errno
  
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
  
TEatree takes RepeatMasker output files as input, either in .out or .align format (i.e. genome.fa.out or genome.fa.align file). If the input is in .align format (e.g. example.align), TEaTree will automatically created a cleaned up file containing only the relevant information in the proper format (e.g. clean_example.align).

## Output files
A normal TEAtree run will generate a`.bed` file with non-overlapping and defragmented TE annotations. This file can for example be used to extract TE content.

- `.bed`  
The 4th column (name field) in the bed file will be `RM_n`.`repeat`.`repeat_class`.`original_id`, which is the same naming convention as the `gene_id` in `.gtf.gz` file.  
The 5th column (score field) contains the SW-score (1st column of the input `.fa.out` file).  
The 6th column shows the strandedness of the TE annotation. 
Zero based base pair positions are used.

If `-mode alignment` is specifed, three files will be created. These can be used for example for the re-alignment of the new annotations to the concensus sequences.
- `.gff`  
Same information as the bed file but in different format and where all overlapping annotations are resolved except for chimers (resolving these chimers would lead to odd re-alignments). 
- `.label.gff`  
The `.gff` but with labels for groups that will be merged/defragmented.
- `.merged.gff`  
The final `.gff` where annotations defragmented/merged.

Notes:
- The BED output is the primary nonoverlapping annotation set intended for counting and summarisation.
- The GFF outputs in alignment mode are intended for alignment to consensus and defragmentation workflows. The plain GFF output does not resolve overlapping chimeras in the same way as the BED output, so prefer the BED for counting.
  
## Options that affect output and their explanations
- `-remove_simple_repeat`  
By default, the script will keep "Simple_repeat" and "Low_complexity" in output files.
If you want to remove such annotations, please add the `-remove_simple_repeat` option.  

- `-min_length` (default = 50)  
Specify the minimum size that fragments must have after being cut to resolve potential overlap, to prevent very short fragments from being formed. This is also the minimum framgent length of elements if the `-remove_short_fragments` option is specified.

- `-annotation_overlap` (default = 80)  
 Sometimes certain parts of a sequence are annotated as belonging to multiple repeats. With this option, you can set the percentage of overlap allowed for two annotations. TEaTree will remove annotations that are contained for more than 80% in an annotation with a higher SW score. 

- `-mode` (default = TEcontent)  
Choose the output you want to get. Choose `TEcontent" `(default) for output where all overlaps are resolved; best suited for looking at the raw TE contents themselves. Choose `alignment` for output in .gff format that is best suited The annotations for re-alignment of annotations to their concensus sequence. When the mode is `alignment`, annotations will be merged (i.e. `-mergemode` will be set to `both` by default).

- `-mergemode` (default = none)  
 Annotations can be fragmented by RepeatMasker or simply not detected as belonging to the same annotation. By specifying a 'mergemode' (or by choosing 'alignment mode'), TEaTree will defragment some of these annotations. How this is done can be determined with this option. First of all, mergers are also required to belong to the same consensus family. If 'ID' is specified, the program will only consider for merging those elements that have the same RepeatMasker id. If 'threshold' is specified, the program will consider annotations that have a distance between them lower than a certain threshold. This threshold can be determined with the option -gapsize. If 'both' is specified, the program will use both options to consider potential mergers. Finally, potential mergers are checked for their position on the family consensus sequence. If the repeats could form a clean whole defragmented repeat, also considering the consensus sequence, the repeats are merged into one. 

- `-gapsize` (default = 150)  
Specify the threshold used for considering annotations for merging when alignment files are being made and -mergemode is 'threshold' or 'both'.

- `-remove_short_fragments` (default = 50)  
Short fragments can often be unreliable or not relevant to your research question. Call this option if you want to remove short fragments. Fragments will be removed at the start of the pipeline before checking overlaps. Remove fragments shorter than this length (default: 50). Use 0 to disable."

- `-leave_overlap` (default = False)  
Call if you do not want to resolve overlapping TE annotations. Can be used for simple reformatting or for TE defragmentation if mergemode is specified.

- `-allowed_consensus_overlap` (default = False)
Merging annotations is an art and a bit subjective. In theory, annotations that have overlapping corresponding sequence on their consensus family are never to be merged as they cannot possibly be from the same TE insertion. However, small mistakes can be made anywhere in the TE annotation or consensus family identification pipelines. Therefore, this option allows for more lenient behaviour towards TEs with overlapping consensus family sequences by either allowing a fixed amount of overlapping basepairs or allowing a proportion of overlapping basepairs relative to the length of the consensus family lenght. Provide an integer as absolute overlap length allowed or an integer + "p" (e.g. 2p) when the amount of consensus overlap allowed is a percentage of the consensus family length (recommend not more than 5p).

- `-family_filtering` (default = none)
There are several options to remove annotations based on a minimal absolute length. This option allows for filtering based on the length of an repeat relative to the length of its family consensus sequence. For this you can provide the path to a file that contains the names of the consensus families with their absolute length in tsv format. This file can be created using the script "get_fam_length.sh". This script takes a fasta file containing all consensus family sequences and provides this tsv. 
TEaTree than checks for each annotation whether it is longer or shorter than 30% (30 is default, different percentage can be set with `-family_filtering_length`). Only annotations longer than this family specific threshold are retained. LTR annotations receive an additional exception: any LTR fragment longer than 250 bp is kept even if it falls below the threshold. Simple repeat like entries whose family name starts with "(" are discarded. 

- `-family_filtering_length` (default = 30)
Set the threshold for family length based filtering (see `-family_filtering`).
  

