# *Drosophila* INDEL analysis pipeline
Henry Juho Barton  
Department of Animal and Plant Sciences, The University of Sheffield  

# Introduction

This document outlines the pipeline used to generate and analyse an INDEL dataset from 17 *Drosophila melanogaster* (ref) and x *Drosophila simulans* haplotypes. The document is subdivided by processing steps.

## Programs and versions required

  * Python 2.7.2
  * SAMtools version 1.2



## Scripts used in this pipeline

\* Note \* that most scripts make use of the python qsub wrapper module ```qsub.py``` described here: <https://github.com/henryjuho/python_qsub_wrapper>.


|                            |                            |                             |                             |
|:---------------------------|:---------------------------|:----------------------------|:----------------------------|
| qsub.py                    | split_bams.py              | merge_chromosomal_bams.py   |                             |

## Reference and annotation files required for analysis

  * *D. melanogaster* reference genome: 
  * *D. simulans* reference genome:
  
## BAM files

# Pipeline
## Preparing BAM files

Multi-sample chromosomal BAM files were converted to single sample whole genome BAM files as follows:

```
$ split_bams.py -bam_dir /fastdata/bop15hjb/drosophila_data/dmel/ -out_dir /fastdata/bop15hjb/drosophila_data/dmel/unmerged_bams/ -evolgen
$ merge_chromosomal_bams.py -bam_dir /fastdata/bop15hjb/drosophila_data/dmel/unmerged_bams/ -out_dir /fastdata/bop15hjb/drosophila_data/dmel/individual_bams/ -evolgen

$ split_bams.py -bam_dir /fastdata/bop15hjb/drosophila_data/dsim/ -out_dir /fastdata/bop15hjb/drosophila_data/dsim/unmerged_bams/ -evolgen
$ merge_chromosomal_bams.py -bam_dir /fastdata/bop15hjb/drosophila_data/dsim/unmerged_bams/ -out_dir /fastdata/bop15hjb/drosophila_data/dsim/individual_bams/ -evolgen
```