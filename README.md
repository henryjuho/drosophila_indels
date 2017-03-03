# *Drosophila* INDEL analysis pipeline
Henry Juho Barton  
Department of Animal and Plant Sciences, The University of Sheffield  

# Introduction

This document outlines the pipeline used to generate and analyse an INDEL dataset from 17 *Drosophila melanogaster* individuals (ref) and 42 *Drosophila simulans* individuals (ref). The document is subdivided by processing steps.

## Programs and versions required

  * Python 2.7.2
  * SAMtools version 1.2
  * GATK version 3.7

## Scripts used in this pipeline

\* Note \* that most scripts make use of the python qsub wrapper module ```qsub.py``` described here: <https://github.com/henryjuho/python_qsub_wrapper>.


|                            |                            |                             |                             |
|:---------------------------|:---------------------------|:----------------------------|:----------------------------|
| qsub.py                    | split_bams.py              | merge_chromosomal_bams.py   | haplotype_caller.py         |
| samtools_calling.py        |                            |                             |                             |

## Reference and annotation files required for analysis

  * *D. melanogaster* reference genome: 
  * *D. simulans* reference genome:
  
## BAM files

| Region                     | _Drosophila melanogaster_   | _Drosophila simulans_     |
|:---------------------------|:---------------------------:|:-------------------------:|
| 2LHet                      | 2LHet.merged.realigned.bam  |    NA                     |
| 2L                         | 2L.merged.realigned.bam     | 2L_merged.realigned.bam   |
| 2RHet                      | 2RHet.merged.realigned.bam  |    NA                     |
| 2R                         | 2R.merged.realigned.bam     | 2R_merged.realigned.bam   |
| 3LHet                      | 3LHet.merged.realigned.bam  |    NA                     |
| 3L                         | 3L.merged.realigned.bam     | 3L_merged.realigned.bam   |
| 3RHet                      | 3RHet.merged.realigned.bam  |    NA                     |
| 3R                         | 3R.merged.realigned.bam     | 3R_merged.realigned.bam   |
| 4                          | 4.merged.realigned.bam      | 4_merged.realigned.bam    |
| XHet                       | XHet.merged.realigned.bam   |    NA                     |
| X                          | X.merged.realigned.bam      | X_merged.realigned.bam    |


# Data preparation pipeline
## Preparing BAM files

Multi-sample chromosomal BAM files were converted to single sample whole genome BAM files as follows:

```
$ split_bams.py -bam_dir /fastdata/bop15hjb/drosophila_data/dmel/ -out_dir /fastdata/bop15hjb/drosophila_data/dmel/unmerged_bams/ -evolgen
$ merge_chromosomal_bams.py -bam_dir /fastdata/bop15hjb/drosophila_data/dmel/unmerged_bams/ -out_dir /fastdata/bop15hjb/drosophila_data/dmel/individual_bams/ -evolgen
$ cd /fastdata/bop15hjb/drosophila_data/dmel/bams/individual_bams
$ ls *bam | while read i; do samtools index -b $i; done
$ ls $PWD/*bam > dmel_bam_list.txt
$
$ split_bams.py -bam_dir /fastdata/bop15hjb/drosophila_data/dsim/ -out_dir /fastdata/bop15hjb/drosophila_data/dsim/unmerged_bams/ -evolgen
$ merge_chromosomal_bams.py -bam_dir /fastdata/bop15hjb/drosophila_data/dsim/unmerged_bams/ -out_dir /fastdata/bop15hjb/drosophila_data/dsim/individual_bams/ -evolgen
$ cd /fastdata/bop15hjb/drosophila_data/dsim/bams/individual_bams
$ ls *bam | while read i; do samtools index -b $i; done
$ ls $PWD/*bam > dsim_bam_list.txt
```

## Variant calling

Raw SNPs and INDELs were called using both GATK's Haplotype Caller and SAMtools mpileup. GATK calling proceeded as follows:

```
$ haplotype_caller.py -ref /fastdata/bop15hjb/drosophila_data/dmel_ref/dmel-all-chromosome-r5.34.fa -bam_dir /fastdata/bop15hjb/drosophila_data/dmel/bams/individual_bams/ -out_dir /fastdata/bop15hjb/drosophila_data/dmel/gatk_calling/gvcf/
$ haplotype_caller.py -ref /fastdata/bop15hjb/drosophila_data/dsim_ref/dsimV2-Mar2012.fa -bam_dir /fastdata/bop15hjb/drosophila_data/dsim/bams/individual_bams/ -out_dir /fastdata/bop15hjb/drosophila_data/dsim/gatk_calling/gvcf/ -evolgen
```

SAMtools calling proceeded as follows:

```
$ samtools_calling.py -bam_list /fastdata/bop15hjb/drosophila_data/dmel/bams/individual_bams/dmel_bam_list.txt -ref /fastdata/bop15hjb/drosophila_data/dmel_ref/dmel-all-chromosome-r5.34.fa -out /fastdata/bop15hjb/drosophila_data/dmel/samtools_calling/dmel_17flys -evolgen
$ samtools_calling.py -bam_list /fastdata/bop15hjb/drosophila_data/dsim/bams/individual_bams/dsim_bam_list.txt -ref /fastdata/bop15hjb/drosophila_data/dsim_ref/dsimV2-Mar2012.fa -out /fastdata/bop15hjb/drosophila_data/dsim/samtools_calling/dsim_42flys -evolgen
```
