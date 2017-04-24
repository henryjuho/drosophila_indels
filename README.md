# *Drosophila* INDEL analysis pipeline
Henry Juho Barton  
Department of Animal and Plant Sciences, The University of Sheffield  

# Introduction

This document outlines the pipeline used to generate and analyse an INDEL dataset from 17 *Drosophila melanogaster* individuals (ref) and 42 *Drosophila simulans* individuals (ref). The document is subdivided by processing steps.

## Programs and versions required

  * Python 2.7.2
  * SAMtools version 1.2
  * GATK version 3.7
  * bcftools version 1.3
  * VCFtools version 0.1.14
  * RepeatMasker version open-4.0.6
  * bedtools version 2.26.0

## Scripts used in this pipeline

\* Note \* that most scripts make use of the python qsub wrapper module ```qsub.py``` described here: <https://github.com/henryjuho/python_qsub_wrapper>.


|                            |                            |                                |                             |
|:---------------------------|:---------------------------|:-------------------------------|:----------------------------|
| qsub.py                    | split_bams.py              | merge_chromosomal_bams.py      | haplotype_caller.py         |
| samtools_calling.py        | genotypeGVCFs.py           | get_consensus_vcf.py           | get_mean_depth.py           |
| depth_filter.py            | filter_length_biallelic.py | rename_dsim_headers.py         | repeat_masking.py           |
| rm_out2bed.py              | repeat_filtering.py        | hardfilter.py                  | VQSR.py                     |
| exclude_snp_in_indel.py    | fasta_add_header_prefix.py | wholegenome_lastz_chain_net.py |     |

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
## Reference preparation

Reference chromosome order files were generate for correct position sorting for GATK as follows:

```
$ grep ^'>' dmel-all-chromosome-r5.34.fa | cut -d ' ' -f 1 | cut -d '>' -f 2 > dmel_chr_order.txt
$ grep ^'>' dsimV2-Mar2012.fa | cut -d ' ' -f 1 | cut -d '>' -f 2 > dsim_chr_order.txt
```
## Preparing BAM files

Multi-sample chromosomal BAM files were converted to single sample whole genome BAM files as follows:

```
$ split_bams.py -bam_dir /fastdata/bop15hjb/drosophila_data/dmel/ -out_dir /fastdata/bop15hjb/drosophila_data/dmel/unmerged_bams/ -evolgen
$ merge_chromosomal_bams.py -bam_dir /fastdata/bop15hjb/drosophila_data/dmel/unmerged_bams/ -out_dir /fastdata/bop15hjb/drosophila_data/dmel/individual_bams/ -evolgen
$ cd /fastdata/bop15hjb/drosophila_data/dmel/bams/individual_bams
$ ls *bam | while read i; do samtools index -b $i; done
$ ls $PWD/*bam > dmel_bam_list.txt

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
$ ls /fastdata/bop15hjb/drosophila_data/dmel/gatk_calling/gvcf/*vcf > /fastdata/bop15hjb/drosophila_data/dmel/gatk_calling/gvcf/dmel_gvcf.list
$ genotypeGVCFs.py -ref /fastdata/bop15hjb/drosophila_data/dmel_ref/dmel-all-chromosome-r5.34.fa -vcf_list /fastdata/bop15hjb/drosophila_data/dmel/gatk_calling/gvcf/dmel_gvcf.list -out /fastdata/bop15hjb/drosophila_data/dmel/gatk_calling/allsites/dmel_17flys.gatk.allsites.vcf -evolgen

$ haplotype_caller.py -ref /fastdata/bop15hjb/drosophila_data/dsim_ref/dsimV2-Mar2012.fa -bam_dir /fastdata/bop15hjb/drosophila_data/dsim/bams/individual_bams/ -out_dir /fastdata/bop15hjb/drosophila_data/dsim/gatk_calling/gvcf/ -evolgen
$ ls /fastdata/bop15hjb/drosophila_data/dsim/gatk_calling/gvcf/*vcf > /fastdata/bop15hjb/drosophila_data/dsim/gatk_calling/gvcf/dsim_gvcf.list
$ genotypeGVCFs.py -ref /fastdata/bop15hjb/drosophila_data/dsim_ref/dsimV2-Mar2012.fa -vcf_list /fastdata/bop15hjb/drosophila_data/dsim/gatk_calling/gvcf/dsim_gvcf.list -out /fastdata/bop15hjb/drosophila_data/dsim/gatk_calling/allsites/dsim_42flys.gatk.allsites.vcf -evolgen
```

SAMtools calling proceeded as follows:

```
$ samtools_calling.py -bam_list /fastdata/bop15hjb/drosophila_data/dmel/bams/individual_bams/dmel_bam_list.txt -ref /fastdata/bop15hjb/drosophila_data/dmel_ref/dmel-all-chromosome-r5.34.fa -out /fastdata/bop15hjb/drosophila_data/dmel/samtools_calling/dmel_17flys -evolgen
$ samtools_calling.py -bam_list /fastdata/bop15hjb/drosophila_data/dsim/bams/individual_bams/dsim_bam_list.txt -ref /fastdata/bop15hjb/drosophila_data/dsim_ref/dsimV2-Mar2012.fa -out /fastdata/bop15hjb/drosophila_data/dsim/samtools_calling/dsim_42flys -per_chr
```

## VQSR

### Generating sets of 'known variants'

In order to run GATKs variant quality score recalibration (VQSR) a set of high confidence variants was generated through the intersection of GATK and SAMtools raw variant calls. This 'truth set' was generated as follows:

```
$ get_consensus_vcf.py -vcf_I /fastdata/bop15hjb/drosophila_data/dmel/gatk_calling/allsites/dmel_17flys.gatk.allsites.vcf -vcf_II /fastdata/bop15hjb/drosophila_data/dmel/samtools_calling/dmel_17flys.samtools.allsites.vcf -ref /fastdata/bop15hjb/drosophila_data/dmel_ref/dmel-all-chromosome-r5.34.fa -out /fastdata/bop15hjb/drosophila_data/dmel/consensus/
$ get_consensus_vcf.py -vcf_I /fastdata/bop15hjb/drosophila_data/dsim/gatk_calling/allsites/dsim_42flys.gatk.allsites.vcf -vcf_II /fastdata/bop15hjb/drosophila_data/dsim/samtools_calling/dsim_42flys.samtools.allsites.vcf -ref /fastdata/bop15hjb/drosophila_data/dsim_ref/dsimV2-Mar2012.fa -out /fastdata/bop15hjb/drosophila_data/dsim/consensus/
```

Mean depth was calculated from the GATK allsites vcf using vcftools (ref) as follows:

```
$ get_mean_depth.py -vcf /fastdata/bop15hjb/drosophila_data/dmel/gatk_calling/allsites/dmel_17flys.gatk.allsites.vcf
$ cat /fastdata/bop15hjb/drosophila_data/dmel/gatk_calling/allsites/*idepth | grep -v ^I | cut -f 3 | awk '{sum+=$1} END {print sum/NR}'

$ get_mean_depth.py -vcf /fastdata/bop15hjb/drosophila_data/dsim/gatk_calling/allsites/dsim_42flys.gatk.allsites.vcf
$ cat /fastdata/bop15hjb/drosophila_data/dsim/gatk_calling/allsites/*idepth | grep -v ^I | cut -f 3 | awk '{sum+=$1} END {print sum/NR}'
```

| Species           | Mean depth  |
|:------------------|:-----------:|
| _D. melanogaster_ | 20x         |
| _D. simulans_     | 46x         |

Consensus INDEL and SNP vcfs were then hardfiltered to remove sites with less than half or more than twice the mean depth, multiallelic sites and INDELs greater than 50bp.

```
$ depth_filter.py -vcf /fastdata/bop15hjb/drosophila_data/dmel/consensus/dmel_17flys.consensus.raw.indels.vcf -mean_depth 20 -N 17
$ depth_filter.py -vcf /fastdata/bop15hjb/drosophila_data/dmel/consensus/dmel_17flys.consensus.raw.snps.vcf -mean_depth 20 -N 17
$ ls /fastdata/bop15hjb/drosophila_data/dmel/consensus/*dpfiltered.vcf | while read i; do filter_length_biallelic.py -vcf $i -ref /fastdata/bop15hjb/drosophila_data/dmel_ref/dmel-all-chromosome-r5.34.fa; done

$ depth_filter.py -vcf /fastdata/bop15hjb/drosophila_data/dsim/consensus/dsim_42flys.consensus.raw.indels.vcf -mean_depth 46 -N 42
$ depth_filter.py -vcf /fastdata/bop15hjb/drosophila_data/dsim/consensus/dsim_42flys.consensus.raw.snps.vcf -mean_depth 46 -N 42
$ ls /fastdata/bop15hjb/drosophila_data/dsim/consensus/*dpfiltered.vcf | while read i; do filter_length_biallelic.py -vcf $i -ref /fastdata/bop15hjb/drosophila_data/dsim_ref/dsimV2-Mar2012.fa; done
```

Bed files of repeat positions were obtained from repeat masker outputs (see Whole genome alignment) and sites falling within repeat regions were excluded as follows: 

```
$ cd /fastdata/bop15hjb/drosophila_data/wga/genomes
$ cat dmel-all-chromosome-r5.34.fa.out | rm_out2bed.py | bedtools sort -faidx ../../dmel_ref/dmel_chr_order.txt > dmel-all-chromosome-r5.34.fa.out.bed
$ grep -v NODE  dsimV2-Mar2012.rename.fa.out | rm_out2bed.py | bedtools sort -faidx ../../dsim_ref/dsim_chr_order.txt > dsimV2-Mar2012.rename.fa.out.bed
$ cd
$ ls /fastdata/bop15hjb/drosophila_data/dmel/consensus/*bial.vcf | while read i; do repeat_filtering.py -vcf $i -ref /fastdata/bop15hjb/drosophila_data/dmel_ref/dmel-all-chromosome-r5.34.fa -bed /fastdata/bop15hjb/drosophila_data/wga/genomes/dmel-all-chromosome-r5.34.fa.out.bed -evolgen ; done
$ ls /fastdata/bop15hjb/drosophila_data/dsim/consensus/*bial.vcf | while read i; do repeat_filtering.py -vcf $i -ref /fastdata/bop15hjb/drosophila_data/dsim_ref/dsimV2-Mar2012.fa -bed /fastdata/bop15hjb/drosophila_data/wga/genomes/dsimV2-Mar2012.rename.fa.out.bed -evolgen ; done
```

Finally sites were hard filtered according to the GATK best practice (QD < 2.0, ReadPosRankSum < -20.0, FS > 200.0, SOR > 10.0 for INDELs, QD < 2.0, MQ < 40.0, FS > 60.0, SOR > 3.0, MQRankSum < -12.5, ReadPosRankSum < -8.0, for SNPs <https://software.broadinstitute.org/gatk/guide/article?id=3225>), additionally, SNPs falling within INDELs were removed.

```
$ hardfilter.py -ref /fastdata/bop15hjb/drosophila_data/dmel_ref/dmel-all-chromosome-r5.34.fa -vcf /fastdata/bop15hjb/drosophila_data/dmel/consensus/dmel_17flys.consensus.raw.indels.dpfiltered.50bp_max.bial.rfiltered.vcf -mode INDEL -evolgen
$ hardfilter.py -ref /fastdata/bop15hjb/drosophila_data/dmel_ref/dmel-all-chromosome-r5.34.fa -vcf /fastdata/bop15hjb/drosophila_data/dmel/consensus/dmel_17flys.consensus.raw.snps.dpfiltered.50bp_max.bial.rfiltered.vcf -mode SNP -evolgen
$ exclude_snp_in_indel.py -vcf /fastdata/bop15hjb/drosophila_data/dmel/consensus/dmel_17flys.consensus.raw.snps.dpfiltered.50bp_max.bial.rfiltered.hfiltered.vcf

$ hardfilter.py -ref /fastdata/bop15hjb/drosophila_data/dsim_ref/dsimV2-Mar2012.fa -vcf /fastdata/bop15hjb/drosophila_data/dsim/consensus/dsim_42flys.consensus.raw.indels.dpfiltered.50bp_max.bial.rfiltered.vcf -mode INDEL -evolgen
$ hardfilter.py -ref /fastdata/bop15hjb/drosophila_data/dsim_ref/dsimV2-Mar2012.fa -vcf /fastdata/bop15hjb/drosophila_data/dsim/consensus/dsim_42flys.consensus.raw.snps.dpfiltered.50bp_max.bial.rfiltered.vcf -mode SNP -evolgen
$ exclude_snp_in_indel.py -vcf /fastdata/bop15hjb/drosophila_data/dsim/consensus/dsim_42flys.consensus.raw.snps.dpfiltered.50bp_max.bial.rfiltered.hfiltered.vcf
```

### 'Truth set' summary

| Filtering step     | _D. mel_ INDELs  | _D. sim_ INDELs | _D. mel_ SNPs  | _D. sim_ SNPs  |
|:-------------------|:----------------:|:---------------:|:--------------:|:--------------:|
| raw GATK           | 798107           | 2066449         | 6161265        | 17259553       |
| raw SAMtools       | 791236           | 1734626         | 3418572        | 10805464       |
| raw consensus      | 550354           | 1113349         | 3316428        | 10174694       |
| depth              | 437145           | 1098950         | 2628415        | 10079794       |
| length, allele no. | 331846           | 784782          | 2471023        | 8258471        |
| repeats            | 286450           | 713285          | 2319409        | 7889884        |
| hardfilters        | 286177           | 712647          | 2036210        | 6732750        |
| no snps in indels  | -                | -               | 2017080        | 6664891        |

### VQSR

VQSR was performed for SNP (for SNPs, variants in INDELs were first exluded from the raw data) and INDELs separately for both species as follows:

```
$ VQSR.py -vcf /fastdata/bop15hjb/drosophila_data/dmel/consensus/dmel_17flys.gatk.raw.indels.vcf -ref /fastdata/bop15hjb/drosophila_data/dmel_ref/dmel-all-chromosome-r5.34.fa -truth_set /fastdata/bop15hjb/drosophila_data/dmel/consensus/dmel_17flys.consensus.raw.indels.dpfiltered.50bp_max.bial.rfiltered.hfiltered.vcf -mode INDEL -out /fastdata/bop15hjb/drosophila_data/dmel/vqsr/
$ exclude_snp_in_indel.py -vcf /fastdata/bop15hjb/drosophila_data/dmel/consensus/dmel_17flys.gatk.raw.snps.vcf
$ VQSR.py -vcf /fastdata/bop15hjb/drosophila_data/dmel/consensus/dmel_17flys.gatk.raw.snps.exsnpindel.vcf -ref /fastdata/bop15hjb/drosophila_data/dmel_ref/dmel-all-chromosome-r5.34.fa -truth_set /fastdata/bop15hjb/drosophila_data/dmel/consensus/dmel_17flys.consensus.raw.snps.dpfiltered.50bp_max.bial.rfiltered.hfiltered.exsnpindel.vcf -mode SNP -out /fastdata/bop15hjb/drosophila_data/dmel/vqsr/

$ VQSR.py -vcf /fastdata/bop15hjb/drosophila_data/dsim/consensus/dsim_42flys.gatk.raw.indels.vcf -ref /fastdata/bop15hjb/drosophila_data/dsim_ref/dsimV2-Mar2012.fa -truth_set /fastdata/bop15hjb/drosophila_data/dsim/consensus/dsim_42flys.consensus.raw.indels.dpfiltered.50bp_max.bial.rfiltered.hfiltered.vcf -mode INDEL -out /fastdata/bop15hjb/drosophila_data/dsim/vqsr/
$ exclude_snp_in_indel.py -vcf /fastdata/bop15hjb/drosophila_data/dsim/consensus/dsim_42flys.gatk.raw.snps.vcf
$ VQSR.py -vcf /fastdata/bop15hjb/drosophila_data/dsim/consensus/dsim_42flys.gatk.raw.snps.exsnpindel.vcf -ref /fastdata/bop15hjb/drosophila_data/dsim_ref/dsimV2-Mar2012.fa -truth_set /fastdata/bop15hjb/drosophila_data/dsim/consensus/dsim_42flys.consensus.raw.snps.dpfiltered.50bp_max.bial.rfiltered.hfiltered.exsnpindel.vcf -mode SNP -out /fastdata/bop15hjb/drosophila_data/dsim/vqsr/
```

### Post VQSR filters

Variants more than twice and half the mean coverage, longer than 50bp, multiallelic or located in repeat regions were then removed from the post VQSR (95% cutoff) data.

```
$ cd /fastdata/bop15hjb/drosophila_data/dmel/vqsr/
$ ls $PWD/*t95.0.pass.vcf | while read i; do depth_filter.py -vcf $i -mean_depth 20 -N 17; done
$ ls $PWD/*dpfiltered.vcf | while read i; do filter_length_biallelic.py -vcf $i -ref /fastdata/bop15hjb/drosophila_data/dmel_ref/dmel-all-chromosome-r5.34.fa; done 
$ ls $PWD/*bial.vcf | while read i; do repeat_filtering.py -vcf $i -ref /fastdata/bop15hjb/drosophila_data/dmel_ref/dmel-all-chromosome-r5.34.fa -bed /fastdata/bop15hjb/drosophila_data/wga/genomes/dmel-all-chromosome-r5.34.fa.out.bed -evolgen ; done

$ cd /fastdata/bop15hjb/drosophila_data/dsim/vqsr/
$ ls $PWD/*t95.0.pass.vcf | while read i; do depth_filter.py -vcf $i -mean_depth 46 -N 42; done
$ ls $PWD/*dpfiltered.vcf | while read i; do filter_length_biallelic.py -vcf $i -ref /fastdata/bop15hjb/drosophila_data/dsim_ref/dsimV2-Mar2012.fa; done
$ ls $PWD/*bial.vcf | while read i; do repeat_filtering.py -vcf $i -ref /fastdata/bop15hjb/drosophila_data/dsim_ref/dsimV2-Mar2012.fa -bed /fastdata/bop15hjb/drosophila_data/wga/genomes/dsimV2-Mar2012.rename.fa.out.bed; done
```

### VQSR summary

| Filtering step     | _D. mel_ INDELs  | _D. sim_ INDELs | _D. mel_ SNPs  | _D. sim_ SNPs  |
|:-------------------|:----------------:|:---------------:|:--------------:|:--------------:|
| raw GATK           | 798107           | 2066449         | 6161265        | 17259553       |
| SNPs in INDELs     | -                | -               | 3357471        |  9786903       |
| post VQSR (95%)    | 651767           | 1576167         | 2416349        |  7703545       |
| depth              | 538762           | 1567823         | 2098303        |  7688121       |
| length, allele no. | 453181           | 1194893         | 2066044        |  7285990       |
| repeats            | 401692           | 1104572         | 1989033        |  7047712       |

## Whole genome alignment

Whole genome alignments were performed between _D. melanogaster_, _D. simulans_ and _D. yakuba_ using MultiZ (ref), following the UCSC pipeline (described here: ref).

First _D. simulans_ fasta headers were truncated to come within the required 50bp max length for RepeatMasker (ref).

```
$ cat dsimV2-Mar2012.fa | rename_dsim_headers.py > dsimV2-Mar2012.rename.fa 
```
Genomes were softmasked using RepeatMasker in the following script:

```
$ ls /fastdata/bop15hjb/drosophila_data/wga/genomes/*fa | while read i; do repeat_masking.py -fa $i -evolgen; done
```

The resulting soft masked fasta file headers were editted to contain species information:

```
$ cd /fastdata/bop15hjb/drosophila_data/wga/genomes/
$ fasta_add_header_prefix.py -fa dmel-all-chromosome-r5.34.fa.masked -pre 'dmel.' -truncate
$ fasta_add_header_prefix.py -fa dsimV2-Mar2012.rename.fa.masked -pre 'dsim.' -truncate
$ fasta_add_header_prefix.py -fa dyak-all-chromosome-r1.3.fa.masked -pre 'dyak.' -truncate
```

The resulting files were teh used to generate pairwise alignments with lastz (ref), which were then chainned and netted using x and y respectively.

```
$ wholegenome_lastz_chain_net.py -ref_fa /fastdata/bop15hjb/drosophila_data/wga/genomes/dmel-all-chromosome-r5.34.fa.masked.rename.fa -ref_name dmel -query_fa /fastdata/bop15hjb/drosophila_data/wga/genomes/dsimV2-Mar2012.rename.fa.masked.rename.fa -query_name dsim -out /fastdata/bop15hjb/drosophila_data/wga/pairwise_alignments/
$ wholegenome_lastz_chain_net.py -ref_fa /fastdata/bop15hjb/drosophila_data/wga/genomes/dmel-all-chromosome-r5.34.fa.masked.rename.fa -ref_name dmel -query_fa /fastdata/bop15hjb/drosophila_data/wga/genomes/dyak-all-chromosome-r1.3.fa.masked.rename.fa -query_name dyak -out /fastdata/bop15hjb/drosophila_data/wga/pairwise_alignments/
```

Single coverage was then ensured for the reference genome:

```
$ single_cov.py -dir /fastdata/bop15hjb/drosophila_data/wga/pairwise_alignments/ -ref_name dmel 
```
