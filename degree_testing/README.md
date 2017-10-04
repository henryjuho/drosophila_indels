# anavar degree analysis

anavar was run with the following degree settings: 25, 50, 75, 100, 150, 200 and 300.

Each run followed the below command line:

```
cds_vs_neutral_anavar_snps.py -vcf /fastdata/bop15hjb/drosophila_data/dmel/analysis_ready_data/dmel_17flys.gatk.raw.snps.exsnpindel.recalibrated.filtered_t95.0.pass.dpfiltered.50bp_max.bial.rmarked.polarised.annotated.ar.degen.vcf.gz -n 17 -c 1 -dfe continuous -call_csv /fastdata/bop15hjb/drosophila_data/dmel_ref/dmel.callablesites.summary_with_degen.csv -neu_type 4fold -out_pre <out> -degree <degree>
```

The jobs were generated and submitted as follows:

```
$ drosophila_indels/degree_testing/anavar_degree_test.py
```