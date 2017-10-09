# Title     : alpha calulation for dmel INDELs
# Objective : calculates alpha
# Created by: henryjuho
# Created on: 09/10/2017

calc_alpha <- function(dn, ds, pn, ps) {
  alpha_val = 1 - ((ds * pn) / (dn * ps))
  return(alpha_val)
}

cds_indel_subs = read.delim('~/sharc_fastdata/drosophila_data/dmel/indel_divergence/dmel_cds_indel_divergence.txt')
non_coding_indel_subs = read.delim('~/sharc_fastdata/drosophila_data/dmel/indel_divergence/dmel_noncoding_indel_divergence.txt')

indel_poly = read.delim('~/sharc_fastdata/drosophila_data/dmel/summary_stats/dmel_17flys_indel_summary_no_bs_split_ar_nc.txt')

dn_data = subset(cds_indel_subs, chromo!='XHet' & chromo!='X' & chromo!='YHet' &
    chromo!='Y' & chromo!='U' & chromo!='dmel_mitochondrion_genome')

ds_data = subset(non_coding_indel_subs, chromo!='XHet' & chromo!='X' & chromo!='YHet' &
    chromo!='Y' & chromo!='U' & chromo!='dmel_mitochondrion_genome')

cds_poly = subset(indel_poly, bin=='CDS' & type=='indel')
non_coding_poly = subset(indel_poly, bin=='non-coding' & type=='indel')

d_n = sum(dn_data$indels)
d_s = sum(ds_data$indels)

p_n =  cds_poly$seg_sites
p_s = non_coding_poly$seg_sites

#print(c(d_n, d_s, p_n, p_s))

indel_alpha = calc_alpha(d_n, d_s, p_n, p_s)

print(indel_alpha)