library(ggplot2)
library(dplyr)

cds_div = read.delim('/Users/henryjuho/sharc_fastdata/drosophila_data/dmel/indel_divergence/dmel_cds_indel_divergence.txt')
non_coding_div = read.delim('/Users/henryjuho/sharc_fastdata/drosophila_data/dmel/indel_divergence/dmel_noncoding_indel_divergence.txt')

cds_div$type = 'cds'
non_coding_div$type = 'non-coding'

div = rbind(cds_div, non_coding_div)

all_chr = summarise(group_by(subset(div, chromo != 'X' & chromo != 'XHet' & chromo != 'YHet'), type),
    indels=sum(indels), callable=sum(callable))


all_chr$div = all_chr$indels / all_chr$callable

div_plot = ggplot(all_chr, aes(x=type, y=div)) +
    geom_bar(stat='identity', position='dodge') +
    theme_bw(base_size = 10) +
    xlab('') + ylab('dvergence')

ggsave('indel_divergence.pdf', plot=div_plot, width=3, height=3)