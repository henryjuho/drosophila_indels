library(ggplot2)
library(dplyr)

cds_div = read.delim('divergence_data/dmel_cds_indel_divergence.txt')
non_coding_div = read.delim('divergence_data/dmel_noncoding_indel_divergence.txt')
ar_div = read.delim('divergence_data/dmel_ar_indel_divergence.txt')
line_div = read.delim('divergence_data/dmel_ar_LINEs_indel_divergence.txt')

cds_div$type = 'cds'
non_coding_div$type = 'non-coding'
ar_div$type = 'ar'
line_div$type = 'LINE'

div = rbind(cds_div, non_coding_div, ar_div, line_div)

div$type = factor(div$type, levels=c('cds', 'non-coding', 'ar', 'LINE'))

all_chr = summarise(group_by(subset(div, chromo != 'X' & chromo != 'XHet' & chromo != 'YHet'), type, indel_type),
    chromo='autosomes', indels=sum(n_variants), callable=sum(callable))

all_chr$divergence = all_chr$indels / all_chr$callable

div_plot = ggplot(all_chr, aes(x=type, y=divergence, fill=indel_type)) +
    geom_bar(stat='identity', position='dodge') +
    theme_bw(base_size = 10) +
    xlab('') + ylab('dvergence') +
    theme(legend.title=element_blank(), legend.position=c(0.2, 0.8), legend.background=element_blank())

ggsave('indel_divergence.pdf', plot=div_plot, width=3.5, height=3)

write.csv(all_chr, file='indel_divergence.csv', row.names=FALSE)