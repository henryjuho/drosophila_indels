# Title     : dmel indel lengths
# Objective : generates summary plot of indel lengths
# Created by: henryjuho
# Created on: 28/09/2017

library(ggplot2)
library(plyr)

length_data = read.delim('./dmel_indel_lengths.txt')

indel_only = subset(length_data, variant=='indel' & region!='non-coding')
indel_only$region = factor(indel_only$region, levels=c('gwide', 'CDS'))
indel_only$region <- revalue(x = indel_only$region, c("gwide" = "Genome wide", "CDS" = "Coding sequence"))

len_plot = ggplot(indel_only, aes(x=length, y=count))+
    geom_bar(stat='identity', position='dodge') +
    theme_bw(base_size = 11) +
    labs(x='INDEL length (bp)', y='Number of variants') +
    facet_wrap(~region, ncol=1, scales='free')

ggsave('dmel_lengths.pdf', height=5, width=5, plot=len_plot)
