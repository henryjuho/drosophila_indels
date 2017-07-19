rm(list=ls())

library(ggplot2)
library(Hmisc)
library(dplyr)

cbPalette <- c("#E69F00", 'tomato 3', 'steel blue', "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

gene_data <- read.delim('/Users/henryjuho/iceberg_fastdata/drosophila_data/dmel/gene_analysis/dmel.dndn.pi0pi4.pi_indel.longest_trans.txt')

genes_300min <- subset(gene_data, dS<=5) #length>=300)

genes_300min$dN_bin <- as.numeric(cut2(genes_300min$dN, g=20))

summary_data_pi0 = summarise(group_by(genes_300min, dN_bin), dN_max = round(max(dN), 3), pi_mean = mean(pi0))
summary_data_pi0$pi_type = 'pi0'
summary_data_pi4 = summarise(group_by(genes_300min, dN_bin), dN_max = round(max(dN), 3), pi_mean = mean(pi4))
summary_data_pi4$pi_type = 'pi4'
summary_data_piI = summarise(group_by(genes_300min, dN_bin), dN_max = round(max(dN), 3), pi_mean = mean(pi_indel))
summary_data_piI$pi_type = 'pi_indel'

ggplot(rbind(summary_data_pi0, summary_data_pi4, summary_data_piI), aes(x=as.factor(dN_max), y=pi_mean, colour=pi_type)) +
  geom_point(stat='identity', size = 2.5) + ylim(0, 0.02) +
  theme_bw(base_size = 15) +
  scale_colour_manual(values=cbPalette) +
  ylab(expression(pi~mean)) + xlab('dN') +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle=90),
        legend.position = c(0.85, 0.85))

#ggsave('/Users/henryjuho/genomics/drosophila_indel_results/dmel.pi_dn.jpg', width=5, height=5)
