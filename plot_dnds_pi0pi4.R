library(ggplot2)
library(Hmisc)
library(dplyr)
library(gridExtra)
library(grid)

cbPalette <- c("#E69F00", 'tomato 3', 'steel blue', "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

gene_data <- read.delim('/Users/henryjuho/sharc_fastdata/drosophila_data/dmel/gene_analysis/dmel.gene_summarystats.longest_trans.txt')

genes_300min <- subset(gene_data, dS<=5) # & tajd0!=0.0 & tajd_indel!=0.0) #length>=300)

# replace 0.0 tajimas D values with NA
genes_300min[which(genes_300min$tajd0 == 0.0), 'tajd0'] = NA
genes_300min[which(genes_300min$tajd_indel == 0.0), 'tajd_indel'] = NA

genes_300min$dN_bin <- as.numeric(cut2(genes_300min$dN, g=20))

# pi plot
summary_data_pi0 = summarise(group_by(genes_300min, dN_bin), dN_max = round(max(dN), 3), pi_mean = mean(pi0))
summary_data_pi0$pi_type = 'zerofold SNPs'
# summary_data_pi4 = summarise(group_by(genes_300min, dN_bin), dN_max = round(max(dN), 3), pi_mean = mean(pi4))
# summary_data_pi4$pi_type = 'fourfold SNPs'
summary_data_piI = summarise(group_by(genes_300min, dN_bin), dN_max = round(max(dN), 3), pi_mean = mean(pi_indel))
summary_data_piI$pi_type = 'INDELs'

pi_plot = ggplot(rbind(summary_data_pi0, summary_data_piI), aes(x=as.factor(dN_max), y=pi_mean, colour=pi_type)) +
  geom_point(stat='identity', size = 2) + # ylim(0, 0.02) +
  theme_bw(base_size = 11) +
  scale_colour_manual(values=cbPalette) +
  ylab(expression(pi~mean)) + xlab(expression('d'[N])) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle=90),
        legend.position = c(0.25, 0.85),
        legend.background=element_blank())

pi_plot = arrangeGrob(pi_plot, top = textGrob("(A)", x = unit(0, "npc"),
          y = unit(0, "npc"), just=c("left","top"),
          gp=gpar(col="black", fontsize=11)))

# taj plot
summary_data_tajd0 = summarise(group_by(genes_300min, dN_bin), dN_max = round(max(dN), 3), taj_mean = mean(tajd0, na.rm=TRUE))
summary_data_tajd0$tajd_type = 'zerofold SNPs'
# summary_data_tajd4 = summarise(group_by(genes_300min, dN_bin), dN_max = round(max(dN), 3), taj_mean = mean(tajd4))
# summary_data_tajd4$tajd_type = 'fourfold SNPs'
summary_data_tajdI = summarise(group_by(genes_300min, dN_bin), dN_max = round(max(dN), 3), taj_mean = mean(tajd_indel, na.rm=TRUE))
summary_data_tajdI$tajd_type = 'INDELs'

tajd_plot = ggplot(rbind(summary_data_tajd0, summary_data_tajdI), aes(x=as.factor(dN_max), y=taj_mean, colour=tajd_type)) +
  geom_point(stat='identity', size = 2) + #ylim(0, 0.02) +
  theme_bw(base_size = 11) +
  scale_colour_manual(values=cbPalette) +
  ylab("Tajima's D mean") + xlab(expression('d'[N])) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle=90),
        legend.position = "none")

tajd_plot = arrangeGrob(tajd_plot, top = textGrob("(B)", x = unit(0, "npc"),
            y = unit(0, "npc"), just=c("left","top"),
            gp=gpar(col="black", fontsize=11)))

combined = grid.arrange(pi_plot, tajd_plot, ncol=2)

ggsave('dmel.pi_tajd_dn.pdf', width=6, height=3, plot=combined)

# spearmans rank cor
# cor.test(summary_data_tajd0$taj_mean, summary_data_tajdI$taj_mean, method='spearman')