# Title     : pol error plot
# Objective : demonstrate problem of pol error
# Created by: henryjuho
# Created on: 07/11/2017

library(ggplot2)
library(viridis)

cbPalette <- c("#E69F00", 'tomato 3', 'steel blue', "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

err <- 0.1
ins <- 10/(1:9)
del <- 20/(1:9)
ins_err <- numeric()
del_err <- numeric()
for(i in 1:9) {ins_err[i] <- (1 - err) * ins[i] + err * del[10 - i]; del_err[i] <- (1 - err) * del[i] + err * ins[10 - i]}

ins_data <- as.data.frame(cbind(1:9, ins))
colnames(ins_data) <- c('freq', 'count')
ins_data$err_type = 'No error'
ins_data$var_type = 'Insertions'
ins_err_data <- as.data.frame(cbind(1:9, ins_err))
colnames(ins_err_data) <- c('freq', 'count')
ins_err_data$err_type = 'With error'
ins_err_data$var_type = 'Insertions'
all_ins_data = rbind(ins_data, ins_err_data)

del_data <- as.data.frame(cbind(1:9, del))
colnames(del_data) <- c('freq', 'count')
del_data$err_type = 'No error'
del_data$var_type = 'Deletions'
del_err_data <- as.data.frame(cbind(1:9, del_err))
colnames(del_err_data) <- c('freq', 'count')
del_err_data$err_type = 'With error'
del_err_data$var_type = 'Deletions'
all_del_data = rbind(del_data, del_err_data)

plot_data = rbind(all_ins_data, all_del_data)
plot_data$var_type = factor(plot_data$var, levels=c('Insertions', 'Deletions'))

err_plot = ggplot(plot_data, aes(x=as.factor(freq), y=count, fill=err_type)) +
    geom_bar(stat='identity', position='dodge') +
    facet_wrap(~var_type, nrow=2) +
    theme_bw() + scale_fill_manual(values=cbPalette) +
    xlab('Derived allele frequency')  +
    ylab("Number of variants") +
    theme(legend.position=c(0.78, 0.9), legend.title=element_blank())

err_plot2 = ggplot(plot_data, aes(x=as.factor(freq), y=count, fill=err_type)) +
    geom_bar(stat='identity', position='dodge') +
    facet_wrap(~var_type, nrow=1) +
    theme_bw() + scale_fill_manual(values=viridis(3)) +
    xlab('Derived allele frequency')  +
    ylab("Number of variants") +
    theme(legend.position=c(0.9, 0.78), legend.title=element_blank())


pdf('polarisation_error_plot.pdf', width=3, height=6)

err_plot

dev.off()

pdf('pol_error_plot_thesis.pdf', width=6, height=3)

err_plot2

dev.off()
