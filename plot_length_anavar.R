library(ggplot2)
library(gridExtra)

len_data = read.csv('dmel_sel_v_neu_anavar_1class_equal_t_lengths.csv')

sel_data = subset(len_data, sel_type=='sel')

# theta plot
t_plot = ggplot(sel_data, aes(x=length, y=theta, colour=var_type)) +
    geom_point(stat='identity') +
    theme_bw() +
    theme(legend.position='none') +
    ylab(expression(theta)) + xlab(' \n\n\n')

# gamma plot
g_plot = ggplot(sel_data, aes(x=length, y=gamma, colour=var_type)) +
    geom_point(stat='identity') +
    theme_bw() + xlab('length\n\n\n')  +
    ylab(expression(gamma))

# n plot
n_data = read.delim('length_analysis_segsite_numbers.txt')

n_plot = ggplot(n_data, aes(x=length, y=n_segsites, colour=var_type)) +
    geom_point(stat='identity') +
    theme_bw() + xlab('') +
    theme(legend.title=element_blank(), legend.position=c(0.15, 0.80),
    axis.text.x=element_text(angle=45, hjust=1))

png('length_anavar.png', height=3, width=9, units='in', res=320)

grid.arrange(t_plot, g_plot + theme(legend.title=element_blank(), legend.position='none'), n_plot, nrow=1)

dev.off()

pdf('length_anavar_gamma.pdf', width=3, height=3)

g_plot + theme(legend.title=element_blank(), legend.position=c(0.15, 0.82),
legend.box.background = element_rect(), legend.margin=margin(.5, 5, 1, 1)) + ylim(min(sel_data$gamma), 0)


dev.off()