library(ggplot2)
library(gridExtra)

len_data = read.csv('dmel_sel_v_neu_anavar_2class_equal_t_lengths.csv')

sel_data = subset(len_data, sel_type=='sel' & gamma>=-40000)

# theta plot
t_plot = ggplot(sel_data, aes(x=length, y=theta, colour=var_type)) +
    geom_point(stat='identity') +
    theme_bw() +
    theme(legend.position='none') +
    ylab(expression(theta)) + xlab(' \n\n\n')

# gamma plot
g_plot = ggplot(sel_data, aes(x=length, y=gamma, colour=var_type)) +
    geom_point(stat='identity') +
    theme_bw() + xlab('length\n\n\n') +
    theme(legend.title=element_blank(), legend.position='none') +
    ylab(expression(gamma))

# n plot
n_data = read.delim('length_analysis_segsite_numbers.txt')

n_plot = ggplot(n_data, aes(x=length, y=n_segsites, colour=var_type)) +
    geom_point(stat='identity') +
    theme_bw() + xlab('') +
    theme(legend.title=element_blank(), legend.position=c(0.15, 0.80),
    axis.text.x=element_text(angle=45, hjust=1))

png('length_anavar_2class.png', height=3, width=9, units='in', res=320)

grid.arrange(t_plot, g_plot, n_plot, nrow=1)

dev.off()