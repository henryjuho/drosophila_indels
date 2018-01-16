library(ggplot2)

file_dir = 'indel_sfs/'

sfs_files = list.files(file_dir, pattern='*.txt')

sfs_data = data.frame()

for (x in sfs_files){
    file_info = strsplit(x, '_')[[1]]
    var_type = file_info[2]
    region = strsplit(file_info[3], '\\.')[[1]][1]

    sfs = read.delim(paste(file_dir, x, sep=''), sep=' ', header=F)

    sfs$V1 = NULL

    colnames(sfs) = c('count', 'freq')

    total_indel = sum(sfs$count)

    sfs$prop = sfs$count / total_indel

    sfs$var_type = var_type
    sfs$region = as.factor(region)

    sfs_data = rbind(sfs_data, sfs)

}

regional_data = subset(sfs_data, region!='ALL')

pdf('regional_indel_sfs.pdf', width=6, height=3)

ggplot(regional_data, aes(x=freq, y=prop, fill=region)) +
    geom_bar(stat='identity', position='dodge') +
    theme_bw() +
    facet_wrap(~var_type, nrow=1) +
    labs(x='Derived allele frequency', y='Proportion of variants') +
    theme(legend.title=element_blank(), legend.position=c(0.85, 0.6), axis.text=element_text(angle=45, hjust=1))

dev.off()

gwide_data = subset(sfs_data, region=='ALL')

pdf('gwide_indel_sfs.pdf', width=3, height=3)

ggplot(gwide_data, aes(x=freq, y=prop, fill=var_type)) +
    geom_bar(stat='identity', position='dodge') +
    theme_bw() +
    labs(x='Derived allele frequency', y='Proportion of variants') +
    theme(legend.title=element_blank(), legend.position=c(0.85, 0.7))

dev.off()