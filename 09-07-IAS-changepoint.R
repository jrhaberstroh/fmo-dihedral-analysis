library(data.table)
ALL.DF <- read.table("/home/jhaberstroh/data/SCRATCH/2016-09-all_t_dataframe.csv",header=TRUE, sep=',')

cols.cdc <- grepl("RESID", names(ALL.DF)) | grepl("bcl", names(ALL.DF)) | grepl("solvent", names(ALL.DF)) | grepl("ion", names(ALL.DF))
cols.label <- names(ALL.DF) == 'time..ns' | names(ALL.DF) == 'dataset'
cols.dih <- grepl("^d[[:digit:]]", names(ALL.DF))
cols.angle <- grepl("^a[[:digit:]]", names(ALL.DF))

CDC.DF <- ALL.DF[, cols.cdc]

library(changepoint)
library(dplyr)
library(ggplot2)
# In [2581]: cdc_df.columns[cdc_std > 15]
# Out[2581]: Index([u'RESID298', u'RESID305', u'RESID354'], dtype='object')
cdcstd.max <- names(CDC.DF) %in% c('RESID298', 'RESID305', 'RESID354')
t.170 <- filter(ALL.DF, dataset == '170nsX10ps') %>% select(time..ns)

################################################################################
## RESIDUE 298
################################################################################
r298.170 <- filter(ALL.DF, dataset == '170nsX10ps') %>% select(RESID298)
qplot(unlist(t.170), unlist(r298.170))

## Perform the changepoint analysis with only mean data
r298.170.cpt <- cpt.mean(unlist(r298.170), Q= 20, method="BinSeg",penalty='MBIC', minseglen=400)
r298.170.cptime <- t.170[r298.170.cpt@cpts, 1]
r298.170.cpt.segments <- data.frame(x1 = c(0.0, head(r298.170.cptime, -1)), 
                                    x2= r298.170.cptime, 
                                    mean=r298.170.cpt@param.est)
qplot(unlist(t.170), unlist(r298.170), alpha=.3) + 
        geom_segment(aes(x = x1, y = mean, xend = x2, yend=mean, size=2, alpha=1), data= r298.170.cpt.segments)
## Conclusion: Using only mean gives too many changeponits! 

## Perform the changepoint analysis with mean and variance data
r298.170.cpt <- cpt.meanvar(unlist(r298.170), Q= 17, method="BinSeg",penalty='MBIC', minseglen=2)
r298.170.cptime <- t.170[r298.170.cpt@cpts, 1]
r298.170.cpt.segments <- data.frame(x1 = c(0.0, head(r298.170.cptime, -1)), 
                                    x2= r298.170.cptime, 
                                    mean=r298.170.cpt@param.est$mean)
ggplot(filter(ALL.DF, dataset == '170nsX10ps'), aes(x=time..ns, y=RESID298)) + geom_point(alpha=.5, color='pink') +
        geom_segment(aes(x = x1, y = mean, xend = x2, yend=mean), size=5, alpha=.8, data= r298.170.cpt.segments) +
        ggsave("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-09/RESID298_cpa.png", height=4, width=12)

################################################################################
## OTHER RESIDUES CPA
################################################################################
## 305
r305.170 <- filter(ALL.DF, dataset == '170nsX10ps') %>% select(RESID305)
r305.170.cpt <- cpt.meanvar(unlist(r305.170), Q= 17, method="BinSeg",penalty='MBIC', minseglen=2)
r305.170.cptime <- t.170[r305.170.cpt@cpts, 1]
r305.170.cpt.segments <- data.frame(x1 = c(0.0, head(r305.170.cptime, -1)),
                                    x2= r305.170.cptime,
                                    mean=r305.170.cpt@param.est$mean)
ggplot(filter(ALL.DF, dataset == '170nsX10ps'), aes(x=time..ns, y=RESID305)) + geom_point(alpha=.5, color='pink') +
        geom_segment(aes(x = x1, y = mean, xend = x2, yend=mean), color='black', size=5, alpha=.8, data= r305.170.cpt.segments) +
        ggsave("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-09/RESID305_cpa.png", height=4, width=12)

## 354
r354.170 <- filter(ALL.DF, dataset == '170nsX10ps') %>% select(RESID354)
r354.170.cpt <- cpt.meanvar(unlist(r354.170), Q= 17, method="BinSeg",penalty='MBIC', minseglen=2)
r354.170.cptime <- t.170[r354.170.cpt@cpts, 1]
r354.170.cpt.segments <- data.frame(x1 = c(0.0, head(r354.170.cptime, -1)),
                                    x2= r354.170.cptime,
                                    mean=r354.170.cpt@param.est$mean)
#qplot(unlist(t.170), unlist(r354.170), alpha=.3, color='pink') +
ggplot(filter(ALL.DF, dataset == '170nsX10ps'), aes(x=time..ns, y=RESID354)) + geom_point(alpha=.5, color='pink') +
        geom_segment(aes(x = x1, y = mean, xend = x2, yend=mean), size=5, color='black',alpha=.8, data= r354.170.cpt.segments) +
        ggsave("/home/jhaberstroh/Dropbox/GraduateSchool/subgroup/2016-09-09/RESID354_cpa.png", height=4, width=12)








