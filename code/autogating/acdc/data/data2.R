library(plyr)
ref_label2<-mapvalues(ref_label, from=0:20, to=c('unknown','naive b cells','memory b cells','classical monocytes','transitional monocytes',
                                                 'non-classical monocytes','basophils','early nks','late nks','pdc','mdc','cd8 naive','cd8 central memory',
                                                 'cd8 effector memory','cd8 terminal effector','cd4 naive','cd4 central memory','cd4 effector memory',
                                                 'cd4 terminal memory','cd4- mait/nkt','cd4-cd8-gd t cells'))

dt2.df<-data.frame(dt2.transform)
dt2.df<-data.frame(dt2.df,cell_type=ref_label2)
write.csv(dt2.df,"D:/research/Liza/cytof/CyTOF existing methods and challenges/programs/autogating/acdc/data/data2.csv")
