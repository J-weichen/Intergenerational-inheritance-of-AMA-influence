#loading package
rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(ggpubr)
library(methylKit)
library(scatterplot3d)
library(pcaMethods)
library(scales)
library(ggsci)
library(grid)
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

grid.newpage(); #清空画板，开始画新图

AMA_father_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/AMA_file/Father")
YOUNG_father_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/YOUNG_file/Father")
AMA_Mother_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/AMA_file/Mother")
YOUNG_Mother_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/YOUNG_file/Mother")
AMA_Kids_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/AMA_file/Kids")
YOUNG_Kids_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/YOUNG_file/Kids")

AMA_father_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/AMA_file/Father/",AMA_father_files))
YOUNG_father_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/YOUNG_file/Father/",YOUNG_father_files))
AMA_Mother_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/AMA_file/Mother/",AMA_Mother_files))
YOUNG_Mother_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/YOUNG_file/Mother/",YOUNG_Mother_files))
AMA_Kids_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/AMA_file/Kids/",AMA_Kids_files))
YOUNG_Kids_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/YOUNG_file/Kids/",YOUNG_Kids_files))

AMA_father_names<- lapply(strsplit(AMA_father_files,"_"), function(x) x[1])
YOUNG_father_names<- lapply(strsplit(YOUNG_father_files,"_"), function(x) x[1])
AMA_Mother_names<- lapply(strsplit(AMA_Mother_files,"_"), function(x) x[1])
YOUNG_Mother_names<- lapply(strsplit(YOUNG_Mother_files,"_"), function(x) x[1])
AMA_Kids_names<- lapply(strsplit(AMA_Kids_files,"_"), function(x) x[1])
YOUNG_Kids_names<- lapply(strsplit(YOUNG_Kids_files,"_"), function(x) x[1])
#covarage 6
myobj=methRead(c(AMA_Kids_files2,YOUNG_Kids_files2,AMA_Mother_files2,YOUNG_Mother_files2,AMA_father_files2,YOUNG_father_files2),
               sample.id=c(AMA_Kids_names,YOUNG_Kids_names,AMA_Mother_names,YOUNG_Mother_names,AMA_father_names,YOUNG_father_names),
               assembly="hg38",pipeline="bismarkCoverage",#default:'amp'
               treatment=c(rep(0,length(AMA_Kids_names)),rep(1,length(YOUNG_Kids_names)),
                           rep(2,length(AMA_Mother_names)),rep(3,length(YOUNG_Mother_names)),
                           rep(4,length(AMA_father_names)),rep(5,length(YOUNG_father_names))),
               context="CpG", resolution = "base", mincov = 6)
dim(as.data.frame(myobj[1]))
saveRDS(myobj, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_myobj_6X_CpG.rds")


#进一步过滤 设置最低覆盖度
myobj<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_myobj_6X_CpG.rds")
depth<-6
filtered.myobj=filterByCoverage(myobj,lo.count=depth,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)

#基因组划分区域进行计算
tiles_100=tileMethylCounts(filtered.myobj,win.size=100,step.size=100,cov.bases =1)
meth_100bp=methylKit::unite(tiles_100, destrand=FALSE,min.per.group=0L,mc.cores = 40)

data_raw_CpG<-as.data.frame(percMethylation(meth_100bp,rowids=T))
meth_100bp_2<-data_raw_CpG[grep("chrUn_*|*_alt|*random|chrM|chrY",rownames(data_raw_CpG),invert=TRUE),]
head(meth_100bp_2);tail(meth_100bp_2);dim(meth_100bp_2)

colMeans(meth_100bp_2,na.rm=TRUE)
# 11-3-3F   12-3-4A   13-3-1G      14-3   15-3-2C        1Q    2-3-3C       3-3       4-3       5-3   10-3-1E   16-3-4D 
# 76.82441  78.27750  77.78294  77.38372  78.04758  76.14811  77.98943  77.25739  77.32195  77.99291  77.69117  78.20885 
# 17-3-2F      18-3   19-3-4E    20Q-4H       6-3       7-3       8-3      E16C      11-1   12-1-3G   13-1-1F      14-1 
# 77.28143  77.13443  79.65598  78.35181  77.03262  77.23707  78.06532  77.22179  76.44639  77.08684  77.41927  76.31513 
# 15-1        1M       2-1        3M       4-1       5-1 10-1-in-9      16-1   17-1-2D      18-1   19-1-2G    20M-4F 
# 73.87670  75.76546  76.52576  75.66734  76.07267  76.69305  75.95893  75.85804  77.24156  76.37562  78.06511  76.47197 
# 6-1       7-1       8-1      E16M      11-2   12-2-3H   13-2-1H      14-2      15-2        1F       2-2        3F 
# 75.76011  76.29330  76.36769  75.36024  75.93244  76.16837  76.70772  75.83063  74.85613  76.52999  75.11345  76.34693 
# 4-2       5-2   10-2-1D   16-2-4C   17-2-2E      18-2   19-2-2H    20F-4G       6-2        7F       8-2      E16F 
#74.95357  76.09993  77.77129  77.03603  76.10706  74.60186  78.66336  73.41562  75.56652  76.92427  76.26255  76.47306 

colMeans(data_raw_CpG,na.rm=TRUE)
#  11-3-3F   12-3-4A   13-3-1G      14-3   15-3-2C        1Q    2-3-3C       3-3       4-3       5-3   10-3-1E   16-3-4D 
#  76.83196  78.27681  77.78702  77.38314  78.05171  76.15535  77.98838  77.26746  77.33141  78.00191  77.69842  78.20803 
#  17-3-2F      18-3   19-3-4E    20Q-4H       6-3       7-3       8-3      E16C      11-1   12-1-3G   13-1-1F      14-1 
#  77.28077  77.13376  79.65968  78.36014  77.04046  77.23658  78.07505  77.23314  76.44594  77.08616  77.41813  76.31444 
#  15-1        1M       2-1        3M       4-1       5-1 10-1-in-9      16-1   17-1-2D      18-1   19-1-2G    20M-4F 
#  73.87573  75.76499  76.52506  75.66642  76.07190  76.69241  75.95847  75.85757  77.24076  76.37498  78.06461  76.47141 
#  6-1       7-1       8-1      E16M      11-2   12-2-3H   13-2-1H      14-2      15-2        1F       2-2        3F 
#  75.75959  76.29271  76.36726  75.35994  75.94172  76.17666  76.71439  75.84208  74.86689  76.53974  75.12326  76.35591 
#  4-2       5-2   10-2-1D   16-2-4C   17-2-2E      18-2   19-2-2H    20F-4G       6-2        7F       8-2      E16F 
#  74.96276  76.10822  77.78146  77.04397  76.11524  74.61442  78.66889  73.42290  75.57398  76.93347  76.27120  76.48429

##基因组划分区域进行计算
tiles_200=tileMethylCounts(filtered.myobj,win.size=200,step.size=200,cov.bases =1)
meth_200bp=methylKit::unite(tiles_200, destrand=FALSE,min.per.group=0L,mc.cores = 40)
data_raw_CpG<-as.data.frame(percMethylation(meth_200bp,rowids=T))
meth_200bp<-data_raw_CpG[grep("chrUn_*|*_alt|*random|chrM|chrY",rownames(data_raw_CpG),invert=TRUE),]
head(meth_200bp);tail(meth_200bp);dim(meth_200bp)
colMeans(meth_200bp,na.rm=TRUE)
# 11-3-3F   12-3-4A   13-3-1G      14-3   15-3-2C        1Q    2-3-3C       3-3       4-3       5-3   10-3-1E   16-3-4D 
# 78.60202  80.18197  79.24642  79.33206  79.16797  77.97883  79.73117  79.08746  79.09648  79.93883  79.57730  80.09511 
# 17-3-2F      18-3   19-3-4E    20Q-4H       6-3       7-3       8-3      E16C      11-1   12-1-3G   13-1-1F      14-1 
# 79.19577  79.10869  80.87017  80.05619  78.99440  79.20222  79.95259  79.27195  78.08192  78.89438  78.46706  78.23573 
# 15-1        1M       2-1        3M       4-1       5-1 10-1-in-9      16-1   17-1-2D      18-1   19-1-2G    20M-4F 
# 75.53217  77.57444  78.18659  77.31798  77.98921  78.51825  77.91922  77.82222  78.99743  78.24994  79.18940  78.33551 
# 6-1       7-1       8-1      E16M      11-2   12-2-3H   13-2-1H      14-2      15-2        1F       2-2        3F 
# 77.59550  78.23091  78.24578  77.38483  77.67830  77.96595  77.88016  77.76736  76.66477  78.31137  77.01091  78.12036 
# 4-2       5-2   10-2-1D   16-2-4C   17-2-2E      18-2   19-2-2H    20F-4G       6-2        7F       8-2      E16F 
# 76.79974  77.81721  79.61060  78.82770  77.80550  76.59807  79.73559  75.15584  77.39744  78.79976  78.10216  78.44325

colMeans(data_raw_CpG,na.rm=TRUE)
# 11-3-3F   12-3-4A   13-3-1G      14-3   15-3-2C        1Q    2-3-3C       3-3       4-3       5-3   10-3-1E   16-3-4D 
# 78.60748  80.18101  79.24893  79.33124  79.17144  77.98490  79.72964  79.09603  79.10420  79.94603  79.58256  80.09398 
# 17-3-2F      18-3   19-3-4E    20Q-4H       6-3       7-3       8-3      E16C      11-1   12-1-3G   13-1-1F      14-1 
# 79.19476  79.10778  80.87276  80.06260  79.00023  79.20155  79.96045  79.28194  78.08126  78.89350  78.46562  78.23480 
# 15-1        1M       2-1        3M       4-1       5-1 10-1-in-9      16-1   17-1-2D      18-1   19-1-2G    20M-4F 
# 75.53089  77.57374  78.18555  77.31675  77.98821  78.51738  77.91857  77.82171  78.99637  78.24922  79.18871  78.33474 
# 6-1       7-1       8-1      E16M      11-2   12-2-3H   13-2-1H      14-2      15-2        1F       2-2        3F 
# 77.59482  78.23009  78.24515  77.38431  77.68597  77.97196  77.88572  77.77677  76.67352  78.31902  77.01885  78.12730 
# 4-2       5-2   10-2-1D   16-2-4C   17-2-2E      18-2   19-2-2H    20F-4G       6-2        7F       8-2      E16F 
# 76.80739  77.82296  79.61820  78.83345  77.81180  76.60905  79.74004  75.16152  77.40296  78.80697  78.10897  78.45261 

#基因组划分区域进行计算
tiles_300=tileMethylCounts(filtered.myobj,win.size=300,step.size=300,cov.bases =1)
meth_300bp=methylKit::unite(tiles_300, destrand=FALSE,min.per.group=0L,mc.cores = 40)
data_raw_CpG<-as.data.frame(percMethylation(meth_300bp,rowids=T))
meth_300bp<-data_raw_CpG[grep("chrUn_*|*_alt|*random|chrM|chrY",rownames(data_raw_CpG),invert=TRUE),]
head(meth_300bp);tail(meth_300bp);dim(meth_300bp)

colMeans(meth_300bp,na.rm=TRUE)
# 11-3-3F   12-3-4A   13-3-1G      14-3   15-3-2C        1Q    2-3-3C       3-3       4-3       5-3   10-3-1E   16-3-4D 
# 79.51158  81.16635  80.07382  80.35720  79.78182  78.94871  80.65334  80.07467  80.05590  80.96554  80.53864  81.07243 
# 17-3-2F      18-3   19-3-4E    20Q-4H       6-3       7-3       8-3      E16C      11-1   12-1-3G   13-1-1F      14-1 
# 80.19912  80.11347  81.53647  80.97422  80.01522  80.22248  80.95962  80.32640  78.97851  79.82798  79.03453  79.24298 
# 15-1        1M       2-1        3M       4-1       5-1 10-1-in-9      16-1   17-1-2D      18-1   19-1-2G    20M-4F 
# 76.35864  78.53192  79.09577  78.14945  78.98315  79.48197  78.94856  78.83591  79.89878  79.25043  79.79756  79.28610 
# 6-1       7-1       8-1      E16M      11-2   12-2-3H   13-2-1H      14-2      15-2        1F       2-2        3F 
# 78.55871  79.25407  79.21976  78.42972  78.62476  78.89573  78.53636  78.78346  77.59073  79.22848  78.00200  79.04164 
# 4-2       5-2   10-2-1D   16-2-4C   17-2-2E      18-2   19-2-2H    20F-4G       6-2        7F       8-2      E16F 
# 77.75318  78.71966  80.55588  79.77768  78.71475  77.62285  80.30782  76.05469  78.34488  79.76621  79.05178  79.46168 
colMeans(data_raw_CpG,na.rm=TRUE)
# 11-3-3F   12-3-4A   13-3-1G      14-3   15-3-2C        1Q    2-3-3C       3-3       4-3       5-3   10-3-1E   16-3-4D 
# 79.51614  81.16539  80.07566  80.35639  79.78519  78.95396  80.65177  80.08269  80.06284  80.97174  80.54282  81.07118 
# 17-3-2F      18-3   19-3-4E    20Q-4H       6-3       7-3       8-3      E16C      11-1   12-1-3G   13-1-1F      14-1 
# 80.19814  80.11249  81.53867  80.97984  80.01952  80.22181  80.96682  80.33553  78.97782  79.82707  79.03312  79.24188 
# 15-1        1M       2-1        3M       4-1       5-1 10-1-in-9      16-1   17-1-2D      18-1   19-1-2G    20M-4F 
# 76.35722  78.53117  79.09473  78.14821  78.98199  79.48098  78.94787  78.83515  79.89773  79.24948  79.79685  79.28526 
# 6-1       7-1       8-1      E16M      11-2   12-2-3H   13-2-1H      14-2      15-2        1F       2-2        3F 
# 78.55808  79.25334  79.21906  78.42924  78.63155  78.90044  78.54147  78.79206  77.59856  79.23533  78.00919  79.04741 
# 4-2       5-2   10-2-1D   16-2-4C   17-2-2E      18-2   19-2-2H    20F-4G       6-2        7F       8-2      E16F 
# 77.75985  78.72452  80.56266  79.78213  78.72008  77.63311  80.31214  76.05973  78.34946  79.77263  79.05717  79.47022 

#基因组划分区域进行计算
tiles_500=tileMethylCounts(filtered.myobj,win.size=500,step.size=500,cov.bases =1)
meth_500bp=methylKit::unite(tiles_500, destrand=FALSE,min.per.group=0L,mc.cores = 40)
data_raw_CpG<-as.data.frame(percMethylation(meth_500bp,rowids=T))
meth_500bp<-data_raw_CpG[grep("chrUn_*|*_alt|*random|chrM|chrY",rownames(data_raw_CpG),invert=TRUE),]
head(meth_500bp);tail(meth_500bp);dim(meth_500bp)

colMeans(meth_500bp,na.rm=TRUE)
# 11-3-3F   12-3-4A   13-3-1G      14-3   15-3-2C        1Q    2-3-3C       3-3       4-3       5-3   10-3-1E   16-3-4D 
# 80.41917  82.17659  80.97716  81.41329  80.46389  79.94558  81.63113  81.11409  81.07411  82.03770  81.53142  82.07932 
# 17-3-2F      18-3   19-3-4E    20Q-4H       6-3       7-3       8-3      E16C      11-1   12-1-3G   13-1-1F      14-1 
# 81.22874  81.18266  82.29789  81.92888  81.06498  81.30665  82.00420  81.48780  79.96057  80.77119  79.68102  80.30936 
# 15-1        1M       2-1        3M       4-1       5-1 10-1-in-9      16-1   17-1-2D      18-1   19-1-2G    20M-4F 
# 77.20001  79.51729  80.07940  78.99328  80.02017  80.49799  80.03122  79.92983  80.84318  80.31219  80.48414  80.27528 
# 6-1       7-1       8-1      E16M      11-2   12-2-3H   13-2-1H      14-2      15-2        1F       2-2        3F 
# 79.56629  80.33053  80.25516  79.59070  79.63736  79.85087  79.27713  79.84474  78.55758  80.18191  79.04009  79.99258 
# 4-2       5-2   10-2-1D   16-2-4C   17-2-2E      18-2   19-2-2H    20F-4G       6-2        7F       8-2      E16F 
# 78.75177  79.65829  81.55385  80.73529  79.67306  78.75771  80.96579  76.96402  79.33052  80.75675  80.07589  80.56415 

colMeans(data_raw_CpG,na.rm=TRUE)
# 11-3-3F   12-3-4A   13-3-1G      14-3   15-3-2C        1Q    2-3-3C       3-3       4-3       5-3   10-3-1E   16-3-4D 
# 80.42272  82.17525  80.97772  81.41213  80.46613  79.95036  81.62899  81.12105  81.07964  82.04273  81.53412  82.07772 
# 17-3-2F      18-3   19-3-4E    20Q-4H       6-3       7-3       8-3      E16C      11-1   12-1-3G   13-1-1F      14-1 
# 81.22725  81.18122  82.29920  81.93311  81.06848  81.30552  82.01005  81.49582  79.95935  80.77004  79.67927  80.30785 
# 15-1        1M       2-1        3M       4-1       5-1 10-1-in-9      16-1   17-1-2D      18-1   19-1-2G    20M-4F 
# 77.19811  79.51629  80.07793  78.99170  80.01871  80.49651  80.03019  79.92881  80.84169  80.31086  80.48281  80.27412 
# 6-1       7-1       8-1      E16M      11-2   12-2-3H   13-2-1H      14-2      15-2        1F       2-2        3F 
# 79.56499  80.32940  80.25418  79.58988  79.64268  79.85472  79.28140  79.85254  78.56481  80.18788  79.04677  79.99731 
# 4-2       5-2   10-2-1D   16-2-4C   17-2-2E      18-2   19-2-2H    20F-4G       6-2        7F       8-2      E16F 
# 78.75789  79.66185  81.55989  80.73912  79.67738  78.76620  80.96983  76.96813  79.33403  80.76212  80.08064  80.57141

tile_list<-c(list(tiles_100),list(tiles_200),list(tiles_300),list(tiles_500))
meth_list<-c(list(meth_100bp),list(meth_200bp),list(meth_300bp),list(meth_500bp))
saveRDS(tile_list, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_6X_tile_list.rds")
saveRDS(meth_list, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_6X_meth_list.rds")

##final_determination corverage and bins for downstream calculation
#rm(list = ls())
#library(methylKit)
myobj<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_myobj_6X_CpG.rds")
depth<-6
filtered.myobj=filterByCoverage(myobj,lo.count=depth,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
#基因组划分区域进行计算
tiles_200=tileMethylCounts(filtered.myobj,win.size=200,step.size=200,cov.bases =3)
meth_200bp=methylKit::unite(tiles_200, destrand=FALSE,min.per.group=0L,mc.cores = 48)
data_raw_CpG<-as.data.frame(percMethylation(meth_200bp,rowids=T))
meth_200bp_2<-data_raw_CpG[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",rownames(data_raw_CpG),invert=TRUE),]
head(meth_200bp_2);tail(meth_200bp_2);dim(meth_200bp_2)
colMeans(meth_200bp_2,na.rm=TRUE)
# 11-3-3F   12-3-4A   13-3-1G      14-3   15-3-2C        1Q    2-3-3C       3-3       4-3       5-3   10-3-1E   16-3-4D 
# 76.30587  77.42403  76.63477  76.53226  76.98775  75.02250  77.00548  76.19401  76.39127  76.88812  76.98639  77.44106 
# 17-3-2F      18-3   19-3-4E    20Q-4H       6-3       7-3       8-3      E16C      11-1   12-1-3G   13-1-1F      14-1 
# 76.42777  76.29513  78.35406  77.39758  76.18036  76.20507  77.08103  76.12252  75.54060  76.48684  76.63453  75.48266 
# 15-1        1M       2-1        3M       4-1       5-1 10-1-in-9      16-1   17-1-2D      18-1   19-1-2G    20M-4F 
# 73.16294  75.08505  75.69402  75.40831  75.45990  75.92407  75.18212  75.08057  76.48092  75.42804  77.07612  75.89212 
# 6-1       7-1       8-1      E16M      11-2   12-2-3H   13-2-1H      14-2      15-2        1F       2-2        3F 
# 75.02502  75.38179  75.69316  74.43111  75.14638  75.66188  75.90849  75.13554  74.43584  75.99729  74.50678  75.81351 
# 4-2       5-2   10-2-1D   16-2-4C   17-2-2E      18-2   19-2-2H    20F-4G       6-2        7F       8-2      E16F 
# 74.34319  75.45553  76.91631  76.36086  75.35751  73.79001  77.57162  73.10204  74.86593  76.26711  75.45062  75.49482

colMeans(data_raw_CpG,na.rm=TRUE)
# 11-3-3F   12-3-4A   13-3-1G      14-3   15-3-2C        1Q    2-3-3C       3-3       4-3       5-3   10-3-1E   16-3-4D 
# 76.36108  77.46737  76.69972  76.59039  77.04092  75.03506  77.04177  76.25661  76.45573  76.94443  77.03486  77.48569 
# 17-3-2F      18-3   19-3-4E    20Q-4H       6-3       7-3       8-3      E16C      11-1   12-1-3G   13-1-1F      14-1 
# 76.48109  76.32861  78.40187  77.45940  76.23432  76.25657  77.14389  76.18734  75.59409  76.54286  76.67177  75.54182 
# 15-1        1M       2-1        3M       4-1       5-1 10-1-in-9      16-1   17-1-2D      18-1   19-1-2G    20M-4F 
# 73.22399  75.13382  75.76610  75.45453  75.51665  75.98562  75.24149  75.14312  76.53396  75.49159  77.08956  75.96357 
# 6-1       7-1       8-1      E16M      11-2   12-2-3H   13-2-1H      14-2      15-2        1F       2-2        3F 
# 75.10314  75.44475  75.74572  74.51508  75.21153  75.72390  75.97056  75.19274  74.49539  76.05488  74.56354  75.86550 
# 4-2       5-2   10-2-1D   16-2-4C   17-2-2E      18-2   19-2-2H    20F-4G       6-2        7F       8-2      E16F 
# 74.39830  75.51132  76.96586  76.41224  75.41568  73.83993  77.61754  73.12774  74.91224  76.31561  75.50717  75.54627

saveRDS(tiles_200, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_6X_tiles_200.rds")
saveRDS(meth_200bp_2, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_6X_meth_200bp_2.rds")
write.table(as.data.frame(meth_200bp_2), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_200bp_CpG_6X_3site_bismark_no_filter.txt",quote=F, row.names=F, col.names=T) 

#plot for view
#tiles_200_2 <- readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/tiles_200.rds")
meth_200bp_2<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_6X_meth_200bp_2.rds")

huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
head(huahua_meta)
#tiles_200_2[1:6,];meth_200bp_2[1:6,]
rownames(huahua_meta)<-huahua_meta$library_code

#count 200bp bin in each samples
#filtering and arrange data
evaluation.table<-meth_200bp_2
evaluation.table[!(is.na(evaluation.table))]<- 1
evaluation.table[is.na(evaluation.table)]<- 0
bin200_detected_number<-colSums(evaluation.table)
bin200_detected_number<-data.frame(bin200_detected_number)

head(huahua_meta);dim(huahua_meta)
colnames(bin200_detected_number) <- "bin200_number"
rownames(bin200_detected_number);rownames(huahua_meta)
bin200_detected_number2<-merge(bin200_detected_number,huahua_meta,by=0)
head(bin200_detected_number2);dim(bin200_detected_number2)
bin200_detected_number2$group2<-paste(bin200_detected_number2$Age_group,bin200_detected_number2$Sample.types,sep="_")

bin200_detected_number2$Age_group <-factor(bin200_detected_number2$Age_group,level=c("YOUNG","AMA"))
bin200_detected_number2$Sample.types <-factor(bin200_detected_number2$Sample.types,level=c("kids","mother","father"))
bin200_detected_number2$group2 <-factor(bin200_detected_number2$group2,level=c("YOUNG_kids","AMA_kids","YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))
#bin200_detected_number2$Gender_RRBS<-factor(bin200_detected_number2$Gender_RRBS,levels = c("Female","Male"))
head(bin200_detected_number2)
bin200_detected_number2$analysis_name<-as.character(bin200_detected_number2$analysis_name)
bin200_detected_number2<-bin200_detected_number2[order(bin200_detected_number2$group2,bin200_detected_number2$bin200_number,decreasing = T),]
bin200_detected_number2$analysis_name<-factor(bin200_detected_number2$analysis_name,levels=as.character(bin200_detected_number2$analysis_name))

#plot
plot_bin200_detected_number0<-ggplot(data=bin200_detected_number2, mapping=aes(x=analysis_name,y=bin200_number,fill=group2))+geom_bar(stat="identity",width=0.8)
plot_bin200_detected_number1 <-plot_bin200_detected_number0+scale_fill_manual(values=ppCor[c(7:8,2:1,3,5)])+
  geom_text(stat="identity",aes(label=bin200_number), color="black", size=1,position=position_stack(1.01))+
  theme_classic()+labs(x="",y="Number of 200bp bin",title="Number of bin detected")+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin200_detected_number1
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.number_of_200bp_bins_detected_in_all_samples_by_6X_CpG_3L.pdf",plot_bin200_detected_number1,width = 16, height = 8)

#plot
plot_bin200_detected_number0<-ggplot(data=bin200_detected_number2, mapping=aes(x=analysis_name,y=bin200_number,fill=type_Batch))+geom_bar(stat="identity",width=0.8)
plot_bin200_detected_number2 <-plot_bin200_detected_number0+scale_fill_manual(values=ppCor[c(8,3)])+
  geom_text(stat="identity",aes(label=bin200_number), color="black", size=1,position=position_stack(1.01))+
  theme_classic()+labs(x="",y="Number of 200bp bin",title="Number of bin detected")+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin200_detected_number2
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.number_of_200bp_bins_6X_CpG_3L_detected_in_per_group_all_sample_no_filter_two.pdf",plot_bin200_detected_number2,width = 16, height =15)

#plot
plot_bin200_detected_number0<-ggplot(data=bin200_detected_number2, mapping=aes(x=analysis_name,y=bin200_number,fill=seq_Batch))+geom_bar(stat="identity",width=0.8)
plot_bin200_detected_number3 <-plot_bin200_detected_number0+scale_fill_manual(values=ppCor)+
  geom_text(stat="identity",aes(label=bin200_number), color="black", size=1,position=position_stack(1.01))+
  theme_classic()+labs(x="",y="Number of 200bp bin",title="Number of bin detected")+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin200_detected_number3
count_plot<-ggarrange(plot_bin200_detected_number1,plot_bin200_detected_number2,plot_bin200_detected_number3,labels = c("A", "B","C"),ncol = 1, nrow = 3)
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.number_of_200bp_bins_6X_CpG_3L_detected_in_per_group_all_sample_no_filter_three.pdf",count_plot,width = 16, height =15)

#plot_bin200_detected_number12<-ggplot(data=bin200_detected_number2, mapping=aes(x=Row.names,y=bin200_number,fill=Gender_RRBS))+geom_bar(stat="identity",width=0.8)
#plot_bin200_detected_number42 <-plot_bin200_detected_number12+scale_fill_manual(values=ppCor[c(10:9)])+
#  geom_text(stat="identity",aes(label=bin200_number), color="black", size=1.5,position=position_stack(1.05))+
#  theme_classic()+labs(x="",y="Number of genes",title="Number of gene detected")+guides(fill=guide_legend(ncol=1)) +
#  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
#        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
#        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
#plot_bin200_detected_number42
#count_plot<-ggarrange(plot_bin200_detected_number3,plot_bin200_detected_number42,labels = c("A", "B"),ncol = 1, nrow = 2)
#ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.number_of_200bp_bins_detected_in_all_samples_by_6X_CpG_3L_no_filtersamples_library.pdf",count_plot,width = 300, height =150, units = "mm")

##使用平均值进行绘图
data_region1<-as.data.frame(colMeans(meth_200bp_2,na.rm = T))
rownames(data_region1) 
length(na.omit(rownames(data_region1)));length(na.omit(rownames(huahua_meta)))
all(rownames(data_region1) %in% rownames(huahua_meta))

colnames(data_region1)<-"meth_mean"
data_region<-merge(data_region1,huahua_meta,by=0)
dim(data_region)#60 10
colnames(data_region)
data_region$Age_group <-factor(data_region$Age_group,level=c("YOUNG","AMA"))
data_region$Sample.types <-factor(data_region$Sample.types,level=c("kids","mother","father"))

head(data_region);dim(data_region)
#compare_means(meth_mean~Age_group, data=data_region,method = "t.test",group_by="group_age")
compare_means(meth_mean~Age_group,group.by = "Sample.types", data=data_region,method = "t.test")
#Sample.types .y.       group1 group2     p p.adj p.format p.signif method
#1 mother       meth_mean YOUNG  AMA    0.840  1    0.84     ns       T-test
#2 father       meth_mean YOUNG  AMA    0.570  1    0.57     ns       T-test
#3 kids         meth_mean YOUNG  AMA    0.331  0.99 0.33     ns       T-test
compare_means(meth_mean~Age_group,group.by = "Sample.types", data=data_region,method = "wilcox.test")
#Sample.types .y.       group1 group2     p p.adj p.format p.signif method
#1 mother       meth_mean YOUNG  AMA    0.529     1 0.53     ns       Wilcoxon
#2 father       meth_mean YOUNG  AMA    0.529     1 0.53     ns       Wilcoxon
#3 kids         meth_mean YOUNG  AMA    0.739     1 0.74     ns       Wilcoxon

global_meth_plot1 <- ggboxplot(data_region, x = "Age_group", y = "meth_mean",facet.by = "Sample.types",
                               color = "Age_group", palette = ppCor[2:1],
                               add = "jitter")+xlab("") +ylab("DNA methylation level(CpG 5mC%)")+ ylim(c(60,90))
#  Add p-value
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.AMA_group_global_comparison_mean_all_site_by_6X_CpG_200bp_3L_nofilter.pdf")
global_meth_plot1+stat_compare_means(method = "t.test",comparisons= list(c("YOUNG","AMA")))
global_meth_plot1+stat_compare_means(method = "wilcox.test",comparisons= list(c("YOUNG","AMA")))
#p+stat_compare_means(aes(label=..p.signif..), label.x = 1.5, label.y = 1.1)
global_meth_plot1+stat_compare_means(method = "t.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 85)
global_meth_plot1+stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 85)
dev.off()

#grid.newpage(); #清空画板，开始画新图

#plot only for kids
global_meth_plot <- ggboxplot(data_region[which(data_region$Sample.types == "kids"),], x = "Age_group", y = "meth_mean",
                               color = "Age_group", palette = ppCor[2:1],
                               add = "jitter")+xlab("") +ylab("DNA methylation level(CpG 5mC%)")+ ylim(c(60,90))+
  stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 85)
global_meth_plot
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.kids_AMA_group_global_comparison_mean_all_site_by_6X_CpG_200bp_3L.pdf",global_meth_plot,width = 6, height = 6)
#plot only for parents
data_region2<-data_region[which(data_region$Sample.types != "kids"),]
data_region2$group2<-paste(data_region2$Age_group,data_region2$Sample.types,sep="_")
data_region2$group2 <-factor(data_region2$group2,level=c("YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))
my_comparisons <- list(c("YOUNG_mother","AMA_mother"),c("YOUNG_father","AMA_father"))
global_meth_plot <- ggboxplot(data_region2, x = "group2", y = "meth_mean",
                              color = "group2", palette = ppCor[c(2:1,3,5)],
                              add = "jitter")+xlab("") +ylab("DNA methylation level(CpG 5mC%)")+ ylim(c(60,90))+
  stat_compare_means(method = "kruskal.test",aes(label=paste0(..method..,": "," p = ",..p.format..)), label.x = 1,label.y = 85)+
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.format") +labs(title = "wilcox.test")
global_meth_plot
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.parents_AMA_group_global_comparison_mean_all_site_by_6X_CpG_200bp_3L.pdf",global_meth_plot,width = 8, height = 8)