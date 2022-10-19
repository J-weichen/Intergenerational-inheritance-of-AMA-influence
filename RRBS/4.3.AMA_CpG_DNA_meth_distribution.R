rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(grid)
library(ggpubr)
library(scales)

library(ggsci)
library(stringr)
library(RColorBrewer)
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

AMA_father_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/meth_around_gene/AMA_file/Father")
YOUNG_father_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/meth_around_gene/YOUNG_file/Father")
AMA_Mother_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/meth_around_gene/AMA_file/Mother")
YOUNG_Mother_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/meth_around_gene/YOUNG_file/Mother")
AMA_Kids_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/meth_around_gene/AMA_file/Kids")
YOUNG_Kids_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/meth_around_gene/YOUNG_file/Kids")

AMA_father_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/meth_around_gene/AMA_file/Father/",AMA_father_files))
YOUNG_father_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/meth_around_gene/YOUNG_file/Father/",YOUNG_father_files))
AMA_Mother_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/meth_around_gene/AMA_file/Mother/",AMA_Mother_files))
YOUNG_Mother_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/meth_around_gene/YOUNG_file/Mother/",YOUNG_Mother_files))
AMA_Kids_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/meth_around_gene/AMA_file/Kids/",AMA_Kids_files))
YOUNG_Kids_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/meth_around_gene/YOUNG_file/Kids/",YOUNG_Kids_files))
###################for kids
sample_list<-c(AMA_Kids_files2,YOUNG_Kids_files2)
sample_names<-c(paste0("AMA_",c(1:10)),paste0("Young_",c(1:10)))

cnane_name<-c()
data <- NULL
for(i in 1:20){
  # i<-5
  sample<-sample_list[[i]]
  Sample_Name<-sample_names[[i]]
  print(as.character(Sample_Name))
  
  data0 <-read.table(sample, header=F)
  head(data0)
  data0$V1<-NULL
  data1<-apply(data0,2,function(x){mean(x,na.rm=TRUE)*100})
  data <- rbind(data,data1)
  cnane_name<-c(cnane_name,as.character(Sample_Name))
}  
nrow(data);length(cnane_name)
rownames(data) <- cnane_name
data[1:6,1:5]
AMA_mean<-colMeans(data[1:10,],na.rm =TRUE)
Young_mean<-colMeans(data[11:20,],na.rm =TRUE)

plot.data<-rbind(data,AMA_mean,Young_mean)
plot.data[1:22,1:5]

col<-c(brewer.pal(8, 'RdYlGn')[8:7],brewer.pal(8, 'RdYlGn')[6:5],brewer.pal(8, 'GnBu')[2:7],
       brewer.pal(8, 'PuRd')[3:8], brewer.pal(8, 'YlOrRd')[8:5],ppCor[2:1])
length(col);show_col(col)

#for no mean line 
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Kids_20sample_TSS_TES.around.Met.AbsoluteDistance_CpG_6X_1.pdf",height = 8,width = 8)
#layout(matrix(1:23,3,3,byrow=FALSE))
#par(mar=c(2,2,1,2),oma=c(5,5,5,5))
#par(lend=2)
for(i in 1:nrow(plot.data[1:20,])){
  # i=1#test line
  body<-seq(from=151,by=5,length.out=100)
  x<-c(1:150,body,max(body):(max(body+149)))
  labelx<-c(seq(from=1,to=150,length.out=4)[-4],body[seq(from=1,to=length(body),length.out=6)],seq(from=max(body),to=max(body)+149,length.out=4)[-1])
  
  data <- plot.data[i,]
  if (rownames(plot.data)[i] == "AMA_1"){
    plot(x=x,y=data,type="l",xaxs="i",yaxs="i",ylim=c(0,100),axes=FALSE,xlab="",ylab="",col=col[i],lwd=1)
    box(bty="l")
    axis(side=1,at=labelx,labels=FALSE,tcl=-0.2)
    axis(side=2,at=pretty(c(0,100),4),tcl=-0.2,las=2)
    labeltext<-c("-15kb","","","0%","20%","40%","60%","80%","100%","","","15kb")
    mtext(side=1,at=labelx,labeltext,cex=0.7,line=0.3)
    mtext(side=2,"Ave. CG Methylation Level (%)",cex=0.8,line=2.5)
    mtext(side=1,at=labelx[c(4,9)],c("TSS","TES"),cex=0.7,line=1)
    legend("top",col=col[1:20],lty=1,ncol=4,bty="n",legend=c(sample_names) ,cex=0.7,lwd=2)
  }else {
    lines(x=x,y=data,col=col[i],lwd=4)
  }
}

segments(160, 0, 640, 0, col= 'darkblue',lty=1,lwd =10)
mtext(side=1,at=labelx[c(7)],c("Gene_body"),cex=0.7,line=1)
mtext(side=3,at=labelx[c(6)],c("Global DNA methylation level of CpGs around Gene body"),cex=0.7,line=1)
dev.off()

#including means line
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Kids_20sample_TSS_TES.around.Met.AbsoluteDistance_CpG_6X_include_mean_line.pdf",height = 8,width = 8)
#layout(matrix(1:23,3,3,byrow=FALSE))
#par(mar=c(2,2,1,2),oma=c(5,5,5,5))
#par(lend=2)
for(i in 1:nrow(plot.data)){
  # i=1#test line
  body<-seq(from=151,by=5,length.out=100)
  x<-c(1:150,body,max(body):(max(body+149)))
  labelx<-c(seq(from=1,to=150,length.out=4)[-4],body[seq(from=1,to=length(body),length.out=6)],seq(from=max(body),to=max(body)+149,length.out=4)[-1])
  
  data <- plot.data[i,]
  if (rownames(plot.data)[i] == "AMA_1"){
    plot(x=x,y=data,type="l",xaxs="i",yaxs="i",ylim=c(0,100),axes=FALSE,xlab="",ylab="",col=col[i],lwd=2)
    box(bty="l")
    axis(side=1,at=labelx,labels=FALSE,tcl=-0.2)
    axis(side=2,at=pretty(c(0,100),4),tcl=-0.2,las=2)
    labeltext<-c("-15kb","","","0%","20%","40%","60%","80%","100%","","","15kb")
    mtext(side=1,at=labelx,labeltext,cex=0.7,line=0.3)
    mtext(side=2,"Ave. CG Methylation Level (%)",cex=0.8,line=2.5)
    mtext(side=1,at=labelx[c(4,9)],c("TSS","TES"),cex=0.7,line=1)
    legend("right",col=col,lty=1,ncol=3,bty="n",legend=c(sample_names,"Young_mean","AMA_mean") ,cex=0.7,lwd=2)
  }else if(!(rownames(plot.data)[i] %in% c("Young_mean","AMA_mean"))) {
    lines(x=x,y=data,col=col[i],lwd=2)
  }else {
    lines(x=x,y=data,col=col[i],lwd=4)
  }
}

segments(160, 0, 640, 0, col= 'darkblue',lty=1,lwd =10)
mtext(side=1,at=labelx[c(7)],c("Gene_body"),cex=0.7,line=1)
mtext(side=3,at=labelx[c(6)],c("Global DNA methylation level of CpGs around Gene body"),cex=0.7,line=1)

dev.off()

pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Kids_20sample_TSS_TES.around.Met.AbsoluteDistance_CpG_6X_mean_line.pdf",height = 8,width = 8)
plot.data2<-plot.data[21:22,]
for(i in 1:nrow(plot.data2)){
  # i=1#test line
  body<-seq(from=151,by=5,length.out=100)
  x<-c(1:150,body,max(body):(max(body+149)))
  labelx<-c(seq(from=1,to=150,length.out=4)[-4],body[seq(from=1,to=length(body),length.out=6)],seq(from=max(body),to=max(body)+149,length.out=4)[-1])
  
  data <- plot.data2[i,]
  if (rownames(plot.data2)[i] == "AMA_mean"){
    plot(x=x,y=data,type="l",xaxs="i",yaxs="i",ylim=c(0,100),axes=FALSE,xlab="",ylab="",col=ppCor[i],lwd=4)
    box(bty="l")
    axis(side=1,at=labelx,labels=FALSE,tcl=-0.2)
    axis(side=2,at=pretty(c(0,100),4),tcl=-0.2,las=2)
    labeltext<-c("-15kb","","","0%","20%","40%","60%","80%","100%","","","15kb")
    mtext(side=1,at=labelx,labeltext,cex=0.7,line=0.3)
    mtext(side=2,"Ave.CG Methylation Level (%)",cex=0.8,line=2.5)
    mtext(side=1,at=labelx[c(4,9)],c("TSS","TES"),cex=0.7,line=1)
    legend("top",col=ppCor[1:2],lty=1,ncol=2,bty="n",legend=c("AMA_mean","Young_mean") ,cex=0.7,lwd=2)
  }else {
    lines(x=x,y=data,col=ppCor[i],lwd=4)
  }
}

segments(160, 0, 640, 0, col= 'orange',lty=1,lwd =20)
mtext(side=1,at=labelx[c(7)],c("Gene_body"),cex=0.7,line=1)
mtext(side=3,at=labelx[c(6)],c("Global DNA methylation level of CpGs around Gene body"),cex=0.7,line=1)

dev.off()
###################for Mother
sample_list<-c(AMA_Mother_files2,YOUNG_Mother_files2)
sample_names<-c(paste0("AMA_",c(1:10)),paste0("Young_",c(1:10)))

cnane_name<-c()
data <- NULL
for(i in 1:20){
  # i<-5
  sample<-sample_list[[i]]
  Sample_Name<-sample_names[[i]]
  print(as.character(Sample_Name))
  
  data0 <-read.table(sample, header=F)
  head(data0)
  data0$V1<-NULL
  data1<-apply(data0,2,function(x){mean(x,na.rm=TRUE)*100})
  data <- rbind(data,data1)
  cnane_name<-c(cnane_name,as.character(Sample_Name))
}  
nrow(data);length(cnane_name)
rownames(data) <- cnane_name
data[1:6,1:5]
AMA_mean<-colMeans(data[1:10,],na.rm =TRUE)
Young_mean<-colMeans(data[11:20,],na.rm =TRUE)

plot.data<-rbind(data,AMA_mean,Young_mean)
plot.data[1:22,1:5]

col<-c(brewer.pal(8, 'RdYlGn')[8:7],brewer.pal(8, 'RdYlGn')[6:5],brewer.pal(8, 'GnBu')[2:7],
       brewer.pal(8, 'PuRd')[3:8], brewer.pal(8, 'YlOrRd')[8:5],ppCor[2:1])
length(col);show_col(col)

#for no mean line 
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Mother_20sample_TSS_TES.around.Met.AbsoluteDistance_CpG_6X_1.pdf",height = 8,width = 8)
#layout(matrix(1:23,3,3,byrow=FALSE))
#par(mar=c(2,2,1,2),oma=c(5,5,5,5))
#par(lend=2)
for(i in 1:nrow(plot.data[1:20,])){
  # i=1#test line
  body<-seq(from=151,by=5,length.out=100)
  x<-c(1:150,body,max(body):(max(body+149)))
  labelx<-c(seq(from=1,to=150,length.out=4)[-4],body[seq(from=1,to=length(body),length.out=6)],seq(from=max(body),to=max(body)+149,length.out=4)[-1])
  
  data <- plot.data[i,]
  if (rownames(plot.data)[i] == "AMA_1"){
    plot(x=x,y=data,type="l",xaxs="i",yaxs="i",ylim=c(0,100),axes=FALSE,xlab="",ylab="",col=col[i],lwd=1)
    box(bty="l")
    axis(side=1,at=labelx,labels=FALSE,tcl=-0.2)
    axis(side=2,at=pretty(c(0,100),4),tcl=-0.2,las=2)
    labeltext<-c("-15kb","","","0%","20%","40%","60%","80%","100%","","","15kb")
    mtext(side=1,at=labelx,labeltext,cex=0.7,line=0.3)
    mtext(side=2,"Ave. CG Methylation Level (%)",cex=0.8,line=2.5)
    mtext(side=1,at=labelx[c(4,9)],c("TSS","TES"),cex=0.7,line=1)
    legend("top",col=col[1:20],lty=1,ncol=4,bty="n",legend=c(sample_names) ,cex=0.7,lwd=2)
  }else {
    lines(x=x,y=data,col=col[i],lwd=4)
  }
}

segments(160, 0, 640, 0, col= 'darkblue',lty=1,lwd =10)
mtext(side=1,at=labelx[c(7)],c("Gene_body"),cex=0.7,line=1)
mtext(side=3,at=labelx[c(6)],c("Global DNA methylation level of CpGs around Gene body"),cex=0.7,line=1)
dev.off()

#including means line
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Mother_20sample_TSS_TES.around.Met.AbsoluteDistance_CpG_6X_include_mean_line.pdf",height = 8,width = 8)
#layout(matrix(1:23,3,3,byrow=FALSE))
#par(mar=c(2,2,1,2),oma=c(5,5,5,5))
#par(lend=2)
for(i in 1:nrow(plot.data)){
  # i=1#test line
  body<-seq(from=151,by=5,length.out=100)
  x<-c(1:150,body,max(body):(max(body+149)))
  labelx<-c(seq(from=1,to=150,length.out=4)[-4],body[seq(from=1,to=length(body),length.out=6)],seq(from=max(body),to=max(body)+149,length.out=4)[-1])
  
  data <- plot.data[i,]
  if (rownames(plot.data)[i] == "AMA_1"){
    plot(x=x,y=data,type="l",xaxs="i",yaxs="i",ylim=c(0,100),axes=FALSE,xlab="",ylab="",col=col[i],lwd=2)
    box(bty="l")
    axis(side=1,at=labelx,labels=FALSE,tcl=-0.2)
    axis(side=2,at=pretty(c(0,100),4),tcl=-0.2,las=2)
    labeltext<-c("-15kb","","","0%","20%","40%","60%","80%","100%","","","15kb")
    mtext(side=1,at=labelx,labeltext,cex=0.7,line=0.3)
    mtext(side=2,"Ave. CG Methylation Level (%)",cex=0.8,line=2.5)
    mtext(side=1,at=labelx[c(4,9)],c("TSS","TES"),cex=0.7,line=1)
    legend("right",col=col,lty=1,ncol=3,bty="n",legend=c(sample_names,"Young_mean","AMA_mean") ,cex=0.7,lwd=2)
  }else if(!(rownames(plot.data)[i] %in% c("Young_mean","AMA_mean"))) {
    lines(x=x,y=data,col=col[i],lwd=2)
  }else {
    lines(x=x,y=data,col=col[i],lwd=4)
  }
}

segments(160, 0, 640, 0, col= 'darkblue',lty=1,lwd =10)
mtext(side=1,at=labelx[c(7)],c("Gene_body"),cex=0.7,line=1)
mtext(side=3,at=labelx[c(6)],c("Global DNA methylation level of CpGs around Gene body"),cex=0.7,line=1)

dev.off()

pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Mother_20sample_TSS_TES.around.Met.AbsoluteDistance_CpG_6X_mean_line.pdf",height = 8,width = 8)
plot.data2<-plot.data[21:22,]
for(i in 1:nrow(plot.data2)){
  # i=1#test line
  body<-seq(from=151,by=5,length.out=100)
  x<-c(1:150,body,max(body):(max(body+149)))
  labelx<-c(seq(from=1,to=150,length.out=4)[-4],body[seq(from=1,to=length(body),length.out=6)],seq(from=max(body),to=max(body)+149,length.out=4)[-1])
  
  data <- plot.data2[i,]
  if (rownames(plot.data2)[i] == "AMA_mean"){
    plot(x=x,y=data,type="l",xaxs="i",yaxs="i",ylim=c(0,100),axes=FALSE,xlab="",ylab="",col=ppCor[i],lwd=4)
    box(bty="l")
    axis(side=1,at=labelx,labels=FALSE,tcl=-0.2)
    axis(side=2,at=pretty(c(0,100),4),tcl=-0.2,las=2)
    labeltext<-c("-15kb","","","0%","20%","40%","60%","80%","100%","","","15kb")
    mtext(side=1,at=labelx,labeltext,cex=0.7,line=0.3)
    mtext(side=2,"Ave.CG Methylation Level (%)",cex=0.8,line=2.5)
    mtext(side=1,at=labelx[c(4,9)],c("TSS","TES"),cex=0.7,line=1)
    legend("top",col=ppCor[1:2],lty=1,ncol=2,bty="n",legend=c("AMA_mean","Young_mean") ,cex=0.7,lwd=2)
  }else {
    lines(x=x,y=data,col=ppCor[i],lwd=4)
  }
}

segments(160, 0, 640, 0, col= 'orange',lty=1,lwd =20)
mtext(side=1,at=labelx[c(7)],c("Gene_body"),cex=0.7,line=1)
mtext(side=3,at=labelx[c(6)],c("Global DNA methylation level of CpGs around Gene body"),cex=0.7,line=1)

dev.off()

###################for Father
sample_list<-c(AMA_father_files2,YOUNG_father_files2)
sample_names<-c(paste0("AMA_",c(1:10)),paste0("Young_",c(1:10)))

cnane_name<-c()
data <- NULL
for(i in 1:20){
  # i<-5
  sample<-sample_list[[i]]
  Sample_Name<-sample_names[[i]]
  print(as.character(Sample_Name))
  
  data0 <-read.table(sample, header=F)
  head(data0)
  data0$V1<-NULL
  data1<-apply(data0,2,function(x){mean(x,na.rm=TRUE)*100})
  data <- rbind(data,data1)
  cnane_name<-c(cnane_name,as.character(Sample_Name))
}  
nrow(data);length(cnane_name)
rownames(data) <- cnane_name
data[1:6,1:5]
AMA_mean<-colMeans(data[1:10,],na.rm =TRUE)
Young_mean<-colMeans(data[11:20,],na.rm =TRUE)

plot.data<-rbind(data,AMA_mean,Young_mean)
plot.data[1:22,1:5]

col<-c(brewer.pal(8, 'RdYlGn')[8:7],brewer.pal(8, 'RdYlGn')[6:5],brewer.pal(8, 'GnBu')[2:7],
       brewer.pal(8, 'PuRd')[3:8], brewer.pal(8, 'YlOrRd')[8:5],ppCor[2:1])
length(col);show_col(col)

#for no mean line 
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Father_20sample_TSS_TES.around.Met.AbsoluteDistance_CpG_6X_1.pdf",height = 8,width = 8)
#layout(matrix(1:23,3,3,byrow=FALSE))
#par(mar=c(2,2,1,2),oma=c(5,5,5,5))
#par(lend=2)
for(i in 1:nrow(plot.data[1:20,])){
  # i=1#test line
  body<-seq(from=151,by=5,length.out=100)
  x<-c(1:150,body,max(body):(max(body+149)))
  labelx<-c(seq(from=1,to=150,length.out=4)[-4],body[seq(from=1,to=length(body),length.out=6)],seq(from=max(body),to=max(body)+149,length.out=4)[-1])
  
  data <- plot.data[i,]
  if (rownames(plot.data)[i] == "AMA_1"){
    plot(x=x,y=data,type="l",xaxs="i",yaxs="i",ylim=c(0,100),axes=FALSE,xlab="",ylab="",col=col[i],lwd=1)
    box(bty="l")
    axis(side=1,at=labelx,labels=FALSE,tcl=-0.2)
    axis(side=2,at=pretty(c(0,100),4),tcl=-0.2,las=2)
    labeltext<-c("-15kb","","","0%","20%","40%","60%","80%","100%","","","15kb")
    mtext(side=1,at=labelx,labeltext,cex=0.7,line=0.3)
    mtext(side=2,"Ave. CG Methylation Level (%)",cex=0.8,line=2.5)
    mtext(side=1,at=labelx[c(4,9)],c("TSS","TES"),cex=0.7,line=1)
    legend("top",col=col[1:20],lty=1,ncol=4,bty="n",legend=c(sample_names) ,cex=0.7,lwd=2)
  }else {
    lines(x=x,y=data,col=col[i],lwd=4)
  }
}

segments(160, 0, 640, 0, col= 'darkblue',lty=1,lwd =10)
mtext(side=1,at=labelx[c(7)],c("Gene_body"),cex=0.7,line=1)
mtext(side=3,at=labelx[c(6)],c("Global DNA methylation level of CpGs around Gene body"),cex=0.7,line=1)
dev.off()

#including means line
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Father_20sample_TSS_TES.around.Met.AbsoluteDistance_CpG_6X_include_mean_line.pdf",height = 8,width = 8)
#layout(matrix(1:23,3,3,byrow=FALSE))
#par(mar=c(2,2,1,2),oma=c(5,5,5,5))
#par(lend=2)
for(i in 1:nrow(plot.data)){
  # i=1#test line
  body<-seq(from=151,by=5,length.out=100)
  x<-c(1:150,body,max(body):(max(body+149)))
  labelx<-c(seq(from=1,to=150,length.out=4)[-4],body[seq(from=1,to=length(body),length.out=6)],seq(from=max(body),to=max(body)+149,length.out=4)[-1])
  
  data <- plot.data[i,]
  if (rownames(plot.data)[i] == "AMA_1"){
    plot(x=x,y=data,type="l",xaxs="i",yaxs="i",ylim=c(0,100),axes=FALSE,xlab="",ylab="",col=col[i],lwd=2)
    box(bty="l")
    axis(side=1,at=labelx,labels=FALSE,tcl=-0.2)
    axis(side=2,at=pretty(c(0,100),4),tcl=-0.2,las=2)
    labeltext<-c("-15kb","","","0%","20%","40%","60%","80%","100%","","","15kb")
    mtext(side=1,at=labelx,labeltext,cex=0.7,line=0.3)
    mtext(side=2,"Ave. CG Methylation Level (%)",cex=0.8,line=2.5)
    mtext(side=1,at=labelx[c(4,9)],c("TSS","TES"),cex=0.7,line=1)
    legend("right",col=col,lty=1,ncol=3,bty="n",legend=c(sample_names,"Young_mean","AMA_mean") ,cex=0.7,lwd=2)
  }else if(!(rownames(plot.data)[i] %in% c("Young_mean","AMA_mean"))) {
    lines(x=x,y=data,col=col[i],lwd=2)
  }else {
    lines(x=x,y=data,col=col[i],lwd=4)
  }
}

segments(160, 0, 640, 0, col= 'darkblue',lty=1,lwd =10)
mtext(side=1,at=labelx[c(7)],c("Gene_body"),cex=0.7,line=1)
mtext(side=3,at=labelx[c(6)],c("Global DNA methylation level of CpGs around Gene body"),cex=0.7,line=1)

dev.off()

pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Father_20sample_TSS_TES.around.Met.AbsoluteDistance_CpG_6X_mean_line.pdf",height = 8,width = 8)
plot.data2<-plot.data[21:22,]
for(i in 1:nrow(plot.data2)){
  # i=1#test line
  body<-seq(from=151,by=5,length.out=100)
  x<-c(1:150,body,max(body):(max(body+149)))
  labelx<-c(seq(from=1,to=150,length.out=4)[-4],body[seq(from=1,to=length(body),length.out=6)],seq(from=max(body),to=max(body)+149,length.out=4)[-1])
  
  data <- plot.data2[i,]
  if (rownames(plot.data2)[i] == "AMA_mean"){
    plot(x=x,y=data,type="l",xaxs="i",yaxs="i",ylim=c(0,100),axes=FALSE,xlab="",ylab="",col=ppCor[i],lwd=4)
    box(bty="l")
    axis(side=1,at=labelx,labels=FALSE,tcl=-0.2)
    axis(side=2,at=pretty(c(0,100),4),tcl=-0.2,las=2)
    labeltext<-c("-15kb","","","0%","20%","40%","60%","80%","100%","","","15kb")
    mtext(side=1,at=labelx,labeltext,cex=0.7,line=0.3)
    mtext(side=2,"Ave.CG Methylation Level (%)",cex=0.8,line=2.5)
    mtext(side=1,at=labelx[c(4,9)],c("TSS","TES"),cex=0.7,line=1)
    legend("top",col=ppCor[1:2],lty=1,ncol=2,bty="n",legend=c("AMA_mean","Young_mean") ,cex=0.7,lwd=2)
  }else {
    lines(x=x,y=data,col=ppCor[i],lwd=4)
  }
}

segments(160, 0, 640, 0, col= 'orange',lty=1,lwd =20)
mtext(side=1,at=labelx[c(7)],c("Gene_body"),cex=0.7,line=1)
mtext(side=3,at=labelx[c(6)],c("Global DNA methylation level of CpGs around Gene body"),cex=0.7,line=1)

dev.off()


###############
#for all samples
sample_list<-c(AMA_Kids_files2,YOUNG_Kids_files2,AMA_Mother_files2,YOUNG_Mother_files2,AMA_father_files2,YOUNG_father_files2)
sample_names<-c(paste0("AMA_Kids",c(1:10)),paste0("YOUNG_Kids",c(1:10)),
                paste0("AMA_Mother",c(1:10)),paste0("YOUNG_Mother",c(1:10)),
                paste0("AMA_father",c(1:10)),paste0("YOUNG_father",c(1:10)))
cnane_name<-c()
data <- NULL
for(i in 1:60){
  # i<-5
  sample<-sample_list[[i]]
  Sample_Name<-sample_names[[i]]
  print(as.character(Sample_Name))
  
  data0 <-read.table(sample, header=F)
  head(data0)
  data0$V1<-NULL
  data1<-apply(data0,2,function(x){mean(x,na.rm=TRUE)*100})
  data <- rbind(data,data1)
  cnane_name<-c(cnane_name,as.character(Sample_Name))
}  
nrow(data);length(cnane_name)
rownames(data) <- cnane_name
data[1:6,1:5]
AMA_Kids_mean<-colMeans(data[1:10,],na.rm =TRUE)
YOUNG_Kids_mean<-colMeans(data[1:10,],na.rm =TRUE)
AMA_Mother_mean<-colMeans(data[1:10,],na.rm =TRUE)
YOUNG_Mother_mean<-colMeans(data[1:10,],na.rm =TRUE)
AMA_father_mean<-colMeans(data[1:10,],na.rm =TRUE)
YOUNG_father_mean<-colMeans(data[1:10,],na.rm =TRUE)

plot.data<-rbind(data,AMA_Kids_mean,YOUNG_Kids_mean,AMA_Mother_mean,YOUNG_Mother_mean,AMA_father_mean,YOUNG_father_mean)
plot.data[1:60,1:5];dim(plot.data)
col_fun1 <- colorRampPalette(c('beige','orange'))
col_fun2 <- colorRampPalette(c("grey",'purple'))
col_fun3 <- colorRampPalette(c('pink','red'))
col_fun4 <- colorRampPalette(c("lightblue",'blue'))
col_fun5 <- colorRampPalette(c("yellow",'brown'))
col_fun6 <- colorRampPalette(c("lightgreen",'darkgreen'))

#for no mean line 
col<-c(col_fun1(10),col_fun2(10),col_fun3(10),col_fun4(10),col_fun5(10),col_fun6(10),ppCor[c(3,5,1,2,8,4)])
length(col);show_col(col)

pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60sample_TSS_TES.around.Met.AbsoluteDistance_CpG_6X_1.pdf",height = 8,width = 8)
#layout(matrix(1:23,3,3,byrow=FALSE))
#par(mar=c(2,2,1,2),oma=c(5,5,5,5))
#par(lend=2)

for(i in 1:nrow(plot.data[1:60,])){
  # i=1#test line
  body<-seq(from=151,by=5,length.out=100)
  x<-c(1:150,body,max(body):(max(body+149)))
  labelx<-c(seq(from=1,to=150,length.out=4)[-4],body[seq(from=1,to=length(body),length.out=6)],seq(from=max(body),to=max(body)+149,length.out=4)[-1])
  
  data <- plot.data[i,]
  if (rownames(plot.data)[i] == "AMA_Kids1"){
    plot(x=x,y=data,type="l",xaxs="i",yaxs="i",ylim=c(0,100),axes=FALSE,xlab="",ylab="",col=col[i],lwd=1)
    box(bty="l")
    axis(side=1,at=labelx,labels=FALSE,tcl=-0.2)
    axis(side=2,at=pretty(c(0,100),4),tcl=-0.2,las=2)
    labeltext<-c("-15kb","","","0%","20%","40%","60%","80%","100%","","","15kb")
    mtext(side=1,at=labelx,labeltext,cex=0.7,line=0.3)
    mtext(side=2,"Ave. CG Methylation Level (%)",cex=0.8,line=2.5)
    mtext(side=1,at=labelx[c(4,9)],c("TSS","TES"),cex=0.7,line=1)
    legend("top",col=col[1:60],lty=1,ncol=4,bty="n",legend=c(sample_names) ,cex=0.7,lwd=2)
  }else {
    lines(x=x,y=data,col=col[i],lwd=1)
  }
}

segments(160, 0, 640, 0, col= 'darkblue',lty=1,lwd =10)
mtext(side=1,at=labelx[c(7)],c("Gene_body"),cex=0.7,line=1)
mtext(side=3,at=labelx[c(6)],c("Global DNA methylation level of CpGs around Gene body"),cex=0.7,line=1)
dev.off()

#including means line
show_col(col)
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60sample_TSS_TES.around.Met.AbsoluteDistance_CpG_6X_add_mean2.pdf",height = 8,width = 8)
#layout(matrix(1:23,3,3,byrow=FALSE))
#par(mar=c(2,2,1,2),oma=c(5,5,5,5))
#par(lend=2)
for(i in 1:nrow(plot.data)){
  # i=1#test line
  body<-seq(from=151,by=5,length.out=100)
  x<-c(1:150,body,max(body):(max(body+149)))
  labelx<-c(seq(from=1,to=150,length.out=4)[-4],body[seq(from=1,to=length(body),length.out=6)],seq(from=max(body),to=max(body)+149,length.out=4)[-1])
  
  data <- plot.data[i,]
  if (rownames(plot.data)[i] == "AMA_Kids1"){
    plot(x=x,y=data,type="l",xaxs="i",yaxs="i",ylim=c(0,100),axes=FALSE,xlab="",ylab="",col=col[i],lwd=2)
    box(bty="l")
    axis(side=1,at=labelx,labels=FALSE,tcl=-0.2)
    axis(side=2,at=pretty(c(0,100),4),tcl=-0.2,las=2)
    labeltext<-c("-15kb","","","0%","20%","40%","60%","80%","100%","","","15kb")
    mtext(side=1,at=labelx,labeltext,cex=0.7,line=0.3)
    mtext(side=2,"Ave. CG Methylation Level (%)",cex=0.8,line=2.5)
    mtext(side=1,at=labelx[c(4,9)],c("TSS","TES"),cex=0.7,line=1)
    legend("top",col=col,lty=1,ncol=4,bty="n",legend=c(sample_names,"AMA_Kids_mean","YOUNG_Kids_mean","AMA_Mother_mean","YOUNG_Mother_mean","AMA_father_mean","YOUNG_father_mean") ,cex=0.7,lwd=1)
  }else if(!(rownames(plot.data)[i] %in% c("AMA_Kids_mean","YOUNG_Kids_mean","AMA_Mother_mean","YOUNG_Mother_mean","AMA_father_mean","YOUNG_father_mean"))) {
    lines(x=x,y=data,col=col[i],lwd=2)
  }else {
    lines(x=x,y=data,col=col[i],lwd=4)
  }
}

segments(160, 0, 640, 0, col= 'orange',lty=1,lwd =10)
mtext(side=1,at=labelx[c(7)],c("Gene_body"),cex=0.7,line=1)
mtext(side=3,at=labelx[c(6)],c("Global DNA methylation level of CpGs around Gene body"),cex=0.7,line=1)

dev.off()
#only mean
col<-ppCor[c(3,5,1,2,8,4)]
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60sample_TSS_TES.around.Met.AbsoluteDistance_CpG_6X_only_mean3.pdf",height = 8,width = 8)
plot.data2<-plot.data[61:66,]
rownames(plot.data2)
for(i in 1:nrow(plot.data2)){
  # i=1#test line
  body<-seq(from=151,by=5,length.out=100)
  x<-c(1:150,body,max(body):(max(body+149)))
  labelx<-c(seq(from=1,to=150,length.out=4)[-4],body[seq(from=1,to=length(body),length.out=6)],seq(from=max(body),to=max(body)+149,length.out=4)[-1])
  
  data <- plot.data2[i,]
  if (rownames(plot.data2)[i] == "AMA_Kids_mean"){
    plot(x=x,y=data,type="l",xaxs="i",yaxs="i",ylim=c(0,100),axes=FALSE,xlab="",ylab="",col=ppCor[i],lwd=3)
    box(bty="l")
    axis(side=1,at=labelx,labels=FALSE,tcl=-0.2)
    axis(side=2,at=pretty(c(0,100),4),tcl=-0.2,las=2)
    labeltext<-c("-15kb","","","0%","20%","40%","60%","80%","100%","","","15kb")
    mtext(side=1,at=labelx,labeltext,cex=0.7,line=0.3)
    mtext(side=2,"Ave.CG Methylation Level (%)",cex=0.8,line=2.5)
    mtext(side=1,at=labelx[c(4,9)],c("TSS","TES"),cex=0.7,line=1)
    legend("top",col=col,lty=1,ncol=3,bty="n",legend=c("AMA_Kids_mean","YOUNG_Kids_mean","AMA_Mother_mean","YOUNG_Mother_mean","AMA_father_mean","YOUNG_father_mean") ,cex=0.7,lwd=4)
  }else {
    lines(x=x,y=data,col=col[i],lwd=3)
  }
}

segments(160, 0, 640, 0, col= 'black',lty=1,lwd =20)
mtext(side=1,at=labelx[c(7)],c("Gene_body"),cex=0.7,line=1)
mtext(side=3,at=labelx[c(6)],c("Global DNA methylation level of CpGs around Gene body"),cex=0.7,line=1)
dev.off()

###############
#for parents samples
sample_list<-c(AMA_Mother_files2,YOUNG_Mother_files2,AMA_father_files2,YOUNG_father_files2)
sample_names<-c(paste0("AMA_Mother",c(1:10)),paste0("YOUNG_Mother",c(1:10)),
                paste0("AMA_father",c(1:10)),paste0("YOUNG_father",c(1:10)))
cnane_name<-c()
data <- NULL
for(i in 1:40){
  # i<-5
  sample<-sample_list[[i]]
  Sample_Name<-sample_names[[i]]
  print(as.character(Sample_Name))
  
  data0 <-read.table(sample, header=F)
  head(data0)
  data0$V1<-NULL
  data1<-apply(data0,2,function(x){mean(x,na.rm=TRUE)*100})
  data <- rbind(data,data1)
  cnane_name<-c(cnane_name,as.character(Sample_Name))
}  
nrow(data);length(cnane_name)
rownames(data) <- cnane_name
data[1:6,1:5]
AMA_Mother_mean<-colMeans(data[1:10,],na.rm =TRUE)
YOUNG_Mother_mean<-colMeans(data[1:10,],na.rm =TRUE)
AMA_father_mean<-colMeans(data[1:10,],na.rm =TRUE)
YOUNG_father_mean<-colMeans(data[1:10,],na.rm =TRUE)

plot.data<-rbind(data,AMA_Mother_mean,YOUNG_Mother_mean,AMA_father_mean,YOUNG_father_mean)
plot.data[1:40,1:5];dim(plot.data)

col_fun1 <- colorRampPalette(c('pink','red'))
col_fun2 <- colorRampPalette(c("lightblue",'blue'))
#col_fun5 <- colorRampPalette(c("yellow",'brown'))
#col_fun6 <- colorRampPalette(c("lightgreen",'darkgreen'))
col_fun3 <- colorRampPalette(c("grey",'purple'))
col_fun4 <- colorRampPalette(c('beige','orange'))

#for no mean line 
col<-c(col_fun1(10),col_fun2(10),col_fun3(10),col_fun4(10),ppCor[c(1,2,5,3)])
length(col);show_col(col)

pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/parental_40sample_TSS_TES.around.Met.AbsoluteDistance_CpG_6X_1.pdf",height = 8,width = 8)
#layout(matrix(1:23,3,3,byrow=FALSE))
#par(mar=c(2,2,1,2),oma=c(5,5,5,5))
#par(lend=2)

for(i in 1:nrow(plot.data[1:40,])){
  # i=1#test line
  body<-seq(from=151,by=5,length.out=100)
  x<-c(1:150,body,max(body):(max(body+149)))
  labelx<-c(seq(from=1,to=150,length.out=4)[-4],body[seq(from=1,to=length(body),length.out=6)],seq(from=max(body),to=max(body)+149,length.out=4)[-1])
  
  data <- plot.data[i,]
  if (rownames(plot.data)[i] == "AMA_Mother1"){
    plot(x=x,y=data,type="l",xaxs="i",yaxs="i",ylim=c(0,100),axes=FALSE,xlab="",ylab="",col=col[i],lwd=1)
    box(bty="l")
    axis(side=1,at=labelx,labels=FALSE,tcl=-0.2)
    axis(side=2,at=pretty(c(0,100),4),tcl=-0.2,las=2)
    labeltext<-c("-15kb","","","0%","20%","40%","60%","80%","100%","","","15kb")
    mtext(side=1,at=labelx,labeltext,cex=0.7,line=0.3)
    mtext(side=2,"Ave. CG Methylation Level (%)",cex=0.8,line=2.5)
    mtext(side=1,at=labelx[c(4,9)],c("TSS","TES"),cex=0.7,line=1)
    legend("top",col=col[1:40],lty=1,ncol=4,bty="n",legend=c(sample_names) ,cex=0.7,lwd=2)
  }else {
    lines(x=x,y=data,col=col[i],lwd=1)
  }
}

segments(160, 0, 640, 0, col= 'darkblue',lty=1,lwd =10)
mtext(side=1,at=labelx[c(7)],c("Gene_body"),cex=0.7,line=1)
mtext(side=3,at=labelx[c(6)],c("Global DNA methylation level of CpGs around Gene body"),cex=0.7,line=1)
dev.off()

#including means line
show_col(col)
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/parental_40sample_TSS_TES.around.Met.AbsoluteDistance_CpG_6X_add_mean2.pdf",height = 8,width = 8)
#layout(matrix(1:23,3,3,byrow=FALSE))
#par(mar=c(2,2,1,2),oma=c(5,5,5,5))
#par(lend=2)
for(i in 1:nrow(plot.data)){
  # i=1#test line
  body<-seq(from=151,by=5,length.out=100)
  x<-c(1:150,body,max(body):(max(body+149)))
  labelx<-c(seq(from=1,to=150,length.out=4)[-4],body[seq(from=1,to=length(body),length.out=6)],seq(from=max(body),to=max(body)+149,length.out=4)[-1])
  
  data <- plot.data[i,]
  if (rownames(plot.data)[i] == "AMA_Mother1"){
    plot(x=x,y=data,type="l",xaxs="i",yaxs="i",ylim=c(0,100),axes=FALSE,xlab="",ylab="",col=col[i],lwd=2)
    box(bty="l")
    axis(side=1,at=labelx,labels=FALSE,tcl=-0.2)
    axis(side=2,at=pretty(c(0,100),4),tcl=-0.2,las=2)
    labeltext<-c("-15kb","","","0%","20%","40%","60%","80%","100%","","","15kb")
    mtext(side=1,at=labelx,labeltext,cex=0.7,line=0.3)
    mtext(side=2,"Ave. CG Methylation Level (%)",cex=0.8,line=2.5)
    mtext(side=1,at=labelx[c(4,9)],c("TSS","TES"),cex=0.7,line=1)
    legend("top",col=col,lty=1,ncol=4,bty="n",legend=c(sample_names,"AMA_Kids_mean","YOUNG_Kids_mean","AMA_Mother_mean","YOUNG_Mother_mean","AMA_father_mean","YOUNG_father_mean") ,cex=0.7,lwd=1)
  }else if(!(rownames(plot.data)[i] %in% c("AMA_Mother_mean","YOUNG_Mother_mean","AMA_father_mean","YOUNG_father_mean"))) {
    lines(x=x,y=data,col=col[i],lwd=2)
  }else {
    lines(x=x,y=data,col=col[i],lwd=4)
  }
}

segments(160, 0, 640, 0, col= 'orange',lty=1,lwd =10)
mtext(side=1,at=labelx[c(7)],c("Gene_body"),cex=0.7,line=1)
mtext(side=3,at=labelx[c(6)],c("Global DNA methylation level of CpGs around Gene body"),cex=0.7,line=1)

dev.off()
#only mean
col<-ppCor[c(1,2,5,3)]
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/parental_40sample_TSS_TES.around.Met.AbsoluteDistance_CpG_6X_only_mean3.pdf",height = 8,width = 8)
plot.data2<-plot.data[41:44,]
rownames(plot.data2)
for(i in 1:nrow(plot.data2)){
  # i=1#test line
  body<-seq(from=151,by=5,length.out=100)
  x<-c(1:150,body,max(body):(max(body+149)))
  labelx<-c(seq(from=1,to=150,length.out=4)[-4],body[seq(from=1,to=length(body),length.out=6)],seq(from=max(body),to=max(body)+149,length.out=4)[-1])
  
  data <- plot.data2[i,]
  if (rownames(plot.data2)[i] == "AMA_Mother_mean"){
    plot(x=x,y=data,type="l",xaxs="i",yaxs="i",ylim=c(0,100),axes=FALSE,xlab="",ylab="",col=ppCor[i],lwd=3)
    box(bty="l")
    axis(side=1,at=labelx,labels=FALSE,tcl=-0.2)
    axis(side=2,at=pretty(c(0,100),4),tcl=-0.2,las=2)
    labeltext<-c("-15kb","","","0%","20%","40%","60%","80%","100%","","","15kb")
    mtext(side=1,at=labelx,labeltext,cex=0.7,line=0.3)
    mtext(side=2,"Ave.CG Methylation Level (%)",cex=0.8,line=2.5)
    mtext(side=1,at=labelx[c(4,9)],c("TSS","TES"),cex=0.7,line=1)
    legend("top",col=col,lty=1,ncol=4,bty="n",legend=c("AMA_Mother_mean","YOUNG_Mother_mean","AMA_father_mean","YOUNG_father_mean") ,cex=0.7,lwd=4)
  }else {
    lines(x=x,y=data,col=col[i],lwd=3)
  }
}

segments(160, 0, 640, 0, col= 'black',lty=1,lwd =20)
mtext(side=1,at=labelx[c(7)],c("Gene_body"),cex=0.7,line=1)
mtext(side=3,at=labelx[c(6)],c("Global DNA methylation level of CpGs around Gene body"),cex=0.7,line=1)
dev.off()
