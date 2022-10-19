library(scales)
library(ggsci)
library(grid)
library(reshape2)
library(stringr)
library(pheatmap)
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

type<-c("Mom","Kid")
for (tag_name in type){
#tag_name<-"Kid";
Conditions<-"group";sampleA <-"AMA";sampleB <-"Young" 
file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/file9_all_filter_",tag_name,"_verification_",sampleA,"_vs_",sampleB,".log2_normalized_count_add_gene_information.txt")
merge_data<-read.table(file=file, sep="\t",header = T,row.names=1)

head(merge_data) 
merge_data$threshold = as.factor(ifelse(merge_data$pvalue < 0.05 & abs(merge_data$log2FoldChange)>= 0.585, ifelse(merge_data$log2FoldChange>= 0.585 ,'Up','Down'),'NoSignifi'))
table(merge_data$threshold)
#Down  NoSignifi        Up 
# 2029     14348      1128 
merge_data$threshold<-factor(merge_data$threshold,level=c("Up","Down","NoSignifi"))
range(merge_data$log2FoldChange)#-5.598469  4.432278
range(-log10(merge_data$pvalue))#0.00000 20.19249


plot_vocano<-ggplot(data = merge_data, aes(x = log2FoldChange, y = -log10(pvalue), colour=threshold)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("red","blue","grey")) +
 # scale_x_continuous(breaks=seq(-6, 6, 1))+
  xlim(-6,6)+
  ylim(0,ceiling(max(-log10(merge_data$pvalue))))+
  geom_vline(xintercept=c(-0.585,0.585),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(Foldchange(AMA/Young)) in gene expression level",y="-log10 (p value)",
       title=paste(tag_name,"_",sampleA,"_vs_",sampleB,":Up_DEGs: ",length(which(merge_data$threshold == "Up"))," & Down_DEG: ",length(which(merge_data$threshold == "Down"))," [p value<0.05 & FoldChange > 1.5]",sep="")) +
  theme_bw()+
  theme(panel.border = element_rect(colour="black",fill=NA),panel.grid.major =element_blank())+
  theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5), legend.position="right", legend.title = element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
ggsave(plot_vocano,file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_",sampleA,"_vs_",sampleB,"_vocano_DEGs.pdf"), width = 200, height = 200, units = "mm")


#绘制各组特异高表达基因热图
res_de_up <- subset(merge_data, pvalue<0.05&log2FoldChange>=0.585 )
res_de_dw <- subset(merge_data, pvalue<0.05&log2FoldChange<=(-1)*0.585 )
de_id_whole = rbind(res_de_up, res_de_dw)
nrow(de_id_whole);head(de_id_whole)
mat<-de_id_whole[,1:10]
anno1 <- data.frame(sample=colnames(mat))
anno1$group<-ifelse(grepl("A",anno1$sample),"AMA","Young")
rownames(anno1)<-anno1$sample
ann_colors = list(group=c(Young=ppCor[2],AMA=ppCor[1]))
range(mat)# 0.00000 18.32706
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01)) 
# 做热图： 
heat_plot<-pheatmap(mat, scale = "row", 
                    color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), 
                    cluster_row =FALSE,
                    cluster_col =FALSE,                    
                    legend_breaks=seq(-2,2,1), 
                    breaks=bk,annotation_colors = ann_colors, 
                    gaps_row = nrow(res_de_up),gaps_col =5,cutree_col = 2,
                    show_rownames = F,
                    annotation_col = anno1[2],
                    main=paste0(tag_name,":DEGs_pvalue <0.05& FoldChange>=1.5")
)

pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/heatmap_All_filter_",tag_name,"_DEGs_p005_FC1.5_AMA_Young.pdf"),width = 7, height =7)
print(heat_plot)
dev.off()
}
