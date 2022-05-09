library(ggplot2)

data <- data.frame(log2FC=log2FC,p.adj=p.adj)
rownames(data) <- rownames(DEG)
data$sig[data$p.adj > 0.05 | (data$log2FC < 0.5 & data$log2FC > -0.5)] <- "notSig"
data$sig[data$p.adj <= 0.05 & data$log2FC >= 0.5] <- "up"
data$sig[data$p.adj <= 0.05 & data$log2FC <= -0.5] <- "down"

x_lim <- max(abs(log2FC))+0.5
y_lim <- max(-1*log10(p.adj))+1

theme <- theme(panel.grid=element_blank(),
               panel.background = element_rect(fill="white",color="black"),
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 15),
               axis.title.x = element_text(size = 15),
               axis.title.y = element_text(size = 15),
               axis.text.x = element_text(size = 13),
               axis.text.y = element_text(size = 13))

p <- ggplot(data,aes(x=log2FC,y=-log10(p.adj),color = sig))+geom_point(size=0.9)+
  labs(x="log2(fold change)",y="-log10(p.adj)",title="Male ob/ob compared to WT") +
  scale_color_manual("",values =c(down="#0072B5",notSig="grey",up="#BC3C28"))+
  scale_x_continuous(limits = c(-2,2),breaks=seq(-2,2,1))+
  scale_y_continuous(limits = c(0,y_lim),breaks=seq(0,150,50),expand = c(0,0))+
  geom_hline(yintercept=-log10(0.05),linetype=3)+
  geom_vline(xintercept=c(-0.5,0.5),linetype=3) + theme+
  theme(legend.position=c(0.85,0.85),
        plot.title = element_text(hjust = 0.5,size = 15))+
  #nudge_x nudge_y 更改标签与点的距离
  geom_text_repel(data = data[gene,],
                  aes(x=data[gene,]$log2FC+0.05,y=-log10(data[gene,]$p.adj)+0.05,label = rownames(data[gene,])),
                  size = 3.5,color="black",nudge_x=+0.5,nudge_y=3)
pdf("male_volcano_plot.pdf",6,6)
p
dev.off()




DEG_data <- read.xlsx("ob_female_2_allcels_DEGS.xlsx", sheet = 3, colNames = T, rowNames = T)

DEG_data$padj_manual <- DEG_data$padj
#DEG_data$padj_manual[DEG_data$padj<2.2e-16] <- 2.2e-16
DEG_data$logP <- -log10(DEG_data$padj) # 对差异基因矫正后p-value进行log10()转换
dim(DEG_data)
## [1] 9021    7
#将基因分为三类：not-siginficant，up，dowm
#将adj.P.value小于0.05，logFC大于2的基因设置为显著上调基因
#将adj.P.value小于0.05，logFC小于-2的基因设置为显著上调基因
DEG_data$Group <- "notSig"
DEG_data$Group[which((DEG_data$padj < 0.05) & DEG_data$log2FoldChange > 0.1)] = "up"
DEG_data$Group[which((DEG_data$padj < 0.05) & DEG_data$log2FoldChange < -0.1)] = "down"
table(DEG_data$Group)
## 
# down notSig     up 
#587   5818   2616
#火山图中添加点(数据构建)
label <- read.xlsx("ob_female_2_allcels_DEGS.xlsx", sheet = 4, colNames = T, rowNames = T)
up_label <- label[label$`Up/Down`=="Up" & label$log2FoldChange > 0.1 , ]
down_label <- label[label$`Up/Down`=="Down"& label$log2FoldChange < -0.1,]
deg_label_gene <- data.frame(gene = c(rownames(up_label),rownames(down_label)),
                             label = c(rownames(up_label),rownames(down_label)))
DEG_data$gene <- rownames(DEG_data)
DEG_data <- merge(DEG_data,deg_label_gene,by = 'gene',all = T)

#不添加label
library(ggpubr)
library(ggplot2)
library(ggthemes)
library(ggrepel)


#添加特定基因label
#   #nudge_x nudge_y 更改标签与点的距离
p1 <- ggscatter(DEG_data,x = "log2FoldChange",y = "logP",
                color = "Group",
                palette = c("#377EB8","gray","#E41A1C"),
                label = DEG_data$label,
                #label.select = deg_label_gene$gene,
                repel = T,
                ylab = "-log10(adjusted P-value)",
                title="Female ob/ob compared to WT",
                size = 1, nudge_x=+0.5,nudge_y=3) + 
  theme_base()+
  theme(element_line(size = 0),element_rect(size = 1.5))+ #坐标轴线条大小设置
  #scale_y_continuous(limits = c(0,8))+
  scale_x_continuous(limits = c(-1,1) ,breaks=seq(-1,1,0.2))+
  geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
  geom_vline(xintercept = c(-0.1,0.1),linetype = "dashed")+ 
  theme(plot.title = element_text(hjust = 0.5)) 

