rm(list = ls())
library(data.table)
library(survminer)
library(survival)
library(openxlsx)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(stringr)
library(Mime1)
library(patchwork) 
library(tidyverse)
library(ggpubr)
library(compareGroups)
library(reshape2)
library(tidyr)
library(ggpubr)
library(GSVA)
library(ggsci)

#读入expr
#分组信息
expr <- read.csv("data/cohort/TCGA_COAD_expr.csv", row.names = 1)
colnames(expr) <- gsub("\\.", "-", colnames(expr))
TCGA_RS <- read.csv("data/revise mime/TCGA_RS.csv")
RS <- TCGA_RS[,c(2,6),drop=F]
# 确保样本顺序匹配
common_samples <- intersect(colnames(expr), RS$ID)
expr <- expr[, common_samples]
RS <- RS[match(common_samples, RS$ID), ]


####### A.Cibersort
source('data/input/CIBERSORT.R')
df <- expr
df[is.na(df)] <- 0
df$gene<-rownames(df)
df <- df[,c(434,1:433)]
write.table(df,"data/revise group/Cibersort_exprMat.txt",sep = "\t",row.names = F)

mixture_file = 'data/revise group/Cibersort_exprMat.txt' # 表达矩阵
sig_matrix <- "data/input/LM22.txt" ## LM22 基因特征集

my_cibersort <- CIBERSORT(sig_matrix, ## LM22基因集 对应22种immnue cells
                          mixture_file, ## 表达矩阵(是否log2不重要;无需normalized;若RNA-seq:推荐TPM)
                          perm=1000, ##置换检验1000次
                          QN=TRUE) ##芯片=T;RNA-seq=F; QN:Quntile Normalization

write.csv(my_cibersort,file = "data/revise group/CIBERSORT.csv",row.names = T,quote = F)

#作图数据
my_cibersort <- read.csv("data/revise group/CIBERSORT.csv")
my_cibersort <- my_cibersort[my_cibersort[,'P.value'] < 0.05,]
my_cibersort <- my_cibersort[,1:23]
colnames(my_cibersort)<-gsub("\\ ","\\_",colnames(my_cibersort))

dat <- my_cibersort
rownames(dat) <- dat$X
dat <-dat[,-1]
dat <- dat %>% as.data.frame() %>%
  rownames_to_column("ID") %>% 
  gather(key = Cell_type,value = Proportion,-ID) 

TCGA_RS <- TCGA_RS[TCGA_RS$ID %in% dat$ID, ]
rownames(TCGA_RS) <- TCGA_RS$ID
Group <- TCGA_RS[,6,drop=F]
dat <- merge(Group,dat,by.x=0,by.y = 1) 

#将细胞按照比例大小排序
a = dat %>% 
  group_by(Cell_type) %>% 
  summarise(m = median(Proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(Cell_type)
dat$Cell_type = factor(dat$Cell_type,levels = a)

#筛选前10个细胞作图 
a <- a[1:10]
dat <- dat[dat$Cell_type %in% a,]
a = dat %>% 
  group_by(Cell_type) %>% 
  summarise(m = median(Proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(Cell_type)

ggplot(dat,aes(Cell_type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = NULL, y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text = element_text(size = 11),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,hjust = 1),
        plot.title = element_text(size = 12,face = "bold",hjust = 0.3),
        legend.title = element_text(size = 12,face = "bold"))+
  scale_fill_manual(values = c("#FF7F0E", "#66A61E"))+ 
  stat_compare_means(aes(group = Group,label = ..p.signif..))

ggsave(filename = "output/revise group/A.Immune_Cibersort.pdf",width = 6,height = 4)



###### B.APS Antigen presentation
#提取参考文献中APS基因
imsig <- read.csv("data/input/NIHMS958212-supplement-7.csv")
imsig <- imsig[1:78,1:8]
ap <- imsig[imsig$Super.Category == "Antigen presentation",]
ap_gene <- ap$Gene
#16个APS基因的表达apExpr
apExpr <- expr[rownames(expr) %in% ap_gene,]

APS <- list('APS'= ap_gene)
GSVA <- gsva(expr = as.matrix(expr),
             gset.idx.list = APS,
             method="gsva")
APS <- data.frame(APS=t(GSVA))%>%rownames_to_column('ID')

dat <- merge(APS, RS, by="ID")
p <- 
  ggplot(dat,aes(Group,APS,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black",width = 0.5) + 
  theme_bw() + 
  labs(x = NULL, y = "Relative APS") +
  theme(legend.position = "top") + 
  theme(axis.text = element_text(size = 11),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 0.9),
        plot.title = element_text(size = 12,face = "bold",hjust = 0.3),
        legend.title = element_text(size = 12,face = "bold"))+
  scale_fill_manual(values = c("#FF7F0E", "#66A61E"))+ 
  stat_compare_means(aes(group = Group,label = ..p.signif..))
p

ggsave(p,filename = "output/revise group/B.APS.pdf",width = 3,height = 3)



###### C.TIS
TISG <- c("CCL5", "CD27", "CD274", "CD276", 
          "CD8A", "CMKLR1", "CXCL9", "HLA-DQA1",
          "HLA-DRB1", "HLA-E", "IDO1",
          "LAG3", "NKG7", "PDCD1LG2", 
          "PSMB10", "STAT1", "TIGIT","CXCR6")
GSVA <- gsva(expr=as.matrix(expr), 
             gset.idx.list=list(TISG), method="gsva")
TIS <- data.frame(TIS=t(GSVA))%>%rownames_to_column('ID')

dat <-  merge(TIS, RS, by="ID")%>% 
  column_to_rownames("ID")

p <- 
  ggplot(dat,aes(Group,TIS,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black",width = 0.5) + 
  theme_bw() + 
  labs(x = NULL, y = "Relative TIS") +
  theme(legend.position = "top") + 
  theme(axis.text = element_text(size = 11),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 0.9),
        plot.title = element_text(size = 12,face = "bold",hjust = 0.3),
        legend.title = element_text(size = 12,face = "bold"))+
  scale_fill_manual(values = c("#FF7F0E", "#66A61E"))+ 
  stat_compare_means(aes(group = Group,label = ..p.signif..))
p

ggsave(p,filename = "output/revise group/C.TIS.pdf",width = 3,height = 3)




###### D.immune molecule
### 免疫分子
mm <- data.table::fread("data/input/immune_related_molecular_human.txt")

## 分组信息
rownames(RS) <- RS$ID

gene <- intersect(rownames(expr),mm$GENE.SYMBLE)
mm <- mm[mm$GENE.SYMBLE %in% gene,]
mm_df <- expr[gene, ] %>%
  t()%>%
  as.data.frame() %>%
  merge(RS,.,by.x = 0,by.y = 0) %>%
  column_to_rownames("Row.names")

#--比较三组的差异-----------------------------------------
mm_df$Group <- factor(mm_df$Group)
dd <- descrTable(Group~.,data = mm_df,method = 2)
print(dd)
pvals <- getResults(dd, "p.overall")
p.adjust(pvals, method = "BH")

#---绘图boxplot-----------------------------------------------
gene <-  colnames(mm_df)[which(pvals<0.5)+1]

mm1 <- mm[mm$GENE.SYMBLE %in% gene,]
mm_df1 <- cbind(Group=mm_df$Group,mm_df[,gene])
ggdata <- melt(mm_df1,id.vars = "Group",
               variable.name = "Gene",value.name = "expr") %>% 
  merge(.,mm1,by.x = 2,by.y = 1,all = T) 
ggdata <- ggdata[order(ggdata$Type,ggdata$Gene),]
ggdata$Gene <- factor(ggdata$Gene)

ggdata_ca <- ggdata[ggdata$Type == "Co-stimulatory immune checkpoint targets",]
ggdata_ci <- ggdata[ggdata$Type == "Co-inhibitory immune checkpoint targets",]
ggdata_mhc <- ggdata[ggdata$Type %in% c("MHC Class II" ,"MHC Class I"),]

#共刺激分子
gg_ca <- ggplot(ggdata_ca,aes(Gene,expr,fill=Group))+
  geom_boxplot()+
  scale_fill_manual(values = c("#FF7F0E", "#66A61E"))+
  theme_bw()+
  labs(title="Co-stimulatory ICTs",x=NULL,y="Expression")+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,hjust = 0.9),
        plot.title = element_text(size = 12,face = "bold",hjust = 0.3),
        legend.title = element_text(size = 12,face = "bold"))+
  stat_compare_means(aes(label=..p.signif..))
#共抑制分子
gg_ci <- ggplot(ggdata_ci,aes(Gene,expr,fill=Group))+
  geom_boxplot()+
  scale_fill_manual(values = c("#FF7F0E", "#66A61E"))+
  theme_bw()+
  labs(title="Co-inhibitory ICTs",x=NULL,y="Expression")+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,hjust = 0.9),
        plot.title = element_text(size = 12,face = "bold",hjust = 0.3),
        legend.title = element_text(size = 12,face = "bold"))+
  stat_compare_means(aes(label=..p.signif..))
#MHC
gg_mhc <- ggplot(ggdata_mhc,aes(Gene,expr,fill=Group))+
  geom_boxplot()+
  scale_fill_manual(values = c("#FF7F0E", "#66A61E"))+
  theme_bw()+
  labs(title="MHC molecules",x=NULL,y="Expression")+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,hjust = 0.9),
        plot.title = element_text(size = 12,face = "bold",hjust = 0.3,vjust=0.8),
        legend.title = element_text(size = 12,face = "bold"))+
  stat_compare_means(aes(label=..p.signif..))

gg <- ggarrange(gg_ca, gg_ci, gg_mhc, ncol = 1,nrow =3)
gg

ggsave(gg,filename = "output/revise group/D.Immune_melocules.pdf",width = 6,height = 8)



###### E.CIC
presurv <- RS
rownames(presurv) <- presurv$ID

rr2 <- presurv[2,drop=F] %>%
  .[order(.$Group),,drop=F]

exp <- expr %>%
  .[,match(rownames(rr2),colnames(.))]

#参考文件
load('data/input/CIC_Marker.rda')
ssGSEA_output <- gsva(as.matrix(exp),CIC_Marker, method = "ssgsea")
ssGSEA_output <- as.data.frame(t(ssGSEA_output))

CIC_input <- merge(ssGSEA_output, rr2, by = 0) %>%
  column_to_rownames('Row.names')

input <- pivot_longer(CIC_input,cols = 1:22, names_to = 'Steps', values_to = 'Value')
input$Steps <- factor(input$Steps,levels = c(colnames(CIC_input)[1],colnames(CIC_input)[2],
                                             colnames(CIC_input)[3],sort(colnames(CIC_input)[4:19]),
                                             colnames(CIC_input)[20],colnames(CIC_input)[21],
                                             colnames(CIC_input)[22]))

col <- c(pal_npg("nrc", alpha =0.9)(9)[1],pal_npg("nrc", alpha =0.9)(9)[2],pal_npg("nrc", alpha =0.9)(9)[3],
         rep(pal_npg("nrc", alpha =0.9)(9)[4],16),pal_npg("nrc", alpha =0.9)(9)[9],
         pal_jama("default", alpha =0.9)(7)[2],pal_jama("default", alpha =0.9)(7)[1])

ggplot(input,aes(Steps,Value,fill = Group))+
  geom_boxplot(outlier.colour = 'grey30',
               outlier.size = 0.2,
               outlier.alpha = 0.5,
               size=0.4)+
  stat_compare_means(label = 'p.signif',method = "kruskal.test",
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", " ")))+
  theme_bw()+ coord_flip()+
  theme(legend.position = 'right',
        panel.grid = element_blank(),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=11),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  scale_fill_manual(values = c("#FF7F0E", "#66A61E"))+
  labs(title = 'Cancer Immunity Cycle',fill='Group',y=NULL,x=NULL)

ggsave(filename = 'output/revise group/E.CIC.pdf',height = 8,width = 6)



###### F.IPS
#自定义函数，用于计算 Immunophenoscore(IPS)
ipsmap<- function (x) {
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=0)
    }
  }
  return(ips)
}

#分配颜色
my_palette <- colorRampPalette(c("#88C4AC", "white", "#E7597C"))(n = 1000)
mapcolors<-function (x) {
  za<-NULL
  if (x>=3) {
    za=1000
  } else {
    if (x<=-3) {
      za=1
    } else {
      za=round(166.5*x+500.5,digits=0)
    }
  }
  return(my_palette[za])
}
my_palette2 <- colorRampPalette(c("#FFC77D", "white"))(n = 1000)
mapbw<-function (x) {
  za2<-NULL
  if (x>=2) {
    za2=1000
  } else {
    if (x<=-2) {
      za2=1
    } else {
      za2=round(249.75*x+500.5,digits=0)
    }
  }
  return(my_palette2[za2])
}


#读取表达矩阵
gene_expression <-expr
## 修改几个基因名
## IPS基因权重文件中的基因名跟常用的gene symbol不同，我们把表达矩阵里的gene symbol修改成跟IPS基因列表一致的基因名。
is.element("CAVIN2",rownames(gene_expression))
#TRUE
is.element("DARS1",rownames(gene_expression))
#TRUE
is.element("IARS1",rownames(gene_expression))
#TRUE
if(is.element("CAVIN2",rownames(gene_expression))) {
  rownames(gene_expression) <- gsub("CAVIN2","SDPR",rownames(gene_expression)) # 根据https://www.genecards.org上的相关信息进行同名转换
}
if(is.element("DARS1",rownames(gene_expression))) {
  rownames(gene_expression) <- gsub("DARS1","DARS",rownames(gene_expression)) # 根据https://www.genecards.org上的相关信息进行同名转换
}
if(is.element("IARS1",rownames(gene_expression))) {
  rownames(gene_expression) <- gsub("IARS1","IARS",rownames(gene_expression)) # 根据https://www.genecards.org上的相关信息进行同名转换
}

## 读取IPS相关基因以及其权重
IPSG <- read.table("data/input/easy_input_IPS_genes.txt",header=TRUE, sep="\t",check.names=FALSE,stringsAsFactors = FALSE,row.names = NULL)
head(IPSG)
unique_ips_genes <- as.vector(unique(IPSG$NAME))

## 初始化数据
IPS <- MHC <- CP <- EC <- SC <- AZ <- NULL

# 获取表达谱里的基因名
GVEC <- row.names(gene_expression)

# 获取IPS基因文件里的基因名
VEC <- as.vector(IPSG$GENE)

# 匹配基因并找到缺失基因
ind <- which(is.na(match(VEC,GVEC)))

MISSING_GENES <- VEC[ind]

dat <- IPSG[ind,]
if (length(MISSING_GENES) > 0) { # 若存在缺失基因报出（如果上面两个基因不修改，这里是要报缺失的，会导致最终结果错误）
  message(paste0("--differently named or missing genes: ",paste(MISSING_GENES,collapse = ", ")))
  print(IPSG[ind,])
  message("please check data and make sure all genes matches!")
} else {
  message("--all genes matched!") # 请确保所有基因都匹配!!!
}


## 循环每个样本
tumsam <- colnames(gene_expression)
sample_names <- tumsam
for (i in 1:length(sample_names)) { 
  GE <- gene_expression[[i]]
  mGE <- mean(GE)
  sGE <- sd(GE)
  Z1 <- (gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
  W1 <- IPSG$WEIGHT
  WEIGHT <- MIG <- NULL
  k <- 1
  for (gen in unique_ips_genes) {
    MIG[k] <- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
    WEIGHT[k] <- mean(W1[which (as.vector(IPSG$NAME)==gen)])
    k <- k + 1
  }
  WG <- MIG * WEIGHT
  MHC[i] <- mean(WG[1:10])
  CP[i] <- mean(WG[11:20])
  EC[i] <- mean(WG[21:24])
  SC[i] <- mean(WG[25:26])
  AZ[i] <- sum(MHC[i],CP[i],EC[i],SC[i])
  IPS[i] <- ipsmap(AZ[i])
  
  # 绘制Immunophenogram
  data_a <- data.frame (start = c(0,2.5,5,7.5,10,15,seq(20,39),0,10,20,30), end = c(2.5,5,7.5,10,15,seq(20,40),10,20,30,40), y1=c(rep(2.6,26),rep(0.4,4)),y2=c(rep(5.6,26),rep(2.2,4)),z=c(MIG[c(21:26,11:20,1:10)],EC[i],SC[i],CP[i],MHC[i]),vcol=c(unlist(lapply(MIG[c(21:26,11:20,1:10)],mapcolors)), unlist(lapply(c(EC[i],SC[i],CP[i],MHC[i]),mapbw))), label = c(unique_ips_genes[c(21:26,11:20,1:10)],"EC","SC","CP","MHC"))
  data_a$label <- factor(data_a$label, levels=unique(data_a$label))
  plot_a1 <- ggplot() + 
    geom_rect(data=data_a, mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label), size=0.5,color="black", alpha=1) +
    coord_polar()  + 
    scale_y_continuous(limits = c(0, 6)) +
    scale_fill_manual(values =as.vector(data_a$vcol),guide=FALSE) +
    theme_bw() + 
    theme(panel.spacing = unit(0, 'mm'), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "white"),
          axis.text=element_blank(), 
          axis.ticks= element_blank()) + 
    geom_text(aes(x=5, y=1.3, label="EC"), size=4) +
    geom_text(aes(x=15, y=1.3, label="SC"), size=4) + 
    geom_text(aes(x=25, y=1.3, label="CP"), size=4) + 
    geom_text(aes(x=35, y=1.3, label="MHC"), size=4)
  plot_a2 <- plot_a1 +
    geom_text(aes(x=1.25, y=4.1, label="+ Act CD4"), angle=78.75, size=4) +
    geom_text(aes(x=3.75, y=4.1, label="+ Act CD8"),angle=56.25, size=4) +
    geom_text(aes(x=6.25, y=4.1, label="+ Tem CD4"), angle=33.75,size=4) +
    geom_text(aes(x=8.75, y=4.1, label="+ Tem CD8"), angle=11.25,size=4) +
    geom_text(aes(x=17.5, y=4.1, label="- MDSC"), angle=-67.5,size=4) +
    geom_text(aes(x=12.5, y=4.1, label="- Treg"), angle=-22.5,size=4)
  plot_a3 <- plot_a2 +
    geom_text(aes(x=20.5, y=4.1, label="PD-1 -"), angle=85.5, size=4) +
    geom_text(aes(x=21.5, y=4.1, label="CTLA4 -"), angle=76.5, size=4) +
    geom_text(aes(x=22.5, y=4.1, label="LAG3 -"), angle=67.5, size=4) +
    geom_text(aes(x=23.5, y=4.1, label="TIGIT -"), angle=58.5, size=4) +
    geom_text(aes(x=24.5, y=4.1, label="TIM3 -"), angle=49.5, size=4) +
    geom_text(aes(x=25.5, y=4.1, label="PD-L1 -"), angle=40.5, size=4) +
    geom_text(aes(x=26.5, y=4.1, label="PD-L2 -"), angle=31.5, size=4) +
    geom_text(aes(x=27.5, y=4.1, label="CD27 +"), angle=22.5, size=4) +
    geom_text(aes(x=28.5, y=4.1, label="ICOS +"), angle=13.5, size=4) +
    geom_text(aes(x=29.5, y=4.1, label="IDO1 -"), angle=4.5, size=4)
  plot_a4 <- plot_a3 +
    geom_text(aes(x=30.5, y=4.1, label="B2M +"), angle=-4.5, size=4) +
    geom_text(aes(x=31.5, y=4.1, label="TAP1 +"), angle=-13.5, size=4) +
    geom_text(aes(x=32.5, y=4.1, label="TAP2 +"), angle=-22.5, size=4) +
    geom_text(aes(x=33.5, y=4.1, label="HLA-A +"), angle=-31.5, size=4) +
    geom_text(aes(x=34.5, y=4.1, label="HLA-B +"), angle=-40.5, size=4) +
    geom_text(aes(x=35.5, y=4.1, label="HLA-C +"), angle=-49.5, size=4) +
    geom_text(aes(x=36.5, y=4.1, label="HLA-DPA1 +"), angle=-58.5, size=4) +
    geom_text(aes(x=37.5, y=4.1, label="HLA-DPB1 +"), angle=-67.5, size=4) +
    geom_text(aes(x=38.5, y=4.1, label="HLA-E +"), angle=-76.5, size=4) +
    geom_text(aes(x=39.5, y=4.1, label="HLA-F +"), angle=-85.5, size=4)
  plot_a5 <- plot_a4 +
    geom_text(aes(x=0, y=6, label=paste("Immunophenoscore: ",IPS[i],sep="")), angle=0,size=6,vjust=-0.5) + 
    theme(axis.title=element_blank())
  plot_a <- plot_a5 + theme(plot.margin=unit(c(0,0,0,0),"mm")) +
    geom_text(vjust=1.15,hjust=0,aes(x=25.5, y=6,label="\n\n\n\n   MHC: Antigen Processing                                 EC: Effector Cells\n   CP: Checkpoints | Immunomodulators              SC: Suppressor Cells\n\n", hjust = 0), size=4)
  
  ## 绘制图例 sample-wise (averaged) z-scores
  data_b <- data.frame (start = rep(0,23), 
                        end = rep(0.7,23), 
                        y1=seq(0,22,by=1), 
                        y2=seq(1,23,by=1),
                        z=seq(-3,3,by=6/22),
                        vcol=c(unlist(lapply(seq(-3,3,by=6/22),mapcolors))), 
                        label = LETTERS[1:23])
  data_b_ticks <- data.frame(x = rep(1.2, 7), 
                             value = seq(-3,3, by=1), 
                             y = seq(0,6, by=1)*(22/6) +0.5)
  legendtheme <- theme(plot.margin = unit(c(2,0,2,0),"inch"), 
                       panel.spacing = unit(0,"null"), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "white"), 
                       axis.text=element_blank(), 
                       axis.ticks= element_blank(), 
                       axis.title.x=element_blank())
  plot_b <- ggplot(hjust=0) + 
    geom_rect(data=data_b, 
              mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label), size=0.5,color="black", alpha=1) +
    scale_x_continuous(limits = c(0, 1.5),expand = c(0,0)) +
    scale_fill_manual(values =as.vector(data_b$vcol),guide=FALSE) +
    geom_text(data=data_b_ticks, aes(x=x, y=y, label=value),hjust="inward", size=4) +
    theme_bw() +
    legendtheme + 
    ylab("Sample-wise (averaged) z-score")
  
  ## 绘制图例 weighted z-scores
  data_c <- data.frame(start = rep(0,23),
                       end = rep(0.7,23),
                       y1=seq(0,22,by=1),
                       y2=seq(1,23,by=1),
                       z=seq(-2,2,by=4/22),
                       vcol=c(unlist(lapply(seq(-2,2,by=4/22),mapbw))), 
                       label = LETTERS[1:23])
  data_c_ticks <- data.frame(x = rep(1.2, 5), 
                             value = seq(-2,2, by=1),
                             y = seq(0,4, by=1)*(22/4) +0.5)
  plot_c <- ggplot() + 
    geom_rect(data=data_c, 
              mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label),size=0.5,color="black", alpha=1) + 
    scale_x_continuous(limits = c(0, 1.5),expand = c(0,0)) + 
    scale_fill_manual(values =as.vector(data_c$vcol),guide=FALSE) + 
    geom_text(data=data_c_ticks, aes(x=x, y=y, label=value),hjust="inward", size=4) +
    theme_bw() + 
    legendtheme + 
    ylab("Weighted z-score")
  
  ## 每一个样本保存图片
  file_name<-paste("IPS_",sample_names[i],".pdf",sep="")
  pdf(file_name, width=10, height=8)
  grid.arrange(plot_a,plot_b,plot_c, ncol=3, widths=c(0.8,0.1,0.1))
  invisible(dev.off())
  
  message(paste0("--analysis of ", tumsam[i]," done..."))
}


# 构建结果，包括各项得分以及IPS
DF <- data.frame(SAMPLE=sample_names,
                 MHC=MHC,
                 EC=EC,
                 SC=SC,
                 CP=CP,
                 AZ=AZ,
                 IPS=IPS,
                 stringsAsFactors = F)

# 输出IPS到文件
write.table(DF,file = "data/revise group/output_IPS.txt", row.names = FALSE, col.names = TRUE, quote=FALSE, sep="\t") 


DF <- read.table(file = "data/revise group/output_IPS.txt", header=TRUE, sep="\t",check.names=FALSE,stringsAsFactors = FALSE,row.names = NULL)
colnames(DF)[1] <- "ID"
DF <- merge(DF, RS, by = "ID")

#绘制箱线图
mycol <- c("#FF7F0E", "#66A61E") 

data_long <- DF %>%
  dplyr::select(ID, Group, MHC, EC, SC, CP, AZ, IPS) %>%
  pivot_longer(
    cols = c(MHC, EC, SC, CP, AZ, IPS),
    names_to = "Category",
    values_to = "Score"
  )

data_long$Category <- factor(data_long$Category, 
                             levels = c("MHC", "EC", "SC", "CP", "AZ", "IPS"))

p <- ggplot(data_long, aes(x = Group, y = Score, fill =  Group)) +
  # 使用 geom_boxplot 绘制箱线图
  geom_boxplot(alpha = 0.8,
               outlier.shape = 21,
               outlier.size = 1.5) +
  # 使用您提供的精美配色
  scale_fill_manual(values = mycol) +
  # 分面
  facet_wrap(~ Category, 
             scales = "free_y") + 
  # 设置标签
  labs(x = "Group", 
       y = "Score") +
  # 使用theme_bw并优化
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    # 分面标题设置
    strip.background = element_rect(fill = "grey95", color = "black"),  # 黑色边框
    strip.text = element_text(face = "bold", size = 10),
    # 移除所有网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 增强面板边框
    panel.border = element_rect(color = "black", linewidth = 0.8),  # 加粗边框
    # 坐标轴标题
    axis.title = element_text(face = "bold")
  )+
  stat_compare_means(
    aes(label = after_stat(p.signif)),   # 显示 p 值，如：p = 0.023
    method = "wilcox.test",
    label.y = Inf,                      # 放在顶部
    vjust = 1.2,                         # 上移一点，避免贴边
    size = 5,
    label.x = 1.5
  )
p

ggsave(p,filename = "output/revise group/F.IPS.pdf",width = 6,height = 4)



###### G.免疫景观
#-读取immune landscape------------------------------------------------
rm(list = ls())
#分组信息
TCGA_RS <- read.csv("data/revise mime/TCGA_RS.csv")
RS <- TCGA_RS[,c(2,6),drop=F]

load("data/input/immune_landscape.RData")
ABSOLUTE_scores$ID <- substr(rownames(ABSOLUTE_scores),1,12)
ABSOLUTE_scores <- ABSOLUTE_scores[!duplicated(ABSOLUTE_scores$ID),-5]
rownames(ABSOLUTE_scores) <- substr(rownames(ABSOLUTE_scores),1,12)

COAD_im <- immune_landscape[immune_landscape$TCGA.Study == "COAD",]
COAD_ab <- ABSOLUTE_scores[rownames(COAD_im),]
COAD_il<- merge(COAD_im,COAD_ab,by=0) %>%  column_to_rownames("Row.names")
rm(immune_landscape,ABSOLUTE_scores,COAD_im,COAD_ab)

rownames(RS) <- RS$ID
rownames(RS) <- substr(rownames(RS), 1, nchar(rownames(RS)) - 4)
RS <- RS[,2,drop=F]
metadata <- merge(RS,COAD_il,by = 0) %>% 
  column_to_rownames("Row.names")

ggdata <- metadata[,-c(3:4)] 

str(ggdata)
ggdata$Group[ggdata$Group %in% "High" ] = '1'
ggdata$Group[ggdata$Group %in% "Low" ] = '2'

ggdata <- apply(ggdata[,-2],2,as.numeric) %>% 
  as.data.frame()

ggdata$Group <- factor(ggdata$Group)
str(ggdata)

#-初步比较两组---------------------------------------------
library(compareGroups)
dd <- descrTable(Group~.,data = ggdata)
summary(dd)
print(dd)
pvals <- getResults(dd, "p.overall")
p.adjust(pvals, method = "BH")

export2csv(dd,file = "data/revise group/immune_landscape.csv")
gg <- read.csv( "data/revise group/immune_landscape.csv",header = T,skip = 1,check.names = F,row.names = 1)
colnames(gg) <- c("High","Low",'p')

#---内源性免疫逃逸机制-----------------------------------------------------------
colnames(ggdata)
img <- c("Group",
         "SNV.Neoantigens","Indel.Neoantigens","CTA.Score",
         
         "Intratumor.Heterogeneity",  
         "Number.of.Segments","Fraction.Altered",
         "LOH_n_seg","LOH_frac_altered",
         
         "Homologous.Recombination.Defects","Aneuploidy.Score",
         
         "TCR.Richness","TCR.Shannon")

dat <- ggdata[,img] 
rownames(dat) <- rownames(metadata)
write.csv(dat,file = "data/revise group/immunogenicity associated indicators.csv")


df <- data.frame()
for (i in 2:13) {
  m1=  mean(dat[dat$Group=="1",i][!is.na(dat[dat$Group=="1",i])])/mean(dat[,i][!is.na(dat[,i])]) ;m1=log2(m1)
  m2=  mean(dat[dat$Group=="2",i][!is.na(dat[dat$Group=="2",i])])/mean(dat[,i][!is.na(dat[,i])]) ;m2=log2(m2)
  
  df <- rbind(df,data.frame("High"=m1,"Low"=m2))
}
df <- df %>% as.matrix()

rownames(df) <- c("SNV Neoantigens","Indel Neoantigens","CTA Score",
                  
                  "Intratumor Heterogeneity",  
                  "Number of Segments","Fraction Altered",
                  "LOH_n_seg","LOH_frac_altered",
                  
                  "Homologous Recombination Defects","Aneuploidy Score",
                  
                  "TCR Richness","TCR Shannon")

colnames(df) <- c("High","Low")

library(ComplexHeatmap)
library(ggsci)
library(circlize)
library(export)
cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.2f", df[i, j]),
            x, y, gp = gpar(fontsize = 10.5))
}
top <- HeatmapAnnotation(Group=c('High','Low'),
                         border = T,
                         annotation_name_side = 'right',
                         annotation_name_gp = gpar(fontsize=10),
                         show_annotation_name = F,
                         na_col = NA,
                         show_legend = F,
                         gap = unit(1, "mm"),
                         col = list(Group=c('High'="#EFC000FF",'Low'="#197EC099")))

plo <- Heatmap(df,
               border = TRUE,
               cluster_rows = F,
               cluster_columns = F,
               heatmap_width = unit(10, "cm"),
               heatmap_height = unit(1, "cm") * nrow(df),
               col = colorRamp2(c(-1,-0.5,0.5, 1), 
                                c("#66A61E", "#98DF8A", "#FFBB78", "#FF7F0E")),
               row_split = rep(c(1, 2, 3, 4), times = c(3, 5, 2, 2)),
               gap = unit(2, "mm"),
               rect_gp = gpar(color = 'grey', lwd = 1.2, size = 0.4),
               show_column_names = T,
               row_names_side = 'right',
               row_names_gp = gpar(fontface = 'bold', color = 'black', fontsize = 12),
               heatmap_legend_param = list(
                 title = "log2(Ratio)",
                 at = c(-1, 0, 1),
                 title_position = "topcenter",
                 title_gp = gpar(fontsize = 10, fontface = "bold"),
                 legend_direction = "vertical"
               ))
plo

pdf("output/revise group/G.immunogenicity_heatmap.pdf", width = 9, height = 6)
draw(plo)
dev.off()
