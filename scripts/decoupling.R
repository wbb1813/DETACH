## ---------------------------
##
## Script name: 
##
## Purpose of script: identify genes that interact with ETl and CTL score
##
## Author: Binbin Wang
##
## Date Created: 
##
## Date Modified: 
##
## Copyright (c) Binbin Wang, 2022
## Email: wangbinbintj@gmail.com
##
## ---------------------------
##
## Notes:
##   
## ---------------------------

## set working directory -----

setwd("./")     

## ---------------------------

options(scipen = 6, digits = 4) #   View outputs in non-scientific notation

## ---------------------------

## load up the packages ------
library(MASS)
library(viridis)
library(ggplot2)
library(reshape)
library(circlize)
library(ComplexHeatmap)
library(ggrepel)
library(ggpubr)
library("psych")  
library(org.Hs.eg.db)
library(clusterProfiler)
#library(ReactomePA)
#library(Seurat)
source('/Users/wangb8/Documents/Research/code/R/functions/enrichment_clusterProfiler.R')
## ---------------------------
## ------- Inputs -------
# skcm_rsem_log2_file='../data/skcm_rsem_log2.rds'
# etl_ctl_score_file='../data/TCGA_ETL_CTL_sig_score.txt'
# outdir='../results/decoupling'
# ctl_id='CytotoxicTcell'
# etl_id='fixed_exhuasted'

## load data
skcm_rsem_log2=readRDS(skcm_rsem_log2_file)

## ------- Interaction regression -------
if (!dir.exists(outdir)){
  dir.create(outdir,recursive = T)
}

## ------- Data pre-process -------
## CTL and PF score
ctl_etl_score=read.delim(etl_ctl_score_file)

tmp_score=merge(ctl_etl_score,as.data.frame(t(skcm_rsem_log2)),by.x='sample_id',by.y=0)
tmp_score=tmp_score[,c(ctl_id,etl_id,rownames(skcm_rsem_log2))]

## correlation between ETL score and proliferation score
sp <- ggscatter(tmp_score, x = etl_id, y = ctl_id,
                add = "reg.line",  # Add regression line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp = sp + stat_cor(method = "pearson") +xlab('Exhuasted T cell score') + ylab('Cytotoxic T cell score')
sp
ggsave(file.path(outdir,'cor_CTL_ETL.png'),sp,width = 4.5,height = 4.5)

## ------- decoupling ------- 
res_interaction=data.frame()
all_genes=rownames(skcm_rsem_log2)
for (j in all_genes){
  tmp_df=tmp_score[,c(ctl_id,etl_id,j)]
  colnames(tmp_df)=c('ctl','ext','gene')
  lm_model=lm(ctl~ext+gene+ext*gene,data = tmp_df)  
  tmp_sum=summary(lm_model)
  tmp_coefficients=tmp_sum$coefficients[,'t value',drop=F]
  tmp_pvalue=tmp_sum$coefficients[,'Pr(>|t|)',drop=F]
  
  if (length(tmp_coefficients)==4){
    tmp_res=rbind(tmp_coefficients,tmp_pvalue)
    rownames(tmp_res)=c('Intercept','ext','gene_t_value','ext_gene_t_value','Intercept_p_t_value','ext_p','gene_p','ext_gene_p')
    tmp_res=as.data.frame(t(tmp_res))
    tmp_res$gene_symbol=j
    
    res_interaction=rbind(res_interaction,tmp_res)
  }
}

## output results
res_interaction$ext_gene_FDR=p.adjust(res_interaction$ext_gene_p,method = 'BH')
write.table(res_interaction,file.path(outdir,'res_interaction.txt'),quote = F,sep = '\t',row.names = F)

## -------- Function annotation of significant genes ------- ##
message('Function annotation')

# MSigDb_gmt_dir='/Users/wangb8/Documents/Research/data/MSigDb/entrez/'
entrizid = data.frame(bitr(res_interaction$gene_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db",drop=T))
df_all=merge(entrizid,res_interaction,by.x='SYMBOL',by.y='gene_symbol')
df_all=df_all[order(df_all$ext_gene_t_value,decreasing = T),]
all_gene_vector=df_all$ext_gene_t_value
names(all_gene_vector)=df_all$ENTREZID

pos_genes_pvalue0.01=res_interaction[which(res_interaction$ext_gene_t_value>0&res_interaction$ext_gene_FDR<=0.01),]
neg_genes_pvalue0.01=res_interaction[which(res_interaction$ext_gene_t_value<0&res_interaction$ext_gene_FDR<=0.01),]

pos_genes_pvalue0.01_id=data.frame(bitr(pos_genes_pvalue0.01$gene_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db",drop=T))
neg_genes_pvalue0.01_id=data.frame(bitr(neg_genes_pvalue0.01$gene_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db",drop=T))

if (nrow(pos_genes_pvalue0.01_id)>=10){
  enrichAll(gene_id=pos_genes_pvalue0.01_id$ENTREZID, all_gene_vector=all_gene_vector, MSigDb_gmt_dir = MSigDb_gmt_dir, out_prefix = 'pos_genes_pvalue0.01',outdir=tmp_dir)
}

if (nrow(neg_genes_pvalue0.01_id)>=10){
  enrichAll(gene_id=neg_genes_pvalue0.01_id$ENTREZID, all_gene_vector=all_gene_vector, MSigDb_gmt_dir = MSigDb_gmt_dir, out_prefix = 'neg_genes_pvalue0.01',outdir=tmp_dir)
}

## T value rank plot
# pos and neg
df=res_interaction
df$T_value=df$ext_gene_t_value
df$Rank_of_T_value=rank(df$T_value)
df$Group='other'
df$Group[which(df$gene_symbol%in%pos_genes_pvalue0.01$gene_symbol)]='Pos'
df$Group[which(df$gene_symbol%in%neg_genes_pvalue0.01$gene_symbol)]='Neg'

if (nrow(pos_genes_pvalue0.01_id)>=10){
  pos10=df[order(df$T_value,decreasing = T),][1:10,'gene_symbol']
}else{
  pos10=pos_genes_pvalue0.01_id$gene_symbol
}
if (nrow(neg_genes_pvalue0.01)>=10){
  neg10=df[order(df$T_value,decreasing = F),][1:10,'gene_symbol']
}else{
  neg10=neg_genes_pvalue0.01$gene_symbol
}
options(ggrepel.max.overlaps = Inf)
p <- ggplot(df, aes(Rank_of_T_value, T_value)) +
  geom_point(aes(color=Group)) +
  theme_classic()
p
p=p + geom_label_repel(data = df[which(df$gene_symbol%in%c(pos10,neg10)),],
                       aes(label = gene_symbol,
                           fill = factor(Group)), color = 'white',
                       size = 3.5) +
  theme(legend.position = "none")+
  scale_fill_manual(values = c("blue","red"))+ 
  scale_color_manual(values = c("blue", "grey","red")) +
  xlab('Rank of T value')+ylab('T value')
p

ggsave(file.path(outdir,'scatter_t_value_rank.pdf'),p,width = 4,height = 3.5)

# neg only 
df=res_interaction
df$T_value=df$ext_gene_t_value
df$Rank_of_T_value=rank(df$T_value)
df$Group='other'
neg_genes_pvalue0.01=res_interaction[which(res_interaction$ext_gene_t_value<0&res_interaction$ext_gene_FDR<=0.01),]

df$Group[which(df$gene_symbol%in%neg_genes_pvalue0.01$gene_symbol)]='Neg'

neg10=df[order(df$T_value,decreasing = F),][1:10,'gene_symbol']

options(ggrepel.max.overlaps = Inf)
p <- ggplot(df, aes(Rank_of_T_value, T_value)) +
  geom_point(aes(color=Group)) +
  theme_classic()
p
p=p + geom_label_repel(data = df[which(df$gene_symbol%in%c(neg10)),],
                       aes(label = gene_symbol,
                           fill = factor(Group)), color = 'white',
                       size = 3.5) +
  theme(legend.position = "none")+
  scale_fill_manual(values = c("blue","red"))+ 
  scale_color_manual(values = c("blue", "grey","red")) +
  xlab('Rank of T value')+ylab('T value')
p

ggsave(file.path(outdir,'scatter_t_value_rank_neg.pdf'),p,width = 3.5,height = 3)

