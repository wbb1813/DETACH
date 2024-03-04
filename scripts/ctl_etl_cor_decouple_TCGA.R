library(ggplot2)
library(reshape)
library(circlize)
library(ggrepel)
library(ggpubr)
library('wesanderson')
library('ggsci')
library(ppcor)

## ------- Inputs -------
# etl_ctl_score_path='../results/ctl_plastic/ctl_etl_score/TCGA/'
# decouple_score_path='../results/ctl_plastic/decouple_score/TCGA/'
etl_ctl_score=read.delim(file.path(etl_ctl_score_path,'gsea_sig_score.txt'))
decouple_score=read.delim(file.path(decouple_score_path,'gsea_sig_score.txt'))

df_score=merge(etl_ctl_score,decouple_score,by='sample_id')

## ------- ETL and CTL partial correlation -------
# partial correlation 
partial_outdir=file.path(outdir,'partial')
if (!dir.exists(partial_outdir)){
  dir.create(partial_outdir,recursive = T)
}

tmp_res=pcor.test(df_score[,ctl_id],df_score[,etl_id],df_score[,decouple_id],method = 'pearson')
tmp_res$Group='Decouple sig. adjusted'

tmp_ctrl_res=pcor.test(df_score[,ctl_id],df_score[,etl_id],df_score[,ctrl_id],method = 'pearson')
tmp_ctrl_res$Group='Control sig. adjusted'

tmp_raw_res=cor.test(df_score[,ctl_id],df_score[,etl_id],method = 'pearson')
tmp_raw_res=data.frame(estimate=tmp_raw_res$estimate,p.value=tmp_raw_res$p.value)
tmp_raw_res$Group='Non-adjusted'

com_col=intersect(colnames(tmp_res),colnames(tmp_raw_res))

df_res=rbind(tmp_res[,com_col],tmp_ctrl_res[,com_col],tmp_raw_res[,com_col])
df_res$Group=factor(df_res$Group,levels = c('Decouple sig. adjusted','Non-adjusted','Control sig. adjusted'))

write.table(df_res,file.path(partial_outdir,'partial_cor.txt'),quote = F,row.names = F,sep = '\t')

## ------- High decoupling group vs low group -------
hl_outdir=file.path(outdir,'high_vs_low')
if (!dir.exists(hl_outdir)){
  dir.create(hl_outdir,recursive = T)
}

df_sort=df_score[order(df_score[,decouple_id],decreasing = T),]

n=0.2
cor_decouple_hl=data.frame()
df_high=df_sort[1:round(nrow(df_sort)*n),]
df_other=df_sort[which(df_sort$sample_id%in%setdiff(df_sort$sample_id,df_high$sample_id)),]
df_low=df_sort[round(nrow(df_sort)*(1-n)):nrow(df_sort),]
set.seed(7)
df_random=df_sort[sample(nrow(df_sort))[1:round(nrow(df_sort)*n)],]

cor_couple_high=cor(df_high[,ctl_id],df_high[,etl_id])
cor_couple_other=cor(df_other[,ctl_id],df_other[,etl_id])
cor_couple_low=cor(df_low[,ctl_id],df_low[,etl_id])
cor_raw=cor(df_score[,ctl_id],df_score[,etl_id])
cor_random=cor(df_random[,ctl_id],df_random[,etl_id])

tmp_res=data.frame(All_samples=cor_raw,decouple_score_low=cor_couple_low,decouple_score_high=cor_couple_high,decouple_score_other=cor_couple_other,cor_random_con=cor_random,high_score_cutoff=n)
cor_decouple_hl=rbind(cor_decouple_hl,tmp_res)

write.table(cor_decouple_hl,file.path(hl_outdir,'cor_decouple_hl.txt'),row.names = F,sep = '\t',quote = F)


