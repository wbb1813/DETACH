## ---------------------------
library(ggplot2)
library(reshape)
library(circlize)
library(ggrepel)
library(ggpubr)
library('wesanderson')
library('ggsci')
library(ppcor)

## ------- Inputs -------
# etl_ctl_score_path=file.path(result_dir,'ctl_etl_score/ICB/')
# decouple_score_path=file.path(result_dir,'decouple_score/ICB/')
# outdir=file.path(result_dir,'ctl_etl_cor_decouple/ICB')
# ctl_id='CytotoxicTcell' # step 2, 4
# etl_id='plastic_exhuasted' # step 2, 4
# decouple_id='neg_genes_fdr0.01' # step 4
# ctrl_id='con_neg_genes_fdr0.01' # step 4

## signature score 
icb_skcm=readRDS(file.path(icb_sig_score_path,'icb_skcm_sig_score.rds'))

## ------ Functions -------
## Correlation of subset patients 
cor_diff=function(df,prefix){
  df_sort=df[order(df[,decouple_id],decreasing = T),]
  cor_decouple_uncouple=data.frame()
  n=0.2
  set.seed(7)
  df_high=df_sort[1:round(nrow(df_sort)*n),]
  df_other=df_sort[which(df_sort$sample_id%in%setdiff(df_sort$sample_id,df_high$sample_id)),]
  df_low=df_sort[round(nrow(df_sort)*(1-n)):nrow(df_sort),]
  df_random=df[sample(nrow(df_sort))[1:round(nrow(df_sort)*n)],]
  
  cor_couple_high=cor(df_high[,ctl_id],df_high[,etl_id])
  cor_couple_other=cor(df_other[,ctl_id],df_other[,etl_id])
  cor_couple_low=cor(df_low[,ctl_id],df_low[,etl_id])
  cor_raw=cor(df[,ctl_id],df[,etl_id])
  cor_couple_random=cor(df_random[,ctl_id],df_random[,etl_id])
  
  tmp_res=data.frame(All_samples=cor_raw,decouple_score_low=cor_couple_low,decouple_score_high=cor_couple_high,decouple_score_other=cor_couple_other,cor_random_con=cor_couple_random,high_score_cutoff=n,Cohort=prefix)
  cor_decouple_uncouple=rbind(cor_decouple_uncouple,tmp_res)
  
  return(cor_decouple_uncouple)
}

## Correlation of subset patients, high resolution 
cor_diff_high_resolution=function(df,prefix){
  df_sort=df[order(df[,decouple_id],decreasing = T),]
  cor_decouple_uncouple=data.frame()
  for (n in seq(0.01,1,0.01)){
    set.seed(7)
    df_high=df_sort[1:round(nrow(df_sort)*n),]
    df_other=df_sort[which(df_sort$sample_id%in%setdiff(df_sort$sample_id,df_high$sample_id)),]
    df_low=df_sort[round(nrow(df_sort)*(1-n)):nrow(df_sort),]
    df_random=df[sample(nrow(df_sort))[1:round(nrow(df_sort)*n)],]
    
    cor_couple_high=cor(df_high[,ctl_id],df_high[,etl_id])
    cor_couple_other=cor(df_other[,ctl_id],df_other[,etl_id])
    cor_couple_low=cor(df_low[,ctl_id],df_low[,etl_id])
    cor_raw=cor(df[,ctl_id],df[,etl_id])
    cor_couple_random=cor(df_random[,ctl_id],df_random[,etl_id])
    
    tmp_res=data.frame(All_samples=cor_raw,decouple_score_low=cor_couple_low,decouple_score_high=cor_couple_high,decouple_score_other=cor_couple_other,cor_random_con=cor_couple_random,high_score_cutoff=n,Cohort=prefix)
    cor_decouple_uncouple=rbind(cor_decouple_uncouple,tmp_res)
  }
  return(cor_decouple_uncouple)
}

# Barplot
gg_bar_icb=function(df){
  p<-ggplot(df, aes(x=Cohort, y=value, fill=variable)) +
    geom_bar(stat="identity", position=position_dodge())+theme_classic()+
    theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))+
    geom_text(aes(label=signif(value,digits = 2)), vjust=-0.35, color="black",
              position = position_dodge(0.9), size=2.5)+
    scale_fill_brewer(palette="Dark2")+xlab('')+ylab('Pearson correlation')
}

## ------- Calculate partial correlation -------
partial_outdir=file.path(outdir,'partial')
if (!dir.exists(partial_outdir)){
  dir.create(partial_outdir,recursive = T)
}

# partial correlation
partal_cor=data.frame()
for (i in names(icb_skcm)){
  tmp_meta=icb_skcm[[i]][['meta']]
  tmp_meta_pre=tmp_meta[which(tmp_meta$Timepoint=='PRE'),]
  tmp_meta_pre_sample=tmp_meta_pre$RNA_ID
  
  sig_score=icb_skcm[[i]][['sig_score']]
  sig_score=sig_score[which(sig_score$sample_id%in%tmp_meta_pre_sample),]
  sel_col=c('estimate','p.value','cohort','Group')
  # decoupling signature adjusted 

  tmp_partial=pcor.test(sig_score[,ctl_id],sig_score[,etl_id],sig_score[,decouple_id],method = 'pearson')
  tmp_partial$cohort=i
  tmp_partial$Group='Decouple sig. adjusted'
  tmp_partial=tmp_partial[,sel_col]
  
  partal_cor=rbind(partal_cor,tmp_partial)
  
  # control signature adjusted 
  tmp_ctrl=pcor.test(sig_score[,ctl_id],sig_score[,etl_id],sig_score[,ctrl_id],method = 'pearson')
  tmp_ctrl$cohort=i
  tmp_ctrl$Group='Control sig. adjusted'
  tmp_ctrl=tmp_ctrl[,sel_col]
  
  partal_cor=rbind(partal_cor,tmp_ctrl)
  
  # Unadjusted 
  raw_res=cor.test(sig_score[,ctl_id],sig_score[,etl_id],method = 'pearson')
  raw_res=data.frame(estimate=raw_res$estimate,p.value=raw_res$p.value)
  raw_res$cohort=i
  raw_res$Group='Non-adjusted'
  raw_res=raw_res[,sel_col]
  
  partal_cor=rbind(partal_cor,raw_res)
}
write.table(partal_cor,file.path(partial_outdir,'cor_partial.txt'),quote = F,sep = '\t',row.names = F)

# barplot 
partal_cor$Group=factor(partal_cor$Group,c('Decouple sig. adjusted','Non-adjusted','Control sig. adjusted'))
p=ggplot(data=partal_cor, aes(x=cohort, y=estimate, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=signif(estimate,digits = 2)), vjust=-0.5, color="black",
            position = position_dodge(0.9), size=2.5)+
  scale_fill_brewer(palette = "Dark2")+theme_classic()+xlab('')+ylab('Pearson correlation')+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust =0.5),legend.position = 'top',legend.title = element_blank())
p
ggsave(file.path(partial_outdir,'partial_vs_raw_cor.pdf'),p,width = 6,height = 4)

## ------- Correlation between high and low decoupling groups -------
hl_outdir=file.path(outdir,'high_vs_low')
if (!dir.exists(hl_outdir)){
  dir.create(hl_outdir,recursive = T)
}
## decouple group vs. coupled group
cor_hl_sum=data.frame()
for (i in names(icb_skcm)){
  tmp_meta=icb_skcm[[i]][['meta']]
  tmp_meta_pre=tmp_meta[which(tmp_meta$Timepoint=='PRE'),]
  tmp_meta_pre_sample=tmp_meta_pre$RNA_ID
  
  sig_score=icb_skcm[[i]][['sig_score']]
  sig_score=sig_score[which(sig_score$sample_id%in%tmp_meta_pre_sample),]
  # high vs low 
  cor_hl=cor_diff(df=sig_score,prefix=i)
  cor_hl_sum=rbind(cor_hl_sum,cor_hl)
}

write.table(cor_hl_sum,file.path(hl_outdir,'cor_decouple_uncouple.txt'),quote = F,sep = '\t',row.names = F)


