library(ggplot2)
library(ggpubr)

## ------- Inputs -------
tcga_etl_ctl_score=read.delim(tcga_etl_ctl_score_file)
tcga_other_sig_score=read.delim(tcga_other_sig_score_file)
tcga_decouple_score=read.delim(decouple_score_file)

icb_skcm=readRDS(icb_sig_score_file)

## TCGA clinical data
tcga_clinic=read.delim(tcga_clinic_file)
tcga_clinic=tcga_clinic[!duplicated(tcga_clinic$case_submitter_id),]

tmp_clinic=tcga_clinic[,c("case_submitter_id","project_id")]
tmp_clinic$project_id=gsub('TCGA-','',tmp_clinic$project_id,)

## TCGA immune cell abundance
skcm_frac=readRDS(tcga_skcm_frac_file)

## COLDFACS immune cell abundance of ICB cohorts
cell_frac=read.delim(cell_frac_file)
cell_frac$id=rownames(cell_frac)
cell_frac$id=sapply(cell_frac$id,function(x)substr(x,start = 1,stop = 12))
cell_frac=cell_frac[!duplicated(cell_frac$id),]
rownames(cell_frac)=cell_frac$id
cell_frac=cell_frac[,-ncol(cell_frac)]

## Other ICB signature score
icb_sig_score=readRDS(icb_response_score)

## TCGA pathology slide TIL data
til_df=read.delim(til_file)

hot_cate=c('Brisk Diffuse','Brisk Band-like')
cold_cate=c("None","Non-Brisk Focal","Non-Brisk Multifocal")

til_hc=til_df[which(til_df$Global_Pattern%in%c(hot_cate,cold_cate)),]
til_hc$Group='cold'
til_hc$Group[which(til_hc$Global_Pattern==hot_cate)]='hot'
rownames(til_hc)=til_hc$ParticipantBarcode; til_hc=til_hc[,-1]

# ## reverse TIDE score, high TIDE score related to worth ICB response
# for (i in names(icb_sig_score)){
#   icb_sig_score[[i]]$TIDE=(-icb_sig_score[[i]]$TIDE)
# }

if (!dir.exists(outdir)){
  dir.create(outdir)
}

## ------- Correlation with reactive T cell score -------
sel_sig=c('CytotoxicTcell','IMPRES','TIDE','CD274','Stroma_EMT','CD8_T_effector_POPLAR','TGFB','neg_genes_fdr0.01')

cal_cor=function(df_score,cohort){
  tmp_df=df_score[,c('CD8_NeoTCR8',sel_sig)]
  cor_df=as.data.frame(cor(tmp_df,method = 'pearson'))
  cor_df$group=cohort
  cor_df=cor_df[which(rownames(cor_df)!='CD8_NeoTCR8'),]
  cor_df=cor_df[,c('CD8_NeoTCR8','group')]
  cor_df$sig_name=rownames(cor_df)
  return(cor_df)
}

## TCGA
tcga_df=merge(tcga_etl_ctl_score,tcga_other_sig_score,by='sample_id')
tcga_df=merge(tcga_df,tcga_decouple_score,by='sample_id')
tcga_df$Group='TCGA SKCM'
tcga_df$TIDE=-tcga_df$TIDE

tcga_cor=cal_cor(df_score=tcga_df,cohort='TCGA_SKCM')

## ICB
cor_res=tcga_cor
for (i in names(icb_skcm)){
  tmp_score1=icb_skcm[[i]]$sig_score
  tmp_score2=icb_sig_score[[i]]
  tmp_score=merge(tmp_score1,tmp_score2,by='sample_id')
  
  tmp_cor=cal_cor(df_score=tmp_score,cohort=i)
  tmp_cor$sig_name=rownames(tmp_cor)
  cor_res=rbind(cor_res,tmp_cor)
}

cor_res$sig_name[which(cor_res$sig_name=='CytotoxicTcell')]='CTL'
cor_res$sig_name[which(cor_res$sig_name=='Stroma_EMT')]='Stroma EMT'
cor_res$sig_name[which(cor_res$sig_name=='CD8_T_effector_POPLAR')]='CD8 T effector'
cor_res$sig_name[which(cor_res$sig_name=='neg_genes_fdr0.01')]='Decoupling score'
cor_res$sig_name=factor(cor_res$sig_name,c('Decoupling score','CTL','IMPRES','TIDE','CD274','Stroma EMT','CD8 T effector','TGFB'))

## Barplot
p=ggplot(data = cor_res, aes(x = sig_name, y = CD8_NeoTCR8)) +
  geom_bar(stat="identity",color="black", fill="#22908C")+theme_classic()+
  xlab('')+ylab('Pearson correlation')+theme(legend.position = 'none',axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5))+
  facet_wrap(vars(group),ncol = 3)
p
ggsave(file.path(outdir,'react_icb_sig_box.pdf'),p,height = 6,width = 9)

# Change the colors manually
p <- ggplot(data=cor_res, aes(x=group, y=CD8_NeoTCR8, fill=sig_name)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()
# Use custom colors
p + scale_fill_manual(values=c('#999999','#E69F00'))
# Use brewer color palettes
p + scale_fill_brewer(palette="Blues")
## ------- Association between decoupling and reactive T cell score -------
cor_scatter=function(dat, prefix,nrow, width = 3.5,height = 3){
  sp <- ggscatter(dat, x = react_id, y = decouple_id,
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE , # Add confidence interval
  )
  # Add correlation coefficient
  sp = sp + stat_cor(method = "pearson") +xlab('Reactive T cell score')+ylab('Decoupling score')+
    facet_wrap(~ Group,nrow = nrow)
    
  sp
  ggsave(file.path(outdir,paste0(prefix,'_react_decouple_scatter.pdf')),sp,width = width,height = height)
  
}

# TCGA skcm
tcga_df=merge(tcga_etl_ctl_score,tcga_other_sig_score,by='sample_id')
tcga_df=merge(tcga_df,tcga_decouple_score,by='sample_id')

tcga_df$Group='TCGA SKCM'

cor_scatter(dat = tcga_df,prefix = 'tcga_skcm',nrow=1)

# ICB cohorts
icb_score_df=data.frame()
for (i in names(icb_skcm)){
  tmp_df=icb_skcm[[i]]$sig_score
  tmp_df$Group=i
  icb_score_df=rbind(icb_score_df,tmp_df)
}

cor_scatter(dat = icb_score_df,prefix = 'comb',nrow=1, width = 15,height = 3)

# combination
com_col=c(react_id,decouple_id,'Group')
comb_df=rbind(tcga_df[,com_col],icb_score_df[,com_col])
comb_df$Group=factor(comb_df$Group,unique(comb_df$Group))

cor_scatter(dat = comb_df,prefix = 'comb',nrow=3, width = 6,height = 9)

## ------- Decoupling score difference between hot and code tumor, TCGA -------
til_df$ParticipantBarcode=gsub('-','.',til_df$ParticipantBarcode)
df_tcga_score=merge(tcga_df,til_df,by.x = 'sample_id',by.y = 'ParticipantBarcode')

## skcm
# all groups
df_skcm=df_tcga_score[which(df_tcga_score$Study=='SKCM'),]
df_skcm$Global_Pattern=factor(df_skcm$Global_Pattern,levels = c("Brisk Diffuse","Brisk Band-like","Non-Brisk Multifocal","Non-Brisk Focal","Indeterminate","None"))

colnames(df_skcm)[which(colnames(df_skcm)==ctl_id)]='ctl'
colnames(df_skcm)[which(colnames(df_skcm)==etl_id)]='etl'
colnames(df_skcm)[which(colnames(df_skcm)==decouple_id)]='decouple'

## 0-1 scale score
df_skcm$decouple=(df_skcm$decouple - min(df_skcm$decouple)) / (max(df_skcm$decouple) - min(df_skcm$decouple))
df_skcm$ctl=(df_skcm$ctl - min(df_skcm$ctl)) / (max(df_skcm$ctl) - min(df_skcm$ctl))
df_skcm$etl=(df_skcm$etl - min(df_skcm$etl)) / (max(df_skcm$etl) - min(df_skcm$etl))

p <- ggboxplot(df_skcm, x = "Global_Pattern", y = "decouple",
               color = "Global_Pattern", palette = "jco",
               add = "jitter")+
  stat_compare_means(method = "anova")+
  xlab('')+ylab('Decoupling score')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),legend.position = 'none')
p

ggsave(file.path(outdir,'pathology_slides_decouple_score_all_group.pdf'),height = 5,width = 6.5)

# filter Indeterminate
df_tmp=df_skcm[which(df_skcm$Global_Pattern!='Indeterminate'),]
df_tmp$Global_Pattern=factor(df_tmp$Global_Pattern,levels = c("Brisk Diffuse","Brisk Band-like","Non-Brisk Multifocal","Non-Brisk Focal","None"))

my_comparisons <- list(c("Brisk Diffuse", "Brisk Band-like"), c("Brisk Diffuse", "Non-Brisk Multifocal"), c("Brisk Diffuse", "Non-Brisk Focal"),c("Brisk Diffuse", "None"),
                       c("Brisk Band-like","Non-Brisk Multifocal"),c("Brisk Band-like","Non-Brisk Focal"),c("Brisk Band-like","None"),
                       c("Non-Brisk Multifocal","Non-Brisk Focal"),c("Non-Brisk Multifocal","None"),c("Non-Brisk Focal","None"))

# decoupling score
p <- ggboxplot(df_tmp, x = "Global_Pattern", y = "decouple",
               color = "Global_Pattern", palette = "jco",
               add = "jitter")+
  #stat_compare_means(method = "anova")+
  xlab('')+ylab('Decoupling score')+
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),legend.position = 'none')
p

ggsave(file.path(outdir,'pathology_slides_decouple_score_sub_group.pdf'),height = 7,width = 4.5)

# hot and cold
df_tmp$Group='Hot'
df_tmp$Group[which(df_tmp$Global_Pattern%in%cold_cate)]='Cold'
df_tmp$Group=factor(df_tmp$Group,levels = c('Hot','Cold'))

# Decoupling score
p <- ggboxplot(df_tmp, x = "Group", y = "decouple",
               color = "Group", palette = "Set1",
               add = "jitter")+
  stat_compare_means(method = "wilcox.test")+
  xlab('')+ylab('Decoupling score')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),legend.position = 'none')
p

ggsave(file.path(outdir,'pathology_slides_decouple_score_hot_cold.pdf'),height = 4,width = 2.5)

## ------- Correlation with immune cell abundance, TCGA -------
df=merge(skcm_frac,tcga_decouple_score,by = 'sample_id')

df=df[,setdiff(c(colnames(skcm_frac),decouple_id),'sample_id')]

cor_sum=data.frame()
for (i in setdiff(colnames(skcm_frac),'sample_id')){
  tmp_df=df[,c(i,decouple_id)]
  tmp_cor=cor.test(tmp_df[,i],tmp_df[,decouple_id],method = 'pearson')
  tmp_res=data.frame(cor=tmp_cor$estimate,pvalue=tmp_cor$p.value,Cell_type=i)
  cor_sum=rbind(cor_sum,tmp_res)
}

cor_sum$FDR=p.adjust(cor_sum$pvalue)

cor_sum$sig='ns'
cor_sum$sig[which(cor_sum$FDR<=0.05)]='*'
cor_sum$sig[which(cor_sum$FDR<=0.01)]='**'
cor_sum$sig[which(cor_sum$FDR<=0.001)]='***'
cor_sum=cor_sum[order(cor_sum$cor,decreasing = T),]
cor_sum$Cell_type=factor(cor_sum$Cell_type,cor_sum$Cell_type)
cor_sum$Color='pos'
cor_sum$Color[which(cor_sum$cor<0)]='neg'


p=ggplot(data = cor_sum, aes(x = Cell_type, y = cor,fill=Color)) +
  geom_bar(stat="identity")+theme_classic()+
  geom_text(aes(label=sig,y = max(cor_sum$cor)), vjust=-0.35, color="black",
            position = position_dodge(0.9), size=3.5)+
  # scale_fill_manual(values = c("#0072B2", "#D55E00"))+
  scale_fill_manual(values = c("#00AFBB", "#FC4E07"))+
  xlab('')+ylab('Pearson correlation')+theme(legend.position = 'none',axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p
ggsave(file.path(outdir,'skcm_frac_decouple_cor.pdf'),p,height = 3.5,width = 4.5)


