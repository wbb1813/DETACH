## ---------------------------
library(ggplot2)
library(reshape)
library(circlize)
library(ggrepel)
library(ggpubr)
library('wesanderson')
library('ggsci')

## ------- Inputs (GSEA results) -------
# tcga_cor_hl_file='../../results/bulk/ctl_etl_cor_decouple/TCGA/gsea/high_vs_low/cor_decouple_uncouple.txt' 
# tcga_cor_partial_file='../../results/bulk/ctl_etl_cor_decouple/TCGA/gsea/cor_sum.txt' 
# icb_cor_hl_file='../../results/bulk/ctl_etl_cor_decouple/ICB/gsea/high_vs_low/cor_decouple_uncouple.txt' 
# icb_cor_partial_file='../../results/bulk/ctl_etl_cor_decouple/ICB/gsea/cor_sum.txt'
# 
# outdir='../../results/bulk/ctl_etl_cor_decouple/TCGA_ICB'

## ------- High vs low -------
tcga_cor_hl=read.delim(tcga_cor_hl_file)
icb_cor_hl=read.delim(icb_cor_hl_file)

tcga_cor_hl$Cohort='TCGA'
tcga_cor_hl$Group='Training'
icb_cor_hl$Group='Testing'

df_cor_hl=rbind(tcga_cor_hl,icb_cor_hl)

if (!dir.exists(outdir)){
  dir.create(outdir,recursive = T)
}

for (cutoff in unique(df_cor_hl$high_score_cutoff)){
  tmp_df=df_cor_hl[which(df_cor_hl$high_score_cutoff==cutoff),]
  tmp_df=melt(tmp_df)
  tmp_df=tmp_df[which(tmp_df$variable%in%c('decouple_score_low','decouple_score_high')),]
  tmp_df$variable=gsub('decouple_score_','',tmp_df$variable)
  tmp_df$variable=factor(tmp_df$variable,c('high','low'))
  tmp_df$Group=factor(tmp_df$Group,levels = c('Training','Testing'))
  
  p=ggplot(data=tmp_df, aes(x=Cohort, y=value, fill=variable)) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=signif(value,digits = 2)), vjust=1.6, color="black",
              position = position_dodge(0.9), size=2.5)+
    scale_fill_brewer(palette = "Dark2")+theme_classic()+xlab('')+ylab('Pearson correlation')+
    theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5),legend.position = 'top',legend.title = element_blank())+
    facet_grid(.~Group,scales = 'free_x',space = "free_x")
  set_palette(p, "jco")
  ggsave(file.path(outdir,paste0('top_',cutoff,'_decouple_score_high_vs_low_cor.pdf')),width = 6,height = 4)
  
  ## High vs others
  tmp_df=df_cor_hl[which(df_cor_hl$high_score_cutoff==cutoff),]
  tmp_df=melt(tmp_df)
  tmp_df=tmp_df[which(tmp_df$variable%in%c('decouple_score_other','decouple_score_high')),]
  tmp_df$variable=gsub('decouple_score_','',tmp_df$variable)
  tmp_df$variable=factor(tmp_df$variable,c('other','high'))
  tmp_df$Group=factor(tmp_df$Group,levels = c('Training','Testing'))
  
  p=ggplot(data=tmp_df, aes(x=Cohort, y=value, fill=variable)) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=signif(value,digits = 2)), vjust=1.6, color="black",
              position = position_dodge(0.9), size=2.5)+
    scale_fill_brewer(palette = "Dark2")+theme_classic()+xlab('')+ylab('Pearson correlation')+
    theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5),legend.position = 'top',legend.title = element_blank())+
    facet_grid(.~Group,scales = 'free_x',space = "free_x")
  set_palette(p, "jco")
  ggsave(file.path(outdir,paste0('top_',cutoff,'_decouple_score_high_vs_other_cor.pdf')),width = 6,height = 4)
  
  ## High vs All
  tmp_df=df_cor_hl[which(df_cor_hl$high_score_cutoff==cutoff),]
  tmp_df=melt(tmp_df)
  tmp_df=tmp_df[which(tmp_df$variable%in%c('All_samples','decouple_score_high')),]
  tmp_df$variable=gsub('decouple_score_','',tmp_df$variable)
  tmp_df$variable[which(tmp_df$variable=='All_samples')]='all'
  tmp_df$variable=factor(tmp_df$variable,c('high','all'))
  tmp_df$Group=factor(tmp_df$Group,levels = c('Training','Testing'))

  p=ggplot(data=tmp_df, aes(x=Cohort, y=value, fill=variable)) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=signif(value,digits = 2)), vjust=1.6, color="black",
              position = position_dodge(0.9), size=2.5)+
    scale_fill_brewer(palette = "Dark2")+theme_classic()+xlab('')+ylab('Pearson correlation')+
    theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5),legend.position = 'top',legend.title = element_blank())+
    facet_grid(.~Group,scales = 'free_x',space = "free_x")
  set_palette(p, "jco")
  ggsave(file.path(outdir,paste0('top_',cutoff,'_decouple_score_high_vs_all_cor.pdf')),width = 6,height = 4)
}

## ------- High vs low, high resolution -------
tcga_cor_hl_high_res=read.delim(tcga_cor_hl_high_resolution_file)
icb_cor_hl_high_res=read.delim(icb_cor_hl_high_resolution_file)

tcga_cor_hl_high_res$Cohort='TCGA'
tcga_cor_hl_high_res$Group='Training'
icb_cor_hl_high_res$Group='Testing'

df_cor_hl_high_res=rbind(tcga_cor_hl_high_res,icb_cor_hl_high_res)

# start percentage 
start_pct=0.1
df=df_cor_hl_high_res[which(df_cor_hl_high_res$high_score_cutoff>start_pct),]
df$decouple_score_high=signif(df$decouple_score_high,2)

for (cohort in unique(df$Cohort)){
  tmp_df=df[which(df$Cohort==cohort),]
  tmp_df=tmp_df[which(tmp_df$decouple_score_high==min(tmp_df$decouple_score_high)),]
  print(cohort)
  print(tmp_df$high_score_cutoff)
}

start_pct=0.1
df=df_cor_hl_high_res[which(df_cor_hl_high_res$high_score_cutoff>start_pct),]
df=df[which(df$high_score_cutoff%in%as.character(seq(0.1,1,0.05))),]
df$decouple_score_high=signif(df$decouple_score_high,2)

for (cohort in unique(df$Cohort)){
  tmp_df=df[which(df$Cohort==cohort),]
  tmp_df=tmp_df[which(tmp_df$decouple_score_high==min(tmp_df$decouple_score_high)),]
  print(cohort)
  print(tmp_df$high_score_cutoff)
}
# 0.16 0.17 gide ok
# 0.20 riaz ok
# 0.17 0.18 Van ok
# 0.18 0.19 puch ok
# 0.18 0.17Liu 
# 0.19 0.25 TCGA 
sum(c(0.19,0.16,0.18,0.19,0.15,0.15,0.27))/7
df_high_res_sel=df_cor_hl_high_res[which((df_cor_hl_high_res$Cohort=='TCGA'&df_cor_hl_high_res$high_score_cutoff==0.25)|
                                          (df_cor_hl_high_res$Cohort=='Liu_2019'&df_cor_hl_high_res$high_score_cutoff==0.15)|
                                          (df_cor_hl_high_res$Cohort=='PUCH'&df_cor_hl_high_res$high_score_cutoff==0.15)|
                                          (df_cor_hl_high_res$Cohort=='VanAllen_2015'&df_cor_hl_high_res$high_score_cutoff==0.15)|
                                          (df_cor_hl_high_res$Cohort=='Riaz_2017'&df_cor_hl_high_res$high_score_cutoff==0.15)|
                                          (df_cor_hl_high_res$Cohort=='Gide_2019'&df_cor_hl_high_res$high_score_cutoff==0.15)),]

## High vs All
tmp_df=melt(df_high_res_sel)
tmp_df=tmp_df[which(tmp_df$variable%in%c('All_samples','decouple_score_high')),]
tmp_df$variable=gsub('decouple_score_','',tmp_df$variable)
tmp_df$variable[which(tmp_df$variable=='All_samples')]='all'
tmp_df$variable=factor(tmp_df$variable,c('high','all'))
tmp_df$Group=factor(tmp_df$Group,levels = c('Training','Testing'))

p=ggplot(data=tmp_df, aes(x=Cohort, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=signif(value,digits = 2)), vjust=1.6, color="black",
            position = position_dodge(0.9), size=2.5)+
  scale_fill_brewer(palette = "Dark2")+theme_classic()+xlab('')+ylab('Pearson correlation')+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5),legend.position = 'top',legend.title = element_blank())+
  facet_grid(.~Group,scales = 'free_x',space = "free_x")
set_palette(p, "jco")
ggsave(file.path(outdir,paste0('high_resolution_individual_cutoff','_decouple_score_high_vs_all_cor.pdf')),width = 6,height = 4)


## ------- Partial -------
## TCGA cor results
tcga_cor_partial=read.delim(tcga_cor_partial_file)
icb_cor_partial=read.delim(icb_cor_partial_file)
tcga_cor_partial$cohort='TCGA'

tcga_cor_partial$Group2='Training'
icb_cor_partial$Group2='Testing'

col_sel=c("estimate","p.value","cohort","Group","Group2" )
df_partial=rbind(tcga_cor_partial[,col_sel],icb_cor_partial[,col_sel])
df_partial$Group=factor(df_partial$Group,levels = c("Decouple sig. adjusted","Control sig. adjusted","Non-adjusted"))
df_partial$Group2=factor(df_partial$Group2,levels = c("Training","Testing"))

p=ggplot(data=df_partial, aes(x=cohort, y=estimate, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=signif(estimate,digits = 2)), vjust=-0.35, color="black",
            position = position_dodge(0.9), size=2.5)+
  scale_fill_brewer(palette = "Dark2")+theme_classic()+xlab('')+ylab('Correlation')+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.5,vjust = 0.5),legend.position = 'top',legend.title = element_blank())+
  facet_grid(.~Group2,scales = 'free_x',space = "free_x")
p
ggsave(file.path(outdir,paste0('partial_tcga_icb.pdf')),width = 7,height = 3.95)

