options(java.parameters = "- Xmx4096m")
library(ggplot2)
library(reshape)
library(circlize)
library(ggrepel)
library(ggpubr)
library("psych")  
library(xlsx)
library(pROC)
library(viridis)
library(caret)
library(wesanderson)
library(cutpointr)
library(psych)
# source('/Users/wangb8/Documents/Research/projects/cold2hot/scripts/functions/general.R')

## ------- Inputs -------
# etl_ctl_score_path=file.path(result_dir,'ctl_etl_score/ICB/')
# decouple_score_path=file.path(result_dir,'decouple_score/ICB/')
# outdir=file.path(result_dir,'ICB_response')
# decouple_id='neg_genes_fdr0.01' # step 4
# ctrl_id='con_neg_genes_fdr0.01' # step 4
## ------- meta data of ICB cohort -------
load(icb_meta)
icb_skcm=readRDS(file.path(icb_sig_score_path,'icb_skcm_sig_score.rds'))

## Curate Gide meta data
meta_gide_pd1$Treatment='aPD1'
meta_gide_ctla4$Treatment='aPD1+CTLA4'
meta_gide=rbind(meta_gide_pd1,meta_gide_ctla4)

col_order=colnames(icb_skcm$Gide_2019$meta)
tmp3=icb_skcm$Gide_2019$meta
colnames(tmp3)[which(colnames(tmp3)=='Treatment')]='old_treat'
tmp4=merge(tmp3,meta_gide[,c('Treatment','Run')],by.x='RNA_ID',by.y='Run')
tmp4=tmp4[,col_order]
icb_skcm$Gide_2019$meta=tmp4

## -------- signature score --------
## Other ICB signature score
icb_sig_score=readRDS(icb_response_score)

## ------- Functions -------
## odd ratio with optimized cutoff by maxmized the sum of spe and sen
response_odd_ratio=function(sig_score,score_id){
  # df=merge(sig_score,meta_df,by.x = 'sample_id',by.y = 'Run')
  df=sig_score
  colnames(df)[which(colnames(df)==score_id)]='score_id'
  
  if (length(unique(df$Response))==2){
    df$Response=as.factor(df$Response)
    
    opt_cut <- cutpointr(df, score_id, Response, direction = ">=", pos_class = "1",
                         neg_class = "0", method = maximize_metric, metric = sum_sens_spec)
    
    cp_sum=summary(opt_cut)
    tmp_matrix=cp_sum$confusion_matrix[[1]]
    
    odd_ratio=(tmp_matrix[,'tp']+tmp_matrix[,'tn'])/(tmp_matrix[,'fn']+tmp_matrix[,'fp'])
    
    tmp_matrix$odd_ratio=odd_ratio
    tmp_matrix$AUC=opt_cut$AUC
    
    return(tmp_matrix)
  }else{
    tmp_res=data.frame(cutpoint=NA,tp=NA, fn=NA, fp=NA, tn=NA, odd_ratio=NA,AUC=NA) 
    return(tmp_res)
  }
}

## odd ratio for multiple cutoffs
mult_odd_ratio=function(sig_score,meta_df,score_id,group,cohort){
  odd_sum=data.frame()
  j=0.2
  df=merge(sig_score,meta_df,by.x = 'sample_id',by.y = 'RNA_ID')
  df_sort=df[order(df[,decouple_id],decreasing = T),]
  df_high=df_sort[1:round(nrow(df_sort)*j),]
  df_other=df_sort[which(df_sort$sample_id%in%setdiff(df_sort$sample_id,df_high$sample_id)),]
  df_low=df_sort[round(nrow(df_sort)*(1-j)):nrow(df_sort),]
  
  ## low group
  res_low=response_odd_ratio(sig_score=df_low,score_id)
  res_low$Decouple='Low'
  res_low$Cutoff_high=j
  odd_sum=rbind(odd_sum,res_low)
  
  ## high group
  res_high=response_odd_ratio(sig_score=df_high,score_id)
  res_high$Decouple='High'
  res_high$Cutoff_high=j
  odd_sum=rbind(odd_sum,res_high)
  
  ## other group
  res_other=response_odd_ratio(sig_score=df_other,score_id)
  res_other$Decouple='Other'
  res_other$Cutoff_high=j
  odd_sum=rbind(odd_sum,res_other)
  
  ## all group
  res_all=response_odd_ratio(sig_score=df,score_id)
  res_all$Decouple='All'
  res_all$Cutoff_high=j
  odd_sum=rbind(odd_sum,res_all)
  
  odd_sum$Group=group
  odd_sum$Cohort=cohort
  return(odd_sum)
}

## barplot for AUC and odd ratio
# percentage
bar_plot=function(df,decouple_group,sub_dir='high_vs_low'){
  sub_dir=file.path(outdir,sub_dir)
  if (!dir.exists(file.path(sub_dir))){
    dir.create(file.path(sub_dir),recursive = T)
  }

  df$Group[which(df$Group==ctl_id)]='CTL'
  df$Group[which(df$Group=='T_effector_IFNG_response')]='CD8 T effector'
  df$Group[which(df$Group=='Stroma_EMT')]='Stroma EMT'
  for (cutoff in unique(df$Cutoff_high)){
    # cohort_sel=c("Gide PD1","Gide CTLA4","Riaz Pre PD1","Liu PD1","Liu CTLA4")
    # all_groups=c("CTL","IMPRES","TIDE","CD274","CXCL9","Proliferation","Stroma EMT","CD8 T effector","CD8_T_effector_POPLAR","TGFB","Myeloid_inflammation" )
    all_groups=c("CTL","IMPRES","TIDE","CD274",'Stroma EMT',"CD8 T effector","TGFB" )
    ## AUC and odd ratio
    df_hl=df[which(df$Cutoff_high%in%cutoff&df$Group=='CTL'),]
    
    df_hl=df_hl[which(df_hl$Decouple%in%decouple_group),]
    df_hl=na.omit(df_hl)
    tmp_count=as.data.frame(table(df_hl$Cohort))
    tmp_count=tmp_count[which(tmp_count$Freq==2),]
    cohort_sel=tmp_count$Var1
    df_hl=df_hl[which(df_hl$Cohort%in%cohort_sel),]
    df_hl$Decouple=factor(df_hl$Decouple,decouple_group)
    
    # AUC
    p <- ggplot(data=df_hl, aes(x=Cohort, y=AUC, fill=Decouple)) +
      geom_bar(stat="identity", color="black", position=position_dodge())+
      theme_classic()+ylab('AUC')+xlab('')+
      theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,colour = 'black'),legend.position = 'top')+
      geom_text(aes(label=signif(AUC,digits = 2)), vjust=-0.5, color="black",
                position = position_dodge(0.9), size=2)
    # Use brewer color palettes
    p = p + scale_fill_brewer(palette="Set1")
    p
    
    ggsave(file.path(sub_dir,paste0(cutoff,'_AUC_CTL_bulk.pdf')),p,width =6,height = 4)
    
    # Odd ratio
    p <- ggplot(data=df_hl, aes(x=Cohort, y=odd_ratio, fill=Decouple)) +
      geom_bar(stat="identity", color="black", position=position_dodge())+
      theme_classic()+ylab('Odd ratio')+xlab('')+
      theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,colour = 'black'),legend.position = 'top')+
      geom_text(aes(label=signif(odd_ratio,digits = 2)), vjust=-0.5, color="black",
                position = position_dodge(0.9), size=2)
    # Use brewer color palettes
    p = p + scale_fill_brewer(palette="Set1")
    p
    
    ggsave(file.path(sub_dir,paste0(cutoff,'_odd_ratio_CTL_bulk.pdf')),p,width = 6,height = 4)
    
    ## compare with impress and tide
    df_all=df[which(df$Cutoff_high==cutoff&df$Group%in%all_groups),]
    df_all=df_all[which(df_all$Decouple=='High'),]
    df_all=na.omit(df_all)
    tmp_count=as.data.frame(table(df_all$Cohort))
    if (nrow(tmp_count)!=0){
      tmp_count=tmp_count[which(tmp_count$Freq==length(all_groups)),]
      cohort_sel=as.character(tmp_count$Var1)
      df_all=df_all[which(df_all$Cohort%in%cohort_sel),]
      df_all$Group=factor(df_all$Group,all_groups)
      # AUC
      p <- ggplot(data=df_all, aes(x=Cohort, y=AUC, fill=Group)) +
        geom_bar(stat="identity", color="black", position=position_dodge())+
        theme_classic()+ylab('AUC')+xlab('')+
        theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,colour = 'black'))+
        geom_text(aes(label=signif(AUC,digits = 2)), vjust=-0.5, color="black",
                  position = position_dodge(0.9), size=2)
      # Use brewer color palettes
      #p=p + scale_fill_viridis(discrete=TRUE, option="viridis")
      # p=p + scale_fill_manual(values = wes_palette("Darjeeling1"))
      colors <- c("#d62728", "#ff7f0e", "#2ca02c","#1f77b4" , "#9467bd",
                  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                  "#aec7e8", "#ffbb78", "#98df8a")
      
      p=p + scale_fill_manual(values=colors)
      p
      ggsave(file.path(sub_dir,paste0(cutoff,'_AUC_all_Riaz_pre.pdf')),p,width = 10,height = 4)
      
      # odd ratio
      p <- ggplot(data=df_all, aes(x=Cohort, y=odd_ratio, fill=Group)) +
        geom_bar(stat="identity", color="black", position=position_dodge())+
        theme_classic()+ylab('Odd ratio')+xlab('')+
        theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,colour = 'black'))+
        geom_text(aes(label=signif(odd_ratio,digits = 2)), vjust=-0.5, color="black",
                  position = position_dodge(0.9), size=2)
      # Use brewer color palettes
      #p=p + scale_fill_viridis(discrete=TRUE, option="viridis")
      #p=p + scale_fill_manual(values = wes_palette("Darjeeling1", n = 3))
      p=p + scale_fill_manual(values=colors)
      p
      ggsave(file.path(sub_dir,paste0(cutoff,'_odd_ratio_all_pre.pdf')),p,width = 10,height = 4)
      
      
      ## Box plot 
      # AUC
      p = ggplot(df_all, aes(x=Group, y=AUC)) +
        geom_boxplot(color="#1f77b4",alpha=0.3,outlier.shape = NA) +
        # Box plot with dot plot
        geom_jitter(aes(colour = Cohort,shape=Cohort), position=position_jitter(0.2),size=2)+
        theme_classic()+theme(axis.text.x = element_text(hjust = 0.5,vjust = 0.5,angle = 45))+
        scale_color_manual(values=colors)+
        scale_shape_manual(values=seq(0,15))+
        xlab('')+
        geom_hline(yintercept = 0.5,linetype='dashed',color="#d62728")
      p
      ggsave(file.path(sub_dir,paste0(cutoff,'_AUC_all_box.pdf')),p,width = 6.5,height = 4.5)
      saveRDS(p,file.path(sub_dir,paste0(cutoff,'_AUC_all_box.rds')))
      
      # add pvalue
      p = p + stat_compare_means(method = "wilcox.test",paired = T,ref.group = "CTL",label='p.format',method.args = list(alternative = "less")) # other groups compare to ref group, the alternative here should be "less"
      p
      
      
      ggsave(file.path(sub_dir,paste0(cutoff,'_AUC_all_Riaz_pre_box_pvalue.pdf')),p,width = 6.5,height = 4.5)
      
      # odd ratio
      p = ggplot(df_all, aes(x=Group, y=odd_ratio)) +
        geom_boxplot(color="#1f77b4",alpha=0.3,outlier.shape = NA) +
        # Box plot with dot plot
        geom_jitter(aes(colour = Cohort,shape=Cohort), position=position_jitter(0.2),size=2)+
        theme_classic()+theme(axis.text.x = element_text(hjust = 0.5,vjust = 0.5,angle = 45))+
        scale_color_manual(values=colors)+
        scale_shape_manual(values=seq(0,15))+
        xlab('')+ylab('Odd ratio')
      p
      ggsave(file.path(sub_dir,paste0(cutoff,'_odd_ratio_all_box.pdf')),p,width = 6.5,height = 4.5)
      saveRDS(p,file.path(sub_dir,paste0(cutoff,'_odd_ratio_all_Riaz_pre_box.rds')))
      
      # add pvalue
      p = p + stat_compare_means(method = "wilcox.test",paired = F,ref.group = "CTL",label='p.format',method.args = list(alternative = "less")) # other groups compare to ref group, the alternative here should be "less"
      ggsave(file.path(sub_dir,paste0(cutoff,'_odd_ratio_all_box_pvalue.pdf')),p,width = 6.5,height = 4.5)
    }
  }
}

## ------- AUC between decoupling high and low groups -------
## Calculate AUC and adds ratio for each signature at cohort level
odd_ratio_sum=data.frame()
for (s in c(ctl_id,'IMPRES','TIDE','CD274','CXCL9','CD8A','Proliferation','Stroma_EMT','T_effector_IFNG_response','CD8_T_effector_POPLAR','TGFB','Myeloid_inflammation','12_cheomokine')){
  print(s)
  # percentage 
  for (i in names(icb_skcm)){
    print(i)
    df_score=icb_skcm[[i]]$sig_score
    icb_score=icb_sig_score[[i]]
    sig_score=merge(df_score,icb_score,by='sample_id')
    
    meta_df=icb_skcm[[i]][['meta']]
    meta_df=meta_df[which(meta_df$Timepoint=='PRE'),]
    meta_df=meta_df[!is.na(meta_df$Response),]
    tmp_meta_pre_sample=meta_df$RNA_ID
    
    sig_score=sig_score[which(sig_score$sample_id%in%tmp_meta_pre_sample),]
    
    for (n in unique(meta_df$Treatment)){
      tmp_meta=meta_df[which(meta_df$Treatment==n),]
      
      # multiple cutoffs points  
      tmp_res=mult_odd_ratio(sig_score=sig_score,meta_df=tmp_meta,score_id=s,group=s,cohort=paste(gsub('_.+','',i),n,sep = ' '))
      odd_ratio_sum=rbind(odd_ratio_sum,tmp_res)
    }
  }
}

## barplot of AUC and Odd ratio
# cohort_sel=c("Gide PD1","Gide CTLA4","Riaz Pre PD1","Riaz Pre CTLA4","Liu PD1","Liu CTLA4" )
# percentage
bar_plot(df=odd_ratio_sum,decouple_group=c('High','All'),sub_dir='high_vs_all')
write.table(odd_ratio_sum,file.path(outdir,'high_vs_all','odd_ratio_sum.txt'),quote = F,sep = '\t')

