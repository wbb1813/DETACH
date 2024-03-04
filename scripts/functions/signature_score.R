## ------- Functions -------
### single gene signature 
# CD8A CD38 CD274 HAVCR2 MEX3B CXCL9
single_gene_sig=function(sig_name,sig_ls,expr_df){
  sig_genes=NULL
  for (i in sig_name){
    sig_genes=c(sig_genes,sig_ls[[i]]$V1)
  }
  
  expr_df_sel=as.data.frame(t(expr_df[intersect(rownames(expr_df),sig_genes),]))
  expr_df_sel$sample_id=rownames(expr_df_sel)
  return(expr_df_sel)
}

### mean expression signature 
mean_expr_sig=function(sig_name,sig_ls,expr_df){
  sig_score=NULL
  for (i in sig_name){
    sig_genes=sig_ls[[i]]$V1
    expr_df_sel=expr_df[intersect(rownames(expr_df),sig_genes),,drop=FALSE]
    if (ncol(expr_df_sel)>=1){
      tmp_sig_score=as.data.frame(apply(expr_df_sel,2,function(x)mean(x)))
      colnames(tmp_sig_score)=i
      tmp_sig_score$sample_id=rownames(tmp_sig_score)
      if (is.null(sig_score)){
        sig_score=tmp_sig_score
      }else{
        sig_score=merge(sig_score,tmp_sig_score,by='sample_id')
      }
    }
  }
  return(sig_score)
}

### geometric mean expression signature 
geometric_mean_expr_sig=function(sig_name,sig_ls,expr_df){
  sig_score=NULL
  for (i in sig_name){
    sig_genes=sig_ls[[i]]$V1
    expr_df_sel=expr_df[intersect(rownames(expr_df),sig_genes),]
    tmp_sig_score=as.data.frame(apply(expr_df_sel,2,function(x)geometric.mean(x)))
    colnames(tmp_sig_score)=i
    tmp_sig_score$sample_id=rownames(tmp_sig_score)
    if (is.null(sig_score)){
      sig_score=tmp_sig_score
    }else{
      sig_score=merge(sig_score,tmp_sig_score,by='sample_id')
    }
  }
  return(sig_score)
}

### weight signature
weight_expr_sig=function(sig_name,sig_ls,expr_df){
  sig_score=NULL
  for (i in sig_name){
    sig_genes=sig_ls[[i]]
    rownames(sig_genes)=sig_genes[,1] # gene name at column 1
    colnames(sig_genes)[2]='weight' # weight at column 2
    
    expr_df_sel=expr_df[intersect(rownames(expr_df),rownames(sig_genes)),]
    tmp_sig_score=weight_score(gene_weight=sig_genes,expr_df)
    colnames(tmp_sig_score)[1]=i
    if (is.null(sig_score)){
      sig_score=tmp_sig_score
    }else{
      sig_score=merge(sig_score,tmp_sig_score,by='sample_id')
    }
  }
  return(sig_score)
}

### cheomokine signature, use the first component of PCA as the signature score
pca_expr_sig=function(sig_name,sig_ls,expr_df){
  sig_score=NULL
  for (i in sig_name){
    sig_genes=sig_ls[[i]]$V1
    expr_df_sel=expr_df[intersect(rownames(expr_df),sig_genes),]
    
    # remove samples with zero var
    sample_var=apply(expr_df_sel,2,function(x)var(x))
    names(sample_var)[which(sample_var==0)]
    expr_df_sel=expr_df_sel[,which(colnames(expr_df_sel)%in%setdiff(colnames(expr_df_sel),names(sample_var)[which(sample_var==0)]))]
    
    res.pca <- prcomp(expr_df_sel, center = TRUE,scale. = TRUE)
    pca_df=as.data.frame(res.pca$rotation)
    pca_df$sample_id=rownames(pca_df)
    pca_df=pca_df[,c('sample_id','PC1')]
    colnames(pca_df)[2]=i
    if (is.null(sig_score)){
      sig_score=pca_df
    }else{
      sig_score=merge(sig_score,pca_df,by='sample_id')
    }
  }
  return(sig_score)
}

### GSVA
library(GSVA)
gsva_expr_sig=function(sig_name,sig_ls,expr_df,gsea_method='gsva'){
  geneset_ls=list()
  for (i in sig_name){
    tmp_df=sig_ls[[i]]
    tmp_ls=list(tmp=tmp_df[,1])
    names(tmp_ls)=i
    geneset_ls=c(geneset_ls,tmp_ls)
  }
  if (gsea_method=='ssgsea'){
    gsva.es <- gsva(expr = as.matrix(expr_df), gset.idx.list = geneset_ls, verbose=FALSE,method='ssgsea')
  }else{
    gsva.es <- gsva(expr = as.matrix(expr_df), gset.idx.list = geneset_ls, verbose=FALSE,method='gsva')
  }
  gsva.es=as.data.frame(t(gsva.es))
  gsva.es$sample_id=rownames(gsva.es)
  return(gsva.es)
}

### IMPRES
IMPRES_sig=function(sig_name,sig_ls,expr_df){
  sig_score=NULL
  for (i in sig_name){
    sig_genes=sig_ls[[i]]
    sig_genes$gene1=0
    sig_genes$gene1[which(sig_genes$V1%in%rownames(expr_df))]=1
    sig_genes$gene2=0
    sig_genes$gene2[which(sig_genes$V2%in%rownames(expr_df))]=1
    sig_genes$pair=sig_genes$gene1+sig_genes$gene2
    sig_genes=sig_genes[which(sig_genes$pair==2),]
    sig_genes$id=paste(sig_genes$V1,sig_genes$V2,sep = '_')
    
    sig_genes=sig_genes[,c('V1','V2','id')]
    
    tmp_expr1=merge(sig_genes,expr_df,by.x = 'V1',by.y = 0)
    tmp_expr2=merge(sig_genes,expr_df,by.x = 'V2',by.y = 0)
    rownames(tmp_expr1)=tmp_expr1$id; tmp_expr1=tmp_expr1[,-c(1,2,3)]
    rownames(tmp_expr2)=tmp_expr2$id; tmp_expr2=tmp_expr2[,-c(1,2,3)]
    
    com_id=intersect(rownames(tmp_expr1),rownames(tmp_expr2))
    tmp_expr1=tmp_expr1[com_id,]
    tmp_expr2=tmp_expr2[com_id,]
    
    tmp_diff=as.matrix(tmp_expr1-tmp_expr2)
    tmp_diff[which(tmp_diff<0)]=1
    tmp_diff[which(tmp_diff>0)]=0
    
    tmp_score=as.data.frame(apply(tmp_diff, 2, function(x)sum(x)))
    colnames(tmp_score)=i
    tmp_score$sample_id=rownames(tmp_score)
    
    if (is.null(sig_score)){
      sig_score=tmp_score
    }else{
      sig_score=merge(sig_score,tmp_score,by='sample_id')
    }
  }
  return(sig_score)
}

### Chiping 45
chiping_sig=function(sig_name,sig_ls,expr_df){
  sig_score=NULL
  for (i in sig_name){
    sig_genes=sig_ls[[i]]
    rownames(sig_genes)=sig_genes$V2
    
    com_genes=intersect(sig_genes$V2,rownames(expr_df))
    
    sig_genes=sig_genes[com_genes,]
    expr_df=expr_df[com_genes,]
    
    tmp_score=expr_df*sig_genes$V1
    tmp_score=as.data.frame(apply(tmp_score,2,function(x)sum(x)))
    colnames(tmp_score)=i
    tmp_score$sample_id=rownames(tmp_score)
    
    if (is.null(sig_score)){
      sig_score=tmp_score
    }else{
      sig_score=merge(sig_score,tmp_score,by='sample_id')
    }
  }
  return(sig_score)
}

### Regression (Peng's NM 2022)
# with background information in the signature 
regression_expr_sig=function(sig_name,sig_ls,expr_df){
  sig_score=NULL
  for (i in sig_name){
    sig_genes=sig_ls[[i]]
    rownames(sig_genes)=sig_genes[,1] # gene name at column 1
    colnames(sig_genes)[2]='value' # binary value at column 2
    
    expr_df_sel=expr_df[intersect(rownames(expr_df),rownames(sig_genes)),]
    sig_genes_sel=sig_genes[intersect(rownames(expr_df),rownames(sig_genes)),]
    
    tmp_score=data.frame()
    for (j in colnames(expr_df_sel)){
      coef=summary(lm(expr_df_sel[,j]~sig_genes_sel[,2]))
      tmp_res=as.data.frame(coef$coefficients)
      #tmp_score=data.frame(sample_id=j,coef=tmp_res$Estimate[2],pvalue=tmp_res$`Pr(>|t|)`[2])
      tmp_res=data.frame(sample_id=j,coef=tmp_res$Estimate[2])
      colnames(tmp_res)[2]=i
      tmp_score=rbind(tmp_score,tmp_res)
    }
    
    if (is.null(sig_score)){
      sig_score=tmp_score
    }else{
      sig_score=merge(sig_score,tmp_score,by='sample_id')
    }
  }
  return(sig_score)
}
# use the whole transcript genes except the signature genes as background 
regression_expr_sig_no_bg=function(sig_name,sig_ls,expr_df){
  sig_score=NULL
  for (i in sig_name){
    
    sig_genes=sig_ls[[i]]
    sig_genes$value=1
    
    other_genes=setdiff(rownames(expr_df),sig_genes$V1)
    tmp_df=data.frame(V1=other_genes,value=0)
    
    sig_genes=rbind(sig_genes,tmp_df)
    sig_genes=sig_genes[!duplicated(sig_genes$V1),]
    
    rownames(sig_genes)=sig_genes[,1] # gene name at column 1
    colnames(sig_genes)[2]='value' # binary value at column 2
    
    expr_df_sel=expr_df[intersect(rownames(expr_df),rownames(sig_genes)),]
    sig_genes_sel=sig_genes[intersect(rownames(expr_df),rownames(sig_genes)),]
    
    tmp_score=data.frame()
    for (j in colnames(expr_df_sel)){
      coef=summary(lm(expr_df_sel[,j]~sig_genes_sel[,2]))
      tmp_res=as.data.frame(coef$coefficients)
      #tmp_score=data.frame(sample_id=j,coef=tmp_res$Estimate[2],pvalue=tmp_res$`Pr(>|t|)`[2])
      tmp_res=data.frame(sample_id=j,coef=tmp_res$Estimate[2])
      colnames(tmp_res)[2]=i
      tmp_score=rbind(tmp_score,tmp_res)
    }
    
    if (is.null(sig_score)){
      sig_score=tmp_score
    }else{
      sig_score=merge(sig_score,tmp_score,by='sample_id')
    }
  }
  return(sig_score)
}
## ------- calculate signature score -------
calcualte_feature_score=function(sig_ls, expr_df,single_sig=NULL, mean_sig=NULL, geometric_mean_sig=NULL, weight_sig=NULL, pca_sig=NULL,impress_sig=NULL,chiping45_sig=NULL,gsva_sig=NULL,reg_sig=NULL,gsea_method='gsva'){
  score_all=NULL
  ### single gene signature 
  # CD8A CD38 CD274 HAVCR2 MEX3B CXCL9
  # sig_name=c('CD8A','CD38','CD274','HAVCR2','MEX3B','CXCL9')
  if (! is.null(single_sig)){
    single_gene_sig_df=single_gene_sig(sig_name=single_sig,sig_ls,expr_df)
    if (is.null(score_all)){
      score_all=single_gene_sig_df
    }else{
      score_all=merge(score_all,single_gene_sig_df,by='sample_id')
    }
  }
  
  ### mean expression signature 
  # sig_name=c('Angiogenesis','Proliferation','Stroma_EMT','T_effector_IFNG_response','TGFB','CD8_T_effector_POPLAR','Myeloid_inflammation','DediffGenes')
  if (! is.null(mean_sig)){
    mean_sig_df=mean_expr_sig(mean_sig,sig_ls,expr_df)
    if (is.null(score_all)){
      score_all=mean_sig_df
    }else{
      score_all=merge(score_all,mean_sig_df,by='sample_id')
    }
  }
  
  ### geometric mean expression signature 
  # sig_name=c('Cytolytic_activity','ExhaustedTcell')
  if (! is.null(geometric_mean_sig)){
    geometric_mean_sig_df=geometric_mean_expr_sig(geometric_mean_sig,sig_ls,expr_df)
    if (is.null(score_all)){
      score_all=geometric_mean_sig_df
    }else{
      score_all=merge(score_all,geometric_mean_sig_df,by='sample_id')
    }
  }
  
  ### weight of each gene
  if (! is.null(weight_sig)){
    weight_sig_df=weight_expr_sig(sig_name=weight_sig,sig_ls,expr_df)
    if (is.null(score_all)){
      score_all=weight_sig_df
    }else{
      score_all=merge(score_all,weight_sig_df,by='sample_id')
    }
  }
  
  ### cheomokine signature, use the first component of PCA as the signature score
  # sig_name=c('12_cheomokine')
  if (! is.null(pca_sig)){
    pca_sig_df=pca_expr_sig(pca_sig,sig_ls,expr_df)
    if (is.null(score_all)){
      score_all=pca_sig_df
    }else{
      score_all=merge(score_all,pca_sig_df,by='sample_id')
    }
  }
  
  ### IMPRES
  # sig_name=c('IMPRES')
  if (! is.null(impress_sig)){
    IMPRES_sig_score=IMPRES_sig(impress_sig,sig_ls,expr_df)
    if (is.null(score_all)){
      score_all=IMPRES_sig_score
    }else{
      score_all=merge(score_all,IMPRES_sig_score,by='sample_id')
    }
  }
  
  ### Chiping 45
  # sig_name=c('Chiping_45_signature')
  if (! is.null(chiping45_sig)){
    chiping_sig_score=chiping_sig(chiping45_sig,sig_ls,expr_df)
    if (is.null(score_all)){
      score_all=chiping_sig_score
    }else{
      score_all=merge(score_all,chiping_sig_score,by='sample_id')
    }
  }
  
  ### GSVA signature
  if (!is.null(gsva_sig)){
    gsva_sig_score=gsva_expr_sig(gsva_sig,sig_ls,expr_df,gsea_method)
    if (is.null(score_all)){
      score_all=gsva_sig_score
    }else{
      score_all=merge(score_all,gsva_sig_score,by='sample_id')
    }
  }
  
  ### Regression signature
  # if (!is.null(reg_sig)){
  #   regression_sig_score=regression_expr_sig(reg_sig,sig_ls,expr_df)
  #   if (is.null(score_all)){
  #     score_all=regression_sig_score
  #   }else{
  #     score_all=merge(score_all,regression_sig_score,by='sample_id')
  #   }
  # }
  
  if (!is.null(reg_sig)){
    regression_sig_score=regression_expr_sig_no_bg(reg_sig,sig_ls,expr_df)
    if (is.null(score_all)){
      score_all=regression_sig_score
    }else{
      score_all=merge(score_all,regression_sig_score,by='sample_id')
    }
  }
  return(score_all)
}
