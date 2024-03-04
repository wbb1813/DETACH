source('./functions/signature_score.R')

## ------- Input -------
## signatures list
sig_ls=readRDS(signature_file)

## Decoupling signature
decouple_sig_ls=list()

regression_res=read.delim(inter_regression_res_file)

pos_genes_fdr0.01=regression_res[which(regression_res$ext_gene_t_value>0&regression_res$ext_gene_FDR<=0.01),]
neg_genes_fdr0.01=regression_res[which(regression_res$ext_gene_t_value<0&regression_res$ext_gene_FDR<=0.01),]

decouple_sig_ls[['neg_genes_fdr0.01']]=data.frame(V1=neg_genes_fdr0.01$gene_symbol)
decouple_sig_ls[['pos_genes_fdr0.01']]=data.frame(V1=pos_genes_fdr0.01$gene_symbol)

pos_genes_fdr0.1=regression_res[which(regression_res$ext_gene_t_value>0&regression_res$ext_gene_FDR<=0.1),]
neg_genes_fdr0.1=regression_res[which(regression_res$ext_gene_t_value<0&regression_res$ext_gene_FDR<=0.1),]

decouple_sig_ls[['neg_genes_fdr0.1']]=data.frame(V1=neg_genes_fdr0.1$gene_symbol)
decouple_sig_ls[['pos_genes_fdr0.1']]=data.frame(V1=pos_genes_fdr0.1$gene_symbol)

set.seed(7)
decouple_sig_ls[['con_neg_genes_fdr0.01']]=data.frame(V1=sample(regression_res$gene_symbol,nrow(neg_genes_fdr0.01)))
decouple_sig_ls[['con_neg_genes_fdr0.1']]=data.frame(V1=sample(regression_res$gene_symbol,nrow(neg_genes_fdr0.1)))

sig_ls=c(sig_ls,decouple_sig_ls)

## Load COLDFACS ICB data
icb_skcm=readRDS(icb_skcm_file)

## calcualte signature score 
icb_score=list()
for (i in names(icb_skcm)){
  print(i)
  tmp_ls=icb_skcm[[i]]
  tmp_tpm=tmp_ls$tpm
  tmp_tpm=tmp_tpm[!duplicated(tmp_tpm$Gene),]
  rownames(tmp_tpm)=tmp_tpm$Gene
  tmp_tpm=tmp_tpm[,-1]
  tmp_tpm=log2(tmp_tpm+1)
  
  sig_score=calcualte_feature_score(sig_ls=sig_ls, expr_df=tmp_tpm, gsva_sig=names(sig_ls),gsea_method = 'ssgsea')
  icb_score[[i]]=icb_skcm[[i]]
  icb_score[[i]]$sig_score=sig_score
}

if (!dir.exists(outdir)){
  dir.create(outdir,recursive = T)
}

saveRDS(icb_score,file.path(outdir,'icb_skcm_sig_score.rds'))
