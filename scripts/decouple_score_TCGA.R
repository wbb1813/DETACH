source('./functions/signature_score.R')

## ------- Input -------
# skcm_rsem_log2_file='../../data/skcm_rsem_log2.txt'
# inter_regression_res_file='../../results/bulk/decoupling/gsea/res_interaction.txt'
# outdir='../../results/bulk/decouple_score'

## TCGA RSEM data
skcm_rsem_log2=readRDS(skcm_rsem_log2_file)

## ------- signatures -------
## Decoupling ctl and fixed signature
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

## ------- Decoupling score for TCGA samples ------- 
if (!dir.exists(outdir)){
  dir.create(outdir,recursive = T)
}

gsea_sig_score=calcualte_feature_score(sig_ls=decouple_sig_ls, expr_df=skcm_rsem_log2, gsva_sig=names(decouple_sig_ls),gsea_method='ssgsea')
write.table(gsea_sig_score,file.path(outdir,'gsea_sig_score.txt'),quote = F,sep = '\t',row.names = F)
