source('./functions/signature_score.R')

## ------- Input -------
## TCGA RSEM data
skcm_rsem_log2=readRDS(skcm_rsem_log2_file)

## signatures list
sig_ls=readRDS(signature_file)

## ------- cytolytic and exhausted score for each TCGA samples ------- 
if (!dir.exists(outdir)){
  dir.create(outdir,recursive = T)
}

gsea_sig_score=calcualte_feature_score(sig_ls=sig_ls, expr_df=skcm_rsem_log2, gsva_sig = names(sig_ls),gsea_method = 'ssgsea')

write.table(gsea_sig_score,file.path(outdir,'gsea_sig_score.txt'),quote = F,sep = '\t',row.names = F)

