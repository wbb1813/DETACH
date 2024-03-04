## ------- Inputs -------
icb_skcm_file='../data/ICB_SKCM.RDS' 
ctl_id='CytotoxicTcell'
etl_id='fixed_exhuasted'
react_id='CD8_NeoTCR8' 
decouple_id='neg_genes_fdr0.01' 
ctrl_id='con_neg_genes_fdr0.01' 
icb_meta='../data/ICB_meta.RData'
icb_response_score='../data/icb_response_score.rds'
signature_file='../data/sig_ls.rds'

result_dir='../results/'

## ------- Step 1: Calculating ETL and CTL score for TCGA samples' -------
message('Step 1: Calculating ETL and CTL score for TCGA samples')
## 1.1 TCGA 
outdir=file.path(result_dir,'ctl_etl_score/TCGA')

source('./ctl_etl_score_TCGA.R')
## -------  Step 1: Decoupling, TCGA SKCM -------
message('Step 1: Decoupling, TCGA SKCM')
skcm_rsem_log2_file='../data/skcm_rsem_log2.rds'
etl_ctl_score_file='../data/TCGA_ETL_CTL_sig_score.txt'
outdir=file.path(result_dir,'decoupling')

source('./decoupling.R')

## ------- Step 2: Decoupling score -------
message('Step 2: Calculating decoupling score')

## 2.1 TCGA
skcm_rsem_log2_file='../data/skcm_rsem_log2.rds'
inter_regression_res_file=file.path(result_dir,'decoupling/res_interaction.txt') 
outdir=file.path(result_dir,'decouple_score/TCGA')

source('./decouple_score_TCGA.R')

## 2.2 ICB
inter_regression_res_file=file.path(result_dir,'decoupling/gsea/res_interaction.txt') # 3.1, 3.2
outdir=file.path(result_dir,'icb_sig_score/ICB')

source('./ICB_sig_score.R')

## ------- Step 3: decouple CTL and ETL correlation -------
message('Step 3: Adjust CTL and ETL correlation with decoupling score')

## TCGA
etl_ctl_score_path=file.path(result_dir,'ctl_etl_score/TCGA/')
decouple_score_path=file.path(result_dir,'decouple_score/TCGA/')
outdir=file.path(result_dir,'ctl_etl_cor_decouple/TCGA')

source('./ctl_etl_cor_decouple_TCGA.R')

## ICB
icb_sig_score_path=file.path(result_dir,'icb_sig_score/ICB/')
outdir=file.path(result_dir,'ctl_etl_cor_decouple/ICB')

source('./ctl_etl_cor_decouple_ICB.R')

## ------- Step 4:  Adjust CTL and ETL correlation with decoupling score -------
message('Step 4: Adjust CTL and ETL correlation with decoupling score')
## TCGA
etl_ctl_score_path=file.path(result_dir,'ctl_etl_score/TCGA/')
decouple_score_path=file.path(result_dir,'decouple_score/TCGA/')
outdir=file.path(result_dir,'ctl_etl_cor_decouple/TCGA')

source('./ctl_etl_cor_decouple_TCGA.R')

## ICB
icb_sig_score_path=file.path(result_dir,'icb_sig_score/ICB/')
outdir=file.path(result_dir,'ctl_etl_cor_decouple/ICB')

source('./ctl_etl_cor_decouple_ICB.R')

## merge TCGA and ICB
tcga_cor_hl_file=file.path(result_dir,'ctl_etl_cor_decouple/TCGA/high_vs_low/cor_decouple_hl.txt')
tcga_cor_partial_file=file.path(result_dir,'ctl_etl_cor_decouple/TCGA/partial/partial_cor.txt' )
icb_cor_hl_file=file.path(result_dir,'ctl_etl_cor_decouple/ICB/high_vs_low/cor_decouple_uncouple.txt' )
icb_cor_partial_file=file.path(result_dir,'ctl_etl_cor_decouple/ICB/partial/cor_partial.txt')
tcga_skcm_frac_file='../data/tcga_skcm_frac.rds'
outdir=file.path(result_dir,'ctl_etl_cor_decouple/TCGA_ICB')

source('./ctl_etl_cor_decouple_merge.R')

## ------- Step 5: ICB response -------
message('Step 5: ICB response')
icb_sig_score_path=file.path(result_dir,'icb_sig_score/ICB/')
outdir=file.path(result_dir,'ICB_response')

source('./ICB_response.R')

## ------- Step 6: TIL infiltration -------
message('Step 6: TIL infiltration')
icb_sig_score_file=file.path(result_dir,'icb_sig_score/ICB/icb_skcm_sig_score.rds')
tcga_etl_ctl_score_file=file.path(result_dir,'ctl_etl_score/TCGA/gsea_sig_score.txt')
decouple_score_file=file.path(result_dir,'decouple_score/TCGA/gsea_sig_score.txt')
tcga_other_sig_score_file='../data/tcga_skcm_sig_score.txt'
tcga_clinic_file='../data/tcga_clinical.tsv'
cell_frac_file='../data/cell_fraction.txt'
til_file='../data/TCGA_pathology_TIL.txt'
outdir=file.path(result_dir,'TIL_infiltration')

source('./TIL_infiltration.R')

