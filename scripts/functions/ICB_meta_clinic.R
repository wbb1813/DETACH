options(java.parameters = "-Xmx4096m")
library(xlsx)
## ------- Inputs -------
## ICB meta data
# Gide, PD1 and PD1 + CTLA4
meta_gide_pd1=read.xlsx('/Users/wangb8/Documents/Research/projects/cold2hot/data/COLDFACS/Gide_antiPD1.xlsx',sheetIndex = 1)
meta_gide_ctla4=read.xlsx('/Users/wangb8/Documents/Research/projects/cold2hot/data/COLDFACS/Gide_antiPD1_antiCTLA4.xlsx',sheetIndex = 1)
colnames(meta_gide_pd1)[which(colnames(meta_gide_pd1)=='Best.RECIST.response')]='BR'
colnames(meta_gide_ctla4)[which(colnames(meta_gide_ctla4)=='Best.RECIST.response')]='BR'

meta_gide_pd1$BR=factor(meta_gide_pd1$BR,c('CR','PR','SD','PD'))
meta_gide_pd1$Response='Non-response'
meta_gide_pd1$Response[which(meta_gide_pd1$BR%in%c('CR','PR'))]='Response'
meta_gide_pd1_pre=meta_gide_pd1[grep('PRE',meta_gide_pd1$Alias),]
meta_gide_pd1_on=meta_gide_pd1[grep('EDT',meta_gide_pd1$Alias),]

meta_gide_ctla4$BR=factor(meta_gide_ctla4$BR,c('CR','PR','SD','PD'))
meta_gide_ctla4$Response='Non-response'
meta_gide_ctla4$Response[which(meta_gide_ctla4$BR%in%c('CR','PR'))]='Response'
meta_gide_ctla4_pre=meta_gide_ctla4[grep('PRE',meta_gide_ctla4$Alias),]
meta_gide_ctla4_on=meta_gide_ctla4[grep('EDT',meta_gide_ctla4$Alias),]

meta_gide_pre=rbind(meta_gide_pd1_pre,meta_gide_ctla4_pre)
meta_gide_on=rbind(meta_gide_pd1_on,meta_gide_ctla4_on)

# Raiz, PD1, pre and on treatmen
meta_raiz_1=readRDS('/Users/wangb8/Documents/Research/projects/cold2hot/data/COLDFACS/Riaz_meta.rds')
meta_raiz_2=read.csv('/Users/wangb8/Documents/Research/projects/cold2hot/data/COLDFACS/timchan_clinical_data.csv')
meta_raiz_1$PatientID=gsub('\\_.+','',meta_raiz_1$title)
meta_raiz_2=meta_raiz_2[,c('PatientID','TRTGRP')]
meta_raiz=merge(meta_raiz_1,meta_raiz_2,by='PatientID',all.x = T)
colnames(meta_raiz)[which(colnames(meta_raiz)=='visit (pre or on treatment):ch1')]='pre_on_treat'
colnames(meta_raiz)[which(colnames(meta_raiz)=='response:ch1')]='BR'

meta_raiz=meta_raiz[which(meta_raiz$BR!='UNK'),] # remove unknown samples
meta_raiz$BR=factor(meta_raiz$BR,c('PRCR','SD','PD'))
meta_raiz$Response='Non-response'
meta_raiz$Response[which(meta_raiz$BR=='PRCR')]='Response'

meta_raiz_pre=meta_raiz[which(meta_raiz$pre_on_treat=='Pre'),]
meta_raiz_on=meta_raiz[which(meta_raiz$pre_on_treat=='On'),]

meta_raiz_pre_pd1=meta_raiz_pre[which(meta_raiz_pre$TRTGRP=='NIV3-NAIVE'),]
meta_raiz_pre_ctla4=meta_raiz_pre[which(meta_raiz_pre$TRTGRP=='NIV3-PROG'),]

meta_raiz_on_pd1=meta_raiz_on[which(meta_raiz_on$TRTGRP=='NIV3-NAIVE'),]
meta_raiz_on_ctla4=meta_raiz_on[which(meta_raiz_on$TRTGRP=='NIV3-PROG'),]

# survival data
meta_riaz_clinical=read.csv('/Users/wangb8/Documents/Research/projects/cold2hot/data/COLDFACS/Riaz_clinical_data.csv')

# Liu, priorCTLA4 and PD1, PD1
meta_liu=read.xlsx('/Users/wangb8/Documents/Research/projects/cold2hot/data/COLDFACS/Liu2019_meta.xlsx',sheetIndex = 1) 
colnames(meta_liu)=meta_liu[1,]
colnames(meta_liu)[1]='Run'
meta_liu=meta_liu[-1,]

colnames(meta_liu)[which(colnames(meta_liu)=="biopsyContext (1=Pre-Ipi; 2=On-Ipi; 3=Pre-PD1; 4=On-PD1)")]='pre_on_treat'
meta_liu$BR=factor(meta_liu$BR,levels = c('CR','PR','MR','SD','PD'))
meta_liu$Response='Non-response'
meta_liu$Response[which(meta_liu$BR%in%c('CR','PR'))]='Response'

meta_liu_ctla4=meta_liu[which(meta_liu$priorCTLA4==1),]
meta_liu_pd1=meta_liu[which(meta_liu$priorCTLA4!=1),]

# Puch
load('/Users/wangb8/Documents/Research/projects/cold2hot/data/COLDFACS/melanoma_meta_org_newcohort.rdata')
meta_puch=meta_org$puch
colnames(meta_puch)[3]='Run'
colnames(meta_puch)[11]='Response'
meta_puch$Response[which(meta_puch$Response==1)]='Response'
meta_puch$Response[which(meta_puch$Response==0)]='Non-response'


## ------- Process clinical data -------
# Gide
gide_clinical_pd1=meta_gide_pd1_pre[c('Run','Age..Years.','Sex','Progression.Free.Survival..Days.','Overall.Survival..Days.','Last.Followup.Status')]
colnames(gide_clinical_pd1)=c('Run','Age','Sex','PFS','OS','status')
gide_clinical_pd1$status[which(gide_clinical_pd1$status=='Alive')]=1 # censored
gide_clinical_pd1$status[which(gide_clinical_pd1$status=='Dead, melanoma')]=2 # dead
gide_clinical_pd1$status[which(gide_clinical_pd1$status=='Dead')]=2 # dead
gide_clinical_pd1$status=as.numeric(gide_clinical_pd1$status)

gide_clinical_ctla4=meta_gide_ctla4_pre[c('Run','Age..Years.','Sex','Progression.Free.Survival..Days.','Overall.Survival..Days.','Last.Followup.Status')]
colnames(gide_clinical_ctla4)=c('Run','Age','Sex','PFS','OS','status')
gide_clinical_ctla4$status[which(gide_clinical_ctla4$status=='Alive')]=1 # censored
gide_clinical_ctla4$status[which(gide_clinical_ctla4$status=='Dead, melanoma')]=2 # dead
gide_clinical_ctla4$status[which(gide_clinical_ctla4$status=='Dead')]=2 # dead
gide_clinical_ctla4$status=as.numeric(gide_clinical_ctla4$status)

# Riaz
riaz_clinical=meta_riaz_clinical[,c("Patient","Dead.Alive..Dead...True.","Time.to.Death..weeks.","M.Stage")]
colnames(riaz_clinical)=c('Patient','raw_status','OS','Stage')
riaz_clinical$status=1 # censored
riaz_clinical$status[which(riaz_clinical$raw_status=='FALSE')]=2 # dead
riaz_clinical$status=as.numeric(riaz_clinical$status)

riaz_clinical_pre=merge(meta_raiz_pre[,c('PatientID','Run')],riaz_clinical,by.x = 'PatientID',by.y='Patient')
riaz_clinical_on=merge(meta_raiz_on[,c('PatientID','Run')],riaz_clinical,by.x = 'PatientID',by.y='Patient')

# Liu
liu_clinical=meta_liu[,c("Run","gender (Male=1, Female=0)","PFS","OS","Mstage (IIIC=0, M1a=1, M1b=2, M1c=3)","dead")]
colnames(liu_clinical)=c('Run','Sex',"PFS","OS",'stage','raw_status')
liu_clinical$status=1 # censored
liu_clinical$status[which(liu_clinical$raw_status==1)]=2 # dead
liu_clinical$Sex[which(liu_clinical$Sex==1)]='M' 
liu_clinical$Sex[which(liu_clinical$Sex==0)]='F'
liu_clinical$status=as.numeric(liu_clinical$status)
liu_clinical$OS=as.numeric(liu_clinical$OS)
liu_clinical$PFS=as.numeric(liu_clinical$PFS)

# Puch
puch_clinical=meta_puch[,c("Run","Response","pfs_status","pfs_days","os_status","os_days")]

