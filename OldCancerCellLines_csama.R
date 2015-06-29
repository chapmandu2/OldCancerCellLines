# Example for CSAMA 2015

#some housekeeping
options(warn=-1)

sourceloc <- 'OldCancerCellLines_functions.R'
gsc_file <- '/Users/pchapman/BigData/GSEA/c5.bp.v4.0.symbols.gmt' #downloaded from Broad MSigDB
workdirtxt <- '/Users/pchapman/Documents/2014 Projects/20140327 CellLineSensitivity Package/test'
dbpathtxt <- '/Users/pchapman/BigData/CellLineData/CellLineData.db'

#load the functions
source(sourceloc)

#define a gene list
genelist <- 'ERBB2 EGFR BRAF NRAS PTEN TP53'

#set up the parameters object
params <- set_params(workdir=workdirtxt,
                     dbpath=dbpathtxt,
                     cls.txt='cl1 cl2 cl3 cl4',
                     gene.symbols.txt=genelist,
                     compounds.txt='Erlotinib Lapatinib AZD6244 ZD-6474',
                     respdata.source='ccle',
                     annmap_version='hs72')

#process the initial params object
params <- setup_params(params)
params <- update_params(params)

#define the cell lines
sampleinfo <- subset(params$allsampleinfo, params$allsampleinfo[,'Site_Primary']=='breast')
params$cls.ids <- sampleinfo$CCLE_name
params <- update_params(params)

#make a data frame
my_df <- make_df(params)
my_df[1:10,1:8]
colnames(my_df)

#make a data frame with only some data types
my_df <- make_df(params, datatype = c('resp', 'affy', 'hybcap'))
my_df[1:10,1:8]
colnames(my_df)

#make a list of two data matrices
my_mat <- make_matrices(params, datatype = c('resp', 'affy', 'hybcap'))
my_mat$featureData[1:10,1:6]
my_mat$responseData[1:10,]

#now make some heatmaps of the data
make_heatmap(params, datatype=c( 'hybcap', 'affy', 'cn'))
make_heatmap(params, datatype=c( 'resp', 'hybcap', 'affy', 'cn'), compound = 'Lapatinib')
make_heatmap(params, datatype=c( 'resp', 'hybcap', 'affy', 'cn'),  order_by = 2, multi_cpds = TRUE)

###################

#now lets do some modeling with a bigger gene set
params$gene.symbols.txt <- 'PARG PARP1 PARP2 PARP3 PARP4 PARP5 MDM2 MDM4 RET PARP6 ATM BRAF CDKN2A RB1 SMAD4 MET ERRB2 EGFR BRCA1 BRCA2 FANCA FANCB FANCC FANCE FANCG FANCF FANCM XRCC1 XRCC2 BLM SMARCA4 SMARC2 EGFR ERBB2 CDKN2A CTNNA1 EGFR AKT1 GATA3 APC HRAS KRAS SMAD4 ATM NRAS PIK3CA PTEN RB1 MAP2K4 BRCA1 BRCA2 BRAF TP53 MUC16 CDH1 KDR PDGFRB CRAF'
params$cls.ids <- params$allsampleinfo$CCLE_name
params <- update_params(params)

#make the data objects
my_df <- make_df(params, datatype=c('affy', 'cn', 'hybcap'))
dim(my_df)
my_mat <- make_matrices(params, datatype=c('affy', 'cn', 'hybcap'))


#feed into a glmnet wrapper function (uses caret)
glmnet_res <- run_glmnet (params, my_df, cpds=1:4)

#view the results for compound 3 - AZD6244, a MEK inhibitor (BRAF and NRAS driven tumours)
plot(glmnet_res[[3]], metric = "RMSE", plotType = "line")
make_rr_betaplot (params, glmnet_res, 3, n=20)
make_rr_heatmap ( params, glmnet_res, 3, n=10)
do_simple_analysis(my_mat, 3,  'BRAF_hybcap', feature_type='discrete', resp_type='continuous')$combo_plot

#view the results for compound 2 - Lapatanib, an ERBB2 inhibitor
plot(glmnet_res[[2]], metric = "RMSE", plotType = "line")
make_rr_betaplot (params, glmnet_res, 2, n=20)
make_rr_heatmap ( params, glmnet_res, 2, n=10)
do_simple_analysis(my_mat, 2,  'ERBB2_affy', feature_type='continuous', resp_type='continuous')$combo_plot

#can also run an analysis of variance/t-test type analysis for mutation status vs response
manova_out <- run_MANOVA (params)
manova_out %>% filter(cpd_id == 3) %>% arrange(pval) %>% dplyr::slice(1:10)

#finally make a volcano plot for AZD6244
make_manova_volcanoplot (params,manova_out, 3, effect_th=0.2, log10FDR_th=1)

