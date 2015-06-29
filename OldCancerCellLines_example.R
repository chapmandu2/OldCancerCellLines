#examples for CellLineSensitivity


  sourceloc <- 'OldCancerCellLines_functions.R'
  gsc_file <- '/Users/pchapman/BigData/GSEA/c5.bp.v4.0.symbols.gmt'
  workdirtxt <- '/Users/pchapman/Documents/2014 Projects/20140327 CellLineSensitivity Package/test'
  dbpathtxt <- '/Users/pchapman/BigData/CellLineData/CellLineData.db'

#source functions
source(sourceloc)

#set user defined parameters
genelist1 <- paste ( getGenesOfBioProcess ('DNA_REPAIR', mode='symbol', gsc_file=gsc_file) , collapse=' ')
genelist2 <- 'PARG PARP1 PARP2 PARP3 BRCA1 BRCA2 FANCA FANCB FANCD FANCE'
NCG40_list <- read.table(NGC40_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
NCG40_list <- paste(unique(NCG40_list$symbol[1:300]), collapse=' ')

params <- set_params(workdir=workdirtxt,
                     dbpath=dbpathtxt,
                     cls.txt='cl1 cl2 cl3 cl4',
                     gene.symbols.txt=NCG40_list,
                     compounds.txt='Erlotinib Lapatinib AZD6244 ZD-6474',
                     respdata.source='ccle',
                     annmap_version='hs72')

#calculate derived parameters
params <- setup_params(params)
params <- update_params(params)

#get cell line info from db
sampleinfo <- subset(params$allsampleinfo, params$allsampleinfo[,'Site_Primary']=='breast')
params$cls.ids <- sampleinfo$CCLE_name
params <- update_params(params)

#get  data individually
rnaseq_data <- get_rnaseq (params)
affy_data <- get_affy (params)
hybcap_data <- get_hybcap (params)
cosmicclp_data <- get_cosmicclp (params)
cn_data <- get_cn (params)

#make heatmap with all data
make_heatmap(params)

#make heatmap with specified data types
make_heatmap(params, datatype=c('rnaseq' , 'hybcap', 'affy', 'cn'))

#make a data matrix
my_mat <- make_matrices(params)


#NOW DO MANOVA
manova_out <- run_MANOVA (params)
for (i in 1:24) {
make_manova_volcanoplot (params,manova_out, i, effect_th=0.2, log10FDR_th=1)
}

#plots for lapatanib
make_manova_volcanoplot (params,manova_out, 9, effect_th=0, log10FDR_th=1)
make_manova_heatmap (params,manova_out, 9, effect_th=0, log10FDR_th=1)

#plots for PLX4720
make_manova_volcanoplot (params,manova_out, 16, effect_th=0, log10FDR_th=1)
make_manova_heatmap (params,manova_out, 16, effect_th=0, log10FDR_th=1)

#plots for AZD6244
make_manova_volcanoplot (params,manova_out, 4, effect_th=0, log10FDR_th=1)
make_manova_heatmap (params,manova_out, 4, effect_th=0, log10FDR_th=1)

#plots for Vandetanib
make_manova_volcanoplot (params,manova_out, 24, effect_th=0, log10FDR_th=1)
make_manova_heatmap (params,manova_out, 24, effect_th=0, log10FDR_th=1)

#plots for Sorafenib
make_manova_volcanoplot (params,manova_out, 20, effect_th=0, log10FDR_th=1)
make_manova_heatmap (params,manova_out, 20, effect_th=0, log10FDR_th=1)

#################
#apply glmnet function instead of stream rr function
#make a data matrix
ct <- proc.time()
my_df <- make_df(params, datatype=c('affy', 'cn', 'hybcap'))
proc.time()-ct
#my_df <- inner_join(my_df, select(params$allsampleinfo, CCLE_name, Site_Primary), by='CCLE_name') %>% mutate(Site_Primary=as.factor(Site_Primary))

#run ridge regression/lasso/elastic net
drugs <- c(4:6,9,11,16)
ct <- proc.time()
glmnet_res <- run_glmnet (params, my_df, cpds=drugs)
proc.time()-ct
#now do various plots
#for AZD6244
make_rr_betaplot (params, glmnet_res, 4, n=20)
plot(glmnet_res[[4]], metric = "RMSE", plotType = "line")
plot(glmnet_res[[4]], metric = "RMSE", plotType = "level") 
make_rr_heatmap ( params, glmnet_res, 4, n=10)

#for Erlotinib
make_rr_betaplot (params, glmnet_res, 5, n=20)
plot(glmnet_res[[5]], metric = "RMSE", plotType = "line")
plot(glmnet_res[[5]], metric = "RMSE", plotType = "level") 
make_rr_heatmap ( params, glmnet_res, 5, n=10)

#for Irinotecan
make_rr_betaplot (params, glmnet_res, 6, n=20)
plot(glmnet_res[[6]], metric = "RMSE", plotType = "line")
plot(glmnet_res[[6]], metric = "RMSE", plotType = "level") 
make_rr_heatmap ( params, glmnet_res, 6, n=10)

#for lapatinib
make_rr_betaplot (params, glmnet_res, 9, n=20)
plot(glmnet_res[[9]], metric = "RMSE", plotType = "line")
plot(glmnet_res[[9]], metric = "RMSE", plotType = "level") 
make_rr_heatmap ( params, glmnet_res, 9, n=10)

#for Nutlin 3
make_rr_betaplot (params, glmnet_res, 11, n=20)
plot(glmnet_res[[11]], metric = "RMSE", plotType = "line")
plot(glmnet_res[[11]], metric = "RMSE", plotType = "level") 
make_rr_heatmap ( params, glmnet_res, 11, n=10)

#for PLX4720
make_rr_betaplot (params, glmnet_res, 16, n=20)
plot(glmnet_res[[16]], metric = "RMSE", plotType = "line")
plot(glmnet_res[[16]], metric = "RMSE", plotType = "level") 
make_rr_heatmap ( params, glmnet_res, 16, n=10)

#put plots in pdf
pdf('glmnet_plots.pdf', width = 8, height=12, onefile = TRUE)
for (i in drugs) {
  print(make_rr_betaplot (params, glmnet_res, i, n=20))
  print(plot(glmnet_res[[i]], metric = "RMSE", plotType = "line")) 
  print(plot(glmnet_res[[i]], metric = "Rsquared", plotType = "line")) 
  print(make_rr_heatmap ( params, glmnet_res, i, n=10))  
}
dev.off()


#########
#run boosting
##########
my_df <- make_df(params, datatype=c('affy', 'cn', 'hybcap'))
drugs <- c(4:6,9,11,16)
#gbm_res <- run_gbm (params, my_df, cpds=drugs)  #quick
ct <- proc.time()
gbm_res <- run_gbm (params, my_df, cpds=drugs, n.trees = (1:30)*100, shrinkage = c(0.001, 0.003, 0.01), interaction.depth = c(1,3,5,7)  )  #slow
proc.time()-ct

#now do various plots
#for AZD6244
make_rr_betaplot (params, gbm_res, 4, n=20)
plot(gbm_res[[4]], metric = "RMSE", plotType = "line")
plot(gbm_res[[4]], metric = "RMSE", plotType = "level") 
make_rr_heatmap ( params, gbm_res, 4, n=10)

#for Erlotinib
make_rr_betaplot (params, gbm_res, 5, n=20)
plot(gbm_res[[5]], metric = "RMSE", plotType = "line")
plot(gbm_res[[5]], metric = "RMSE", plotType = "level") 
make_rr_heatmap ( params, gbm_res, 5, n=10)

#put plots in pdf
pdf('gbm_plots.pdf', width = 8, height=12, onefile = TRUE)
for (i in drugs) {
  print(make_rr_betaplot (params, gbm_res, i, n=20))
  print(plot(gbm_res[[i]], metric = "RMSE", plotType = "line")) 
  print(plot(gbm_res[[i]], metric = "Rsquared", plotType = "line")) 
  print(make_rr_heatmap ( params, gbm_res, i, n=10))
}
dev.off()




