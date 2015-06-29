#get shRNA knockdown data from achilles 

#source functions
source('OldCancerCellLines_functions.R')

#set genes
genelist1 <- paste ( getGenesOfBioProcess ('DNA_REPAIR', mode='symbol', gsc_file='/Users/pchapman/BigData/GSEA/c5.bp.v4.0.symbols.gmt') , collapse=' ')
genelist6 <- 'HNF4A UBE2T RET MDM2 MDM4 MET ERRB2 EGFR EGFR ERBB2 HRAS KRAS ATM NRAS PIK3CA PTEN BRAF TP53 KDR PDGFRB CRAF MAP3K1 SMARCA4 SMARCA2 ARID1A ARID1B ATM'

#set user defined parameters
params <- set_params(workdir='/Users/pchapman/Documents/2014 Projects/20140327 CellLineSensitivity Package/achilles_example/',
                     dbpath='/Users/pchapman/BigData/CellLineData/CellLineData.db',
                     cls.txt='cl1 cl2 cl3 cl4',
                     gene.symbols.txt <- paste(genelist1, genelist6),
                     annmap_version='hs74',
                     compounds.txt='SMARCA2 BRAF PARG DNMT1 HNF4A',
                     respdata.source='achilles')

#calculate derived parameters
params <- setup_params(params)
params <- update_params(params)

#get cell line info from db
params$cls.ids <- unique(params$allsampleinfo[,'CCLE_name'])
params <- update_params(params)

#test the get functinos
#x <- get_response(params)
#y <- get_affy(params)
#y <- get_hybcap(params)
#y <- get_cn(params)
#y <- get_cosmicclp(params)
ct <- proc.time()
my_mat <- make_matrices(params, c('resp','hybcap', 'affy', 'cn', 'rnaseq'))
print(proc.time()-ct)

ct <- proc.time()
my_df <- make_df(params, c('hybcap', 'affy', 'cn'))
print(proc.time()-ct)

#get the data
my_df<- make_df(params, datatype=c('affy' , 'cn', 'hybcap', 'cosmic'))

#use wrapper function to do miltiple glmnets
this_gene <- 'SMARCA2'
cpds = grep(this_gene, params$compounds)
glmnet_res <- run_glmnet (params, my_df, cpds=cpds) 

#plot the basic data from achilles
achilles_data <-dbGetQuery(params$con,  sprintf("select * from achilles_v243 where shRNA_ID IN (%s) ", params$compounds.sql))
plot_data <- achilles_data %>% filter(Symbol==this_gene) 
p <- ggplot(plot_data, aes(x=log2FC, color=Symbol)) + geom_density() + facet_wrap(~ shRNA_ID)
print(p)
#make a pdf
pdf(file='test2.pdf', width = 10, height = 10, onefile = TRUE)
print(p)
for (i in cpds) {
  print(i)
  print(plot(glmnet_res[[i]], metric = "RMSE", plotType = "line"))
  print(plot(glmnet_res[[i]], metric = "Rsquared", plotType = "line"))
  print(make_rr_betaplot (params, glmnet_res, i, n=20))
  print(make_rr_heatmap (params, glmnet_res, i, n=10))
  print(do_simple_analysis(my_mat, i,  'SMARCA4_hybcap', feature_type='discrete', resp_type='continuous')$combo_plot)
}
dev.off()

#try manova
#NOW DO MANOVA
manova_out <- run_MANOVA (params)
for (i in cpds) {
  make_manova_volcanoplot (params,manova_out, i, effect_th=0.2, log10FDR_th=0.5)
}

#finally let's have a look at how the models compare
resamps <- resamples(list(GBM = gbmFit, GLMNET = glmnetFit))
resamps <- resamples(glmnet_res)
summary(resamps)
bwplot(resamps)

