#functions for CellLineSensitivity package

# library import
library('ggplot2')
library('RSQLite')
library('reshape2')
library(plyr)
library(scales)
library(reshape)
library(XLConnect)
library(annmap)
library(piano)
library(gridExtra)
library(dplyr)
library(tidyr)

#set parameters for analysis
set_params <- function(workdir,dbpath,cls.txt,gene.symbols.txt,annmap_version='hs74',cls.ids=NULL,respdata=NULL, respdata.source='custom',compounds.txt=NULL) {
  return ( list (workdir=workdir,
                 dbpath=dbpath,
                 cls.txt=cls.txt,
                 gene.symbols.txt=gene.symbols.txt,
                 annmap_version=annmap_version,
                 cls.ids=cls.ids,
                 respdata=respdata,
                 respdata.source = respdata.source,
                 compounds.txt=compounds.txt) )
}

#initial setup of standard parameters from those given
setup_params <- function(params) {
  

    #setwd
    setwd(params$workdir)  
    #connect to annmap
    annmapConnect(params$annmap_version)
    #sqlite setup
    drv <- dbDriver("SQLite")
    con <- dbConnect(drv, dbname = params$dbpath)


    new <- list(drv=drv,
                con=con)
    
    return ( c(params,new)  )
  
}

#this part can be called at any time to update parameters
update_params <- function(params) {
  
  #sort out genes list
  gene.symbols <- unlist(strsplit(params$gene.symbols.txt, ' '))
  gene.symbols <- unique(gene.symbols)
  gene.symbols.sql <- paste("'", gene.symbols, "'", sep='', collapse=',')
  gene.ids <- symbolToGene(gene.symbols, as.vector=T)
  
  #sort out cell lines list
  if (length(params$cls.ids) == 0) {
    params$cls.ids <- unlist(strsplit(params$cls.txt, ' '))
  }
  
  #make cell lines sql
  cls.ids.sql <- paste("'", params$cls.ids, "'", sep='', collapse=',')
    
  #get cell line info
  allsampleinfo <- dbGetQuery(params$con, 
                        'select * 
                        from ccle_cell_line_info as t1 
                         inner join ccle_id_mapping as t2  ON t1.CCLE_name = t2.CCLE_name ')
  allsampleinfo <- allsampleinfo [ , !grepl('row_names', colnames(allsampleinfo))  ]
  
  #get cosmic cell line ids
  cls.ids.cosmic <- allsampleinfo [ allsampleinfo$CCLE_name %in% params$cls.ids , 'cosmicclp_id']
  cls.ids.cosmic <- cls.ids.cosmic [!is.na(cls.ids.cosmic)]
  
  #get compound ids
  if (params$respdata.source == 'custom') {
    compounds <- unique(params$respdata$CompoundID)
    compounds.sql <- ''
  } else {
    compounds <- unlist(strsplit(params$compounds.txt, ' '))
    compounds <- unique(compounds)
    compounds.sql <- paste("'", compounds, "'", sep='', collapse=',')
  }
  
  #need to deal with 'compounds' differently if using achilles data
  if (params$respdata.source == 'achilles') {
    shRNA_ID.sql <- sprintf('select distinct Symbol, shRNA_ID from achilles_v243 where Symbol IN (%s)', compounds.sql)
    shRNA_ID.df <- dbGetQuery(params$con, shRNA_ID.sql)
    compounds <- unique(shRNA_ID.df$shRNA_ID)
    compounds.sql <- paste("'", compounds, "'", sep='', collapse=',')
  }

  
  #put new and updated variables into params
  params$gene.ids <- gene.ids
  params$gene.symbols <- gene.symbols
  params$gene.symbols.sql <- gene.symbols.sql
  params$cls.ids.sql <- cls.ids.sql
  params$allsampleinfo <- allsampleinfo
  params$cls.ids.cosmic <- cls.ids.cosmic
  params$compounds <- compounds
  params$compounds.sql <- compounds.sql
  
  return ( params  )
  
}

get_rnaseq <- function (params) {
  
  #get affy data for genes and cell lines of interest
  sql <- sprintf('select * from ccle_rnaseq where symbol IN (%s) AND CCLE_name IN (%s)', params$gene.symbols.sql, params$cls.ids.sql) 
  rnaseq.data <- dbGetQuery(params$con, sql)
  
  #select and rename columns, filter and log transform RPKM
  rnaseq.data.hm <- rnaseq.data %>% filter(!is.na(CCLE_name)) %>%
                                    transmute(CCLE_name,
                                              ID=symbol,
                                              Original=as.character(RPKM),
                                              Value=pmax(log2(RPKM), -2),
                                              Type='rnaseq')
  
  #get rid of any duplicates
  rnaseq.data.hm <- rnaseq.data.hm[!duplicated(rnaseq.data.hm[,c('CCLE_name', 'ID')]) , ]
  
  #scale values for heatmap
  rnaseq.data.hm <- rnaseq.data.hm %>% group_by(ID) %>% mutate(zscore=scale(Value), zscore_offset=zscore + 60)

  #make the scalerange and colours
  rnaseq.data.scaleinfo <- data.frame(Type='rnaseq',
                                      Level=c('min', 'mean', 'max'),
                                      Value=c(min(rnaseq.data.hm$zscore_offset), mean(rnaseq.data.hm$zscore_offset), max(rnaseq.data.hm$zscore_offset)),
                                      Colour=c('blue', 'white', 'red'),
                                      stringsAsFactors=FALSE)
  
  return ( list (data_hm = ungroup(rnaseq.data.hm),
                 scaleinfo = rnaseq.data.scaleinfo))
  
}

get_affy <- function(params) {
  
  #get affy data for genes and cell lines of interest
  sql <- sprintf('select * from ccle_exprs_tall where Description IN (%s) and Tumor_Sample_Barcode IN (%s)', params$gene.symbols.sql, params$cls.ids.sql) 
  affy.data <- dbGetQuery(params$con, sql)
  
  #select and rename columns
  affy.data.hm <- affy.data %>% transmute(CCLE_name = Tumor_Sample_Barcode,
                                          ID = Description,
                                          Original = as.character(Signal),
                                          Value = Signal,
                                          Type = 'affy')
  
  #check for duplicates
  dup_check <- affy.data.hm %>% group_by(ID) %>% summarise(N=n()) %>% arrange(desc(N)) %>% filter(N == median(N))
  
  #scale values for heatmap
  affy.data.hm <- affy.data.hm %>% filter(ID %in% dup_check$ID) %>% group_by(ID) %>% 
    mutate(zscore=scale(Value), zscore_offset=zscore + 90)

  #make the scalerange and colours
  affy.data.scaleinfo <- data.frame(Type='affy',
                                      Level=c('min', 'mean', 'max'),
                                      Value=c(min(affy.data.hm$zscore_offset), mean(affy.data.hm$zscore_offset), max(affy.data.hm$zscore_offset)),
                                      Colour=c('blue', 'white', 'red'),
                                    stringsAsFactors=FALSE)
  
  return ( list (data_hm = ungroup(affy.data.hm),
                 scaleinfo = affy.data.scaleinfo))
  
}

get_hybcap <- function (params) {
  
  #get hybrid capture sequencing data for genes and cell lines of interest
  sql <- sprintf('select * from ccle_hybrid_capture where Hugo_symbol IN (%s) and Tumor_Sample_Barcode IN (%s)', params$gene.symbols.sql, params$cls.ids.sql)
  hybcap.data <- dbGetQuery(params$con,sql)
  
  #filter rows that have a protein change and select and rename columns
  hybcap.data.hm <- hybcap.data %>% filter(nchar(Protein_Change) > 1) %>%
                                    transmute(CCLE_name = Tumor_Sample_Barcode,
                                              ID = Hugo_Symbol,
                                              Original = Protein_Change)
  
  #ensure that only one row per gene/cell line and add value column
  hybcap.data.hm <- hybcap.data.hm %>% group_by(CCLE_name, ID) %>% summarise(Original=paste(Original, collapse=',')) %>% mutate(Value=1)
  
  #which samples were actually sequenced?
  hybcap.tested <- dbGetQuery(params$con, 'select distinct Tumor_Sample_Barcode from ccle_hybrid_capture')
  hybcap.tested <- hybcap.tested %>% filter(Tumor_Sample_Barcode %in% params$cls.ids)
  sequenced_ids <- hybcap.tested$Tumor_Sample_Barcode
  notsequenced_ids <- setdiff(params$cls.ids, sequenced_ids)
  
  #generate dataframes for sequenced and not sequenced cell lines
  hybcap.sequenced <- data.frame ( CCLE_name=rep(sequenced_ids, length(params$gene.symbols)) , 
                                   ID=rep(params$gene.symbols, each= length(sequenced_ids) ),
                                   Original='-',
                                   Value=0, stringsAsFactors=FALSE)
  
  if (length(notsequenced_ids) > 0 ) {
    
    hybcap.notsequenced <- data.frame ( CCLE_name=rep(notsequenced_ids, length(params$gene.symbols)) , 
                                        ID=rep(params$gene.symbols, each= length(notsequenced_ids) ),
                                        Original=NA,
                                        Value=NA, stringsAsFactors=FALSE )
    
  } else {
    hybcap.notsequenced <- data.frame ()
  }
  
  #get rid of rows in hybcap.sequenced which are duplicated in hybcap.data.hm
  hybcap.sequenced <- hybcap.sequenced %>% filter(!( paste(CCLE_name, ID) %in% paste(hybcap.data.hm$CCLE_name, hybcap.data.hm$ID) ))
  
  #combine hybcap dataframes and add additional standard columns
  hybcap.data.hm <- bind_rows(hybcap.data.hm, hybcap.sequenced, hybcap.notsequenced) %>%
                                  mutate(Type='hybcap', zscore=Value, zscore_offset=Value+20)
  
  #make the scalerange and colours
  hybcap.data.scaleinfo <- data.frame(Type='hybcap',
                                    Level=c('min', 'max'),
                                    Value=c(20,21),
                                    Colour=c('lightblue', 'darkblue'),
                                    stringsAsFactors=FALSE)
  
  return ( list (data_hm = hybcap.data.hm,
                 scaleinfo = hybcap.data.scaleinfo))
  
}

get_cosmicclp <- function(params) {
  
  #get cosmic cell line exome sequencing for genes of interest
  sql4 <- sprintf('select * from cosmicclp_mutations where Gene_name IN (%s)', params$gene.symbols.sql)  #make cosmic sql
  cosmic.data <- dbGetQuery(params$con, sql4) #get cosmic data
  cosmic.data.filtered <- subset(cosmic.data, cosmic.data$Sample_name %in% params$cls.ids.cosmic ) #filter by cell lines of interest
  cosmic.data.filtered <- subset(cosmic.data.filtered, !grepl('silent|Unknown', cosmic.data.filtered$Mutation_Description) ) #filter by coding change mutations
  cosmic.data.hm <- cosmic.data.filtered [  , c('Sample_name', 'Gene_name', 'Mutation_AA')  ] #just get fields of interest
  colnames(cosmic.data.hm) <- c('CCLE_name', 'ID', 'Original')
  
  #ensure that only one row per gene/cell line
  cosmic.data.hm <- ddply (cosmic.data.hm, .(CCLE_name, ID), function(x) c(Original=paste(x$Original, collapse=',')))
  cosmic.data.hm$Value <- 1
  
  #convert cosmic id back to CCLE id
  cosmic.data.hm$CCLE_name <-  params$allsampleinfo [ match(cosmic.data.hm$CCLE_name, params$allsampleinfo$cosmicclp_id) , 'CCLE_name']
  
  #which samples were actually seqenced?
  sequenced_ids <- params$allsampleinfo [ params$allsampleinfo$cosmicclp_id %in%  params$cls.ids.cosmic , 'CCLE_name']
  notsequenced_ids <-  setdiff(params$cls.ids, sequenced_ids )   

  #generate dataframes for sequenced and not sequenced cell lines
  cosmic.sequenced <- data.frame ( CCLE_name=rep(sequenced_ids, length(params$gene.symbols)) , 
                                   ID=rep(params$gene.symbols, each= length(sequenced_ids) ),
                                   Original='-',
                                   Value=0)

  
  if (length(notsequenced_ids) > 0 ) {
    
    cosmic.notsequenced <- data.frame ( CCLE_name=rep(notsequenced_ids, length(params$gene.symbols)) , 
                                        ID=rep(params$gene.symbols, each= length(notsequenced_ids) ),
                                        Original=NA,
                                        Value=NA )
    
  } else {
    cosmic.notsequenced <- data.frame ()
  }
  
  #get rid of rows in cosmic.sequenced which are duplicated in cosmic.data.hm
  cosmic.sequenced <- cosmic.sequenced [  !( paste(cosmic.sequenced$CCLE_name, cosmic.sequenced$ID) %in% paste(cosmic.data.hm$CCLE_name, cosmic.data.hm$ID) ) ,   ]
  
  
  #combine cosmic dataframes
  cosmic.data.hm <- rbind(cosmic.data.hm, cosmic.sequenced,cosmic.notsequenced)
  cosmic.data.hm$Type <- 'cosmic'
  #rescale to get zscores
  cosmic.data.hm$zscore <- cosmic.data.hm$Value
  cosmic.data.hm$zscore_offset <- cosmic.data.hm$Value +30
  
  #make the scalerange and colours
  cosmic.data.scaleinfo <- data.frame(Type='cosmic',
                                      Level=c('min', 'max'),
                                      Value=c(30,31),
                                      Colour=c('lightgreen', 'darkgreen'))
  
  return ( list (data_hm = cosmic.data.hm,
                 scaleinfo = cosmic.data.scaleinfo))
  
}

get_cn <- function (params) {
  
  #get copy number data for genes of interest
  sql1 <- sprintf('select * from ccle_copynumber_tall where geneName IN (%s) and Tumor_Sample_Barcode IN (%s)', params$gene.symbols.sql, params$cls.ids.sql) 
  cn.data <- dbGetQuery(params$con, sql1)
  
  #select and rename columns
  cn.data.hm <- cn.data %>% transmute(CCLE_name = Tumor_Sample_Barcode, 
                                      ID = geneName,
                                      Original = as.character(log2cn),
                                      Value = log2cn,
                                      Type = 'cn')
  
  #scale values for heatmap
  cn.data.hm <- cn.data.hm %>% group_by(ID) %>% mutate(zscore=scale(Value), zscore_offset=zscore + 120)
  
  #make the scalerange and colours
  cn.data.scalerange <- c(min(cn.data.hm$zscore_offset), 120, max(cn.data.hm$zscore_offset))
  cn.data.scalecols <- c('blue', 'white', 'red')
  
  #make the scalerange and colours
  cn.data.scaleinfo <- data.frame(Type='cn',
                                      Level=c('min', 'mean', 'max'),
                                      Value=c(min(cn.data.hm$zscore_offset), mean(cn.data.hm$zscore_offset), max(cn.data.hm$zscore_offset)),
                                      Colour=c('blue', 'white', 'red'))
  
  return ( list (data_hm = ungroup(cn.data.hm),
                 scaleinfo = cn.data.scaleinfo))
  
}

#get compound response data
get_response <- function (params) {
  
  if (params$respdata.source == 'custom') {
    
    #what to do if response data is null
    if (is.null(params$respdata)) {
      return( list (data_hm = data.frame(),
                    scaleinfo = data.frame()))
    }
    
    resp.data <- params$respdata
    #rename cols
    colnames(resp.data) <- c('CCLE_name', 'ID', 'Original', 'Value')
    
  }
  

  if (params$respdata.source == 'gdsc') {
    sql6 <- sprintf('select * from gdsc_screening where value_type = \'IC_50\' and compound IN (%s)', params$compounds.sql) 
    resp.data <- dbGetQuery(params$con, sql6)
    #just get fields of interest
    resp.data <- resp.data[,c('Cell_Line', 'compound', 'value')]
    colnames(resp.data) <- c('CCLE_name', 'ID', 'Original')
    resp.data$Value <- (log10(exp(resp.data$Original)) - 6 ) * -1
    #convert cosmic id back to CCLE id
    resp.data$CCLE_name <-  params$allsampleinfo [ match(resp.data$CCLE_name, params$allsampleinfo$gdsc_id) , 'CCLE_name']
    #get rid of data where cell line not in CCLE
    resp.data <- subset(resp.data, !is.na(resp.data$CCLE_name))
    #get rid of data where value is NA
    resp.data <- subset(resp.data, !is.na(resp.data$Value))
  }
  
  if (params$respdata.source == 'ccle') {
    sql6 <- sprintf('select * from ccle_drug_data where Compound IN (%s)', params$compounds.sql) 
    resp.data <- dbGetQuery(params$con, sql6)
    #just get fields of interest
    resp.data <- resp.data[,c('CCLE_Cell_Line_Name', 'Compound', 'ActArea')]
    colnames(resp.data) <- c('CCLE_name', 'ID', 'Original')
    resp.data$Value <- resp.data$Original

  }
  
  if (params$respdata.source == 'achilles') {
    sql6 <- sprintf('select * from achilles_v243 where shRNA_ID IN (%s)', params$compounds.sql) 
    resp.data <- dbGetQuery(params$con, sql6)
    #just get fields of interest
    resp.data <- resp.data[,c('CCLE_ID', 'shRNA_ID', 'log2FC')]
    colnames(resp.data) <- c('CCLE_name', 'ID', 'Original')
    resp.data$Value <- resp.data$Original * (-1)
    
  }
  
  #filter by specified cell lines
  resp.data.hm <- subset(resp.data, resp.data$CCLE_name %in% params$cls.ids)
  #define type
  resp.data.hm$Type <- 'resp'
  #rescale the data to get z scores
  resp.data.hm <- ddply(resp.data.hm, .(ID), transform, zscore=scale(Value))
  resp.data.hm$zscore_offset <- resp.data.hm$zscore - 10
  resp.data.hm$Original <- as.character(resp.data.hm$Original)
  
  #make the scalerange and colours
  resp.data.scaleinfo <- data.frame(Type='resp',
                                  Level=c('min', 'mean', 'max'),
                                  Value=c(min(resp.data.hm$zscore_offset, na.rm=TRUE), mean(resp.data.hm$zscore_offset, na.rm=TRUE), max(resp.data.hm$zscore_offset, na.rm=TRUE)),
                                  Colour=c('red', 'yellow', 'green'))
  
  return ( list (data_hm = resp.data.hm,
                 scaleinfo = resp.data.scaleinfo))
  
}

#get custom data
get_custom <- function (params) {
  
  #what to do if custom data is null
  if (is.null(params$customdata)) {
    return( list (data_hm = data.frame(),
                  scaleinfo = data.frame()))
  }
  
  cust.data <- params$customdata
  #rename cols
  colnames(cust.data) <- c('CCLE_name', 'ID', 'Original', 'Value')
  #filter by specified cell lines
  cust.data.hm <- subset(cust.data, cust.data$CCLE_name %in% params$cls.ids)

  #define type
  cust.data.hm$Type <- 'custom'
  
  #rescale the data to get z scores
  cust.data.hm <- ddply(cust.data.hm, .(ID), transform, zscore=scale(Value))
  cust.data.hm$zscore_offset <- cust.data.hm$zscore + 160
  
  #make the scalerange and colours
  cust.data.scaleinfo <- data.frame(Type='custom',
                                    Level=c('min', 'mean', 'max'),
                                    Value=c(min(cust.data.hm$zscore_offset, na.rm=TRUE), mean(cust.data.hm$zscore_offset, na.rm=TRUE), max(cust.data.hm$zscore_offset, na.rm=TRUE)),
                                    Colour=c('blue', 'white', 'red'))
  
  return ( list (data_hm = cust.data.hm,
                 scaleinfo = cust.data.scaleinfo))
  
}

make_heatmap <- function (params, datatype=c('resp' , 'rnaseq', 'affy', 'hybcap', 'cosmic', 'cn', 'custom'), order_by='resp', compound=1, multi_cpds=FALSE) {
  
  #get the data
  resp_data <- get_response (params)
  
  if ('rnaseq' %in% datatype) {
    rnaseq_data <- get_rnaseq (params)
  } else { rnaseq_data <- NULL }
  
  if ('affy' %in% datatype) {
    affy_data <- get_affy (params)
  } else { affy_data <- NULL }
  
  if ('hybcap' %in% datatype) {
    hybcap_data <- get_hybcap (params)
  } else { hybcap_data <- NULL }
  
  if ('cosmic' %in% datatype) {
    cosmicclp_data <- get_cosmicclp (params)
  } else { cosmicclp_data <- NULL }
  
  if ('cn' %in% datatype) {
    cn_data <- get_cn (params)
  } else { cn_data <- NULL }
  
  if ('custom' %in% datatype) {
    custom_data <- get_custom (params)
  } else { custom_data <- NULL }
  
  #combine
  plot.data <- rbind(resp_data$data_hm, hybcap_data$data_hm , cosmicclp_data$data_hm , rnaseq_data$data_hm , affy_data$data_hm, cn_data$data_hm, custom_data$data_hm)
  
  #filter out unneeded
  plot.data <- subset(plot.data , plot.data$Type %in% datatype)
  
  #get compound to order by
  if (is.numeric(compound)) {
    
    if (compound <= length(params$compounds)) {
      compound_name <- params$compounds[compound]
    } else {
      stop(sprintf('There are only %s compounds but you selected compound no %s', length(params$compounds), compound))
    }
    
    
  } else {
    
    if (compound %in% params$compounds) {
      compound_name <- compound
    } else {
      stop(sprintf('Compound %s needs to be specified in params', compound))
    }
    
  }
  
  #get rid of unneeded response data if multi_cpds is FALSE
  if (multi_cpds == FALSE) {
    plot.data <- subset(plot.data, !(plot.data$Type == 'resp' & plot.data$ID != compound_name)  )
  }
  
  #featurename order
  plot.data$FeatureName <- paste(plot.data$Type, plot.data$ID, sep='_')
  plot.data$Type <- factor(plot.data$Type, levels=datatype)
  ordered_feature_names <- unique(plot.data[ order(plot.data$Type, plot.data$ID) ,c('FeatureName')])
  if (order_by == 'resp') {
    order_data <- subset(resp_data$data_hm, resp_data$data_hm$ID == compound_name)
    ordered_cell_lines <- order_data [ order(order_data$Value, decreasing=T), 'CCLE_name'  ]
  } else {
    ordered_cell_lines <- sort(unique(plot.data$CCLE_name), decreasing=T) 
  }
  
  #get original scaleinfo 
  scaleinfo <- rbind(resp_data$scaleinfo, hybcap_data$scaleinfo , cosmicclp_data$scaleinfo, rnaseq_data$scaleinfo, affy_data$scaleinfo, cn_data$scaleinfo, custom_data$scaleinfo)
  
  #make new scaleinfo
  #scaleinfo_new <- scaleinfo
  scaleinfo_new <- ddply ( plot.data , .(Type) , function (x) c(min=min(x$zscore_offset, na.rm=TRUE), mean=mean(x$zscore_offset, na.rm=TRUE), max=max(x$zscore_offset, na.rm=TRUE)) )
  scaleinfo_new <- melt(scaleinfo_new, id.vars='Type')
  colnames(scaleinfo_new) <- c('Type', 'Level', 'Value') 
  
  #get colours
  scaleinfo_new$Colour <- scaleinfo [ match( paste(scaleinfo_new$Type,scaleinfo_new$Level) , paste(scaleinfo$Type,scaleinfo$Level) ) ,'Colour'  ]
  
  #order and get rid of nas
  scaleinfo_new <- subset(scaleinfo_new, !is.na(scaleinfo_new$Colour))
  scaleinfo_new <- scaleinfo_new [ order(scaleinfo_new$Value) , ]
  
  p <- ggplot(plot.data, aes(x=CCLE_name, y=FeatureName)) 
  p <- p + geom_tile(aes(fill = zscore_offset), linetype=0 ) 
  p <- p + scale_fill_gradientn(colours = scaleinfo_new$Colour, values = rescale(scaleinfo_new$Value))
  p <- p + scale_y_discrete(limits=ordered_feature_names) #this orders the y axis as we want it
  p <- p + scale_x_discrete(limits=ordered_cell_lines) #this orders the x axis as we want it
  p <- p + coord_flip()
  p <- p + theme_bw(base_size = 9) 
  p <- p + labs(x = "", y = "")
  p <- p + theme(axis.text.x = element_text(size=rel(1), angle=330, hjust=0, vjust=1), axis.text.y = element_text(size=rel(1)))
  p <- p + theme(panel.background=element_rect(fill="gray90", color='white'), legend.position = "none")
  
  
  return(p)
  
  
  
}

make_matrices <- function (params, datatype=c('rnaseq', 'affy', 'hybcap', 'cosmic', 'cn', 'custom'))  {
  
  #get the data
  if ('rnaseq' %in% datatype) {
    rnaseq_data <- get_rnaseq (params)
  } else { rnaseq_data <- NULL }
    
  if ('affy' %in% datatype) {
    affy_data <- get_affy (params)
  } else { affy_data <- NULL }
  
  if ('hybcap' %in% datatype) {
    hybcap_data <- get_hybcap (params)
  } else { hybcap_data <- NULL }
  
  if ('cosmic' %in% datatype) {
    cosmicclp_data <- get_cosmicclp (params)
  } else { cosmicclp_data <- NULL }
  
  if ('cn' %in% datatype) {
    cn_data <- get_cn (params)
  } else { cn_data <- NULL }
  
  if ('custom' %in% datatype) {
    custom_data <- get_custom (params)
  } else { custom_data <- NULL }
  
  #combine
  all_data <- rbind(hybcap_data$data_hm , cosmicclp_data$data_hm , rnaseq_data$data_hm , affy_data$data_hm, cn_data$data_hm, custom_data$data_hm)
  
  #create matrix
  featureData <- acast (all_data, CCLE_name ~ ID + Type, value.var='Value', fun.aggregate=first) 
  #x <- all_data %>% transmute(CCLE_name, fn=paste(ID,Type), Value) %>% spread(fn, Value)
  
  #get the response data
  resp_data <- get_response (params)
  
  #create matrix
  responseData <- acast (resp_data$data_hm , CCLE_name ~ ID , value.var='Value', fun.aggregate=first)
  
  #order matrix columns so that the same order as params$compounds
  check_cpds <- setdiff(params$compounds, colnames(responseData))
  if (length(check_cpds) == 0) {
    responseData <- responseData [,params$compounds]
  } else {
    warning (sprintf('WARNING: Data not found for the following compounds: %s', paste(check_cpds, collapse=', ')))
  }
    
  #ensure that same cell lines present in both matrices
  notinboth <- c ( setdiff(rownames(featureData),rownames(responseData)) , setdiff(rownames(responseData),rownames(featureData))  )
  inboth <- intersect(rownames(featureData),rownames(responseData))
  
  if (length(notinboth) > 0) {
    warning (sprintf('The following cell lines dont have response and feature data: %s', paste(notinboth, collapse=',')))
    featureData <- featureData [ inboth ,  ]
    responseData <- responseData [inboth , ]
  }
  
  return (list( featureData=featureData,
                responseData=responseData))
  
}

make_caret_df <- function (params, data_mat, cpd=1, filterNZV=TRUE) {
  
  caret_df <- data.frame(resp=data_mat$responseData[,params$compounds[cpd]], data_mat$featureData)  #make a data frame out of the matrix for one cpd
  factor_cols <- grepl('_hybcap|_cosmic', colnames(caret_df))  #get 0/1 variables 
  for (i in which(factor_cols)) {caret_df[,i] <- as.factor(caret_df[,i])} #turn 0/1 variables into factors
  
  #get rid of near zero variance variables
  if (filterNZV == TRUE) {
    nzv <- nearZeroVar(caret_df)
    caret_df <- caret_df[, -nzv]  
  }
  
  return(caret_df)
  
}

make_df <- function (params, datatype=c('rnaseq', 'affy', 'hybcap', 'cosmic', 'cn', 'custom'))  {
  
  #make a data container list
  all_data.list <- list()
  
  #get the data
  if ('rnaseq' %in% datatype) {
    rnaseq_data <- get_rnaseq (params)
    all_data.list[['rnaseq_data']] <- rnaseq_data$data_hm
  } 
  
  if ('affy' %in% datatype) {
    affy_data <- get_affy (params)
    all_data.list[['affy_data']] <- affy_data$data_hm
  } 
  
  if ('hybcap' %in% datatype) {
    hybcap_data <- get_hybcap (params)
    all_data.list[['hybcap_data']] <- hybcap_data$data_hm
  } 
  
  if ('cosmic' %in% datatype) {
    cosmicclp_data <- get_cosmicclp (params)
    all_data.list[['cosmicclp_data']] <- cosmicclp_data$data_hm
  } 
  
  if ('cn' %in% datatype) {
    cn_data <- get_cn (params)
    all_data.list[['cn_data']] <- cn_data$data_hm
  } 
  
  if ('custom' %in% datatype) {
    custom_data <- get_custom (params)
    all_data.list[['custom_data']] <- custom_data$data_hm
  } 
  
  #get the response data
  resp_data <- get_response (params)
  all_data.list[['resp_data']] <- resp_data$data_hm
  
  #check that all compounds in params have been returned, give warning if not
  missing_resps <- setdiff(params$compounds,unique(resp_data$data_hm$ID))
  if (length(missing_resps) != 0) {
    warning(sprintf('No data returned for: %s', paste(missing_resps, collapse=', ')))
    
  }
  
  #combine
  all_data <- bind_rows(all_data.list)
  
  #get rid of cell lines with no response data
  present_cls <- unique(resp_data$data_hm$CCLE_name)
  missing_cls <- setdiff(params$cls.ids, present_cls)
  all_data <- all_data %>% filter(CCLE_name %in% present_cls)
  if (length(missing_cls) != 0) {
    warning(sprintf('No response data for following cell lines: %s', paste(missing_cls, collapse=', ')))
    
  }
  
  #create matrix
  output.df <- all_data %>% transmute(CCLE_name, fn=paste(ID,Type,sep='_'), Value) %>% spread(fn, Value)
  
  #reorder columns
  output.df <- output.df %>% dplyr::select(CCLE_name, ends_with('_resp'), everything())
  
  #make mutation fields into factors
  output.df <- output.df %>% mutate_each(funs(as.factor), ends_with('_hybcap|_cosmic'))
      
  return(output.df)

  
  
}

run_glmnet <- function (params, the_df, cpds=NULL, lambda=NULL, alpha=NULL) {
 
  require(caret);require(foreach);require(parallel);require(doParallel)
  
  result_list <- list()
  #check that there are data in the matrix for the cpd id's requested
  #if (length(setdiff(cpds, 1:ncol(matrices$responseData)) ) > 0 ) {stop('Not all compounds present in supplied dataset')}
  
  #get rid of near zero variance variables
  nzv <- nearZeroVar(the_df)
  the_df <- the_df[, -nzv]  

  
  #allow for custom values of lambda and alpha
  if(is.null(lambda)) {lambda <- 10^seq(8,-5, length=100)}
  if(is.null(alpha)) {alpha <- seq(0,1,0.05)}
  
  for (i in cpds) {
    cpd_name <- params$compounds[i]
    print(sprintf('Fitting glmnet for compound %s: %s', i, cpd_name))
    caret_data <- the_df %>% dplyr::select(everything(), -ends_with('resp'), contains(cpd_name), -CCLE_name)
    ctrl <- trainControl(method = "repeatedcv",number=5, repeats = 5 )  #5 fold x validation
    glmnet_grid <- expand.grid(lambda=lambda, alpha = alpha) #grid of parameters to test 
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    glmnetFit <- train(formula(sprintf("`%s_resp` ~ .", cpd_name)),
                       data = caret_data,
                       method='glmnet',
                       tuneGrid=glmnet_grid,
                       preProcess = c("center", "scale"),
                       trControl = ctrl
                       
    )
    stopCluster(cl)
    result_list[[i]] <- glmnetFit
    names(result_list)[i] <- cpd_name
  }
  
  return(result_list)
  
}

run_gbm <- function (params, the_df, cpds=NULL, n.trees=NULL, interaction.depth=NULL, shrinkage=NULL)  {
  
  require(caret);require(foreach);require(parallel);require(doParallel)
  
  result_list <- list()
  #check that there are data in the matrix for the cpd id's requested
  #if (length(setdiff(cpds, 1:ncol(matrices$responseData)) ) > 0 ) {stop('Not all compounds present in supplied dataset')}
  
  #get rid of near zero variance variables
  nzv <- nearZeroVar(the_df)
  the_df <- the_df[, -nzv]  
  
  #allow for custom values of lambda and alpha
  if(is.null(n.trees)) { n.trees <- (1:30)*15 }
  if(is.null(interaction.depth)) { interaction.depth <- c(1,2,3) }
  if(is.null(shrinkage)) { shrinkage <- c(0.003, 0.01) }
  
  for (i in cpds) {
    cpd_name <- params$compounds[i]
    print(sprintf('Fitting boosted trees for compound %s: %s', i, cpd_name))
    caret_data <- the_df %>% dplyr::select(everything(), -ends_with('resp'), contains(cpd_name), -CCLE_name)
    ctrl <- trainControl(method = "repeatedcv",number=5, repeats = 5 )  #5 fold x validation
    gbmFit_grid <- expand.grid(n.trees=n.trees, interaction.depth=interaction.depth, shrinkage=shrinkage) #grid of parameters to test
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    gbmFit <- train(formula(sprintf("`%s_resp` ~ .", cpd_name)),
                       data = caret_data,
                       method='gbm',
                       tuneGrid=gbmFit_grid,
                       trControl = ctrl
                       
    )
    stopCluster(cl)
    result_list[[i]] <- gbmFit
    names(result_list)[i] <- cpd_name
  }
  
  return(result_list)
  
  
}

process_betas <- function (res, cpd_name) {
  
  if ('train' %in% class(res)) {
    #in this case we have just been given a single caret train object so leave res as is
    print('single train object') 
    warning('Note that in single train object mode no check is done on compound being correct')
  } else if ('list' %in% class(res) & 'train' %in% class(res[[cpd_name]])) {
    print('list of train objects')
    #in this case we have a list of train objects so take the one we need or stop if it isn't there
    if(!(cpd_name %in% names(res))) {stop('No results for supplied compounds')}
    res <- res[[cpd_name]]
  } else if ('list' %in% class(res) & 'matrix' %in% class(res[[1]])) {
    print('rr_list object')
    #in this case we have an rr_list object so need to handle it differently and return output
    output <- res$Betas.bma[[cpd_name]]
    if(is.null(output)) {stop('No results for supplied compounds')}
    output <- as.data.frame(output)
    return(output)
  } else {
    stop('uknown result type supplied')
  }
  
  if (res$method == 'glmnet') {
    #varImp doesn't work very well for glmnet objects so extract the betas ourselves 
    bestModel <- res$finalModel
    bestLambda_idx <- which.min(abs(bestModel$lambda - bestModel$lambdaOpt)) 
    output <- as.data.frame(bestModel$beta [ , bestLambda_idx ])
  } else {
    #most of the time just use the varImp function
    output <- varImp(res, scale = FALSE)$importance
  }
  
  return(output)
  
}

make_rr_heatmap <- function (params, res, cpd, n=20) {
  
  require(caret)
  require(dplyr)
  require(ggplot2)
  
  #resolve the provided compound
  if(is.numeric(cpd) & cpd <= length(params$compounds)) {
    cpd_name <- params$compounds[cpd]
  } else if ( is.character(cpd) & cpd %in% params$compounds ) {
    cpd_name <- cpd
  } else {
    stop('The cpd parameter should either be a number corresponding to the compound in params$compounds or a string matching a compound in params$compounds')
  }
  
  #work out what res is and extract importance/beta values
  betas <- process_betas (res, cpd_name)
  
  #process the betas data frame
  bn <- nrow(betas)
  colnames(betas)[1] <- 'beta_value' 
  betas$beta_name <- gsub('1$', '', rownames(betas))  #get rid of anything appended to the end of cosmic or hybcap
  betas <- betas %>% dplyr::arrange(desc(abs(beta_value))) %>% dplyr::mutate(idx=1:bn)
  
  #generte plot.data dataframe which is the top n betas
  plot_betas <- betas %>% filter(idx <= n) %>% arrange(desc(beta_value))
  
  #get symbols that we'll need and modify params so that we don't have to get data for everything
  params$gene.symbols <- unique(gsub('_hybcap|_cosmic|_custom|_affy|_rnaseq|_cn','',plot_betas$beta_name))
  params$gene.symbols.sql <- paste("'", params$gene.symbols, "'", sep='', collapse=',')
  params$gene.ids <- symbolToGene(params$gene.symbols, as.vector=T)
  
  #get the feature data
  rnaseq_data <- get_rnaseq (params)
  affy_data <- get_affy (params)
  hybcap_data <- get_hybcap (params)
  cosmicclp_data <- get_cosmicclp (params)
  cn_data <- get_cn (params)
  custom_data <- get_custom (params)
  
  #combine
  feature.data <- rbind(hybcap_data$data_hm , cosmicclp_data$data_hm , rnaseq_data$data_hm , affy_data$data_hm, cn_data$data_hm, custom_data$data_hm)
  
  #construct a beta name and  filter so that we only have the ones in the beta list
  feature.data <- feature.data %>% mutate(beta_name = paste(ID,Type,sep='_')) %>% filter(beta_name %in% plot_betas$beta_name)
  
  #now get the response data for just this cpd and construct the beta name
  resp_data <- get_response (params)
  resp.data.filtered <- resp_data$data_hm %>% filter(ID == cpd_name & !is.na(Value)) %>% mutate(beta_name = paste(ID,Type,sep='_')) %>% arrange(desc(Value))
      
  #get the cell line ids in the resp data
  respdata_cls <- resp.data.filtered$CCLE_name
  if(anyDuplicated(respdata_cls)>0) {warning('Duplicate cell lines present - this may cause issues')}
  
  #combine data
  plot.data <- rbind(feature.data , resp.data.filtered)
  
  #finally fiter the plot data by the cell lines in the response data
  plot.data <- plot.data %>% filter(CCLE_name %in% respdata_cls)
  
  #get original scaleinfo 
  scaleinfo <- rbind(resp_data$scaleinfo, hybcap_data$scaleinfo , cosmicclp_data$scaleinfo, rnaseq_data$scaleinfo, affy_data$scaleinfo, cn_data$scaleinfo, custom_data$scaleinfo)
  scaleinfo <- scaleinfo %>% mutate(Type=as.character(Type), Level=as.character(Level), Colour=as.character(Colour)) %>% dplyr::select(-Value)
  
  #make new scaleinfo
  scaleinfo_new <- plot.data %>% group_by(Type) %>% summarise(min=min(zscore_offset, na.rm=TRUE), mean=mean(zscore_offset, na.rm=TRUE), max=max(zscore_offset, na.rm=TRUE)) 
  scaleinfo_new <- melt(as.data.frame(scaleinfo_new), id.vars='Type') %>% mutate(variable=as.character(variable))
  colnames(scaleinfo_new) <- c('Type', 'Level', 'Value') 
  
  #get colours from scaleinfo
  scaleinfo_new <- scaleinfo_new %>% inner_join (scaleinfo, by=c('Type'='Type', 'Level'='Level')) %>% arrange(Value)
  
  #order beta names including resp beta name
  ordered_beta_names <- c(unique(resp.data.filtered$beta_name) , plot_betas$beta_name)
  
  p <- ggplot(plot.data, aes(x=CCLE_name, y=beta_name)) 
  p <- p + geom_tile(aes(fill = zscore_offset), linetype=0 ) 
  p <- p + scale_fill_gradientn(colours = scaleinfo_new$Colour, values = rescale(scaleinfo_new$Value))
  p <- p + scale_y_discrete(limits=ordered_beta_names) #this orders the y axis as we want it
  p <- p + scale_x_discrete(limits=respdata_cls) #this orders the x axis as we want it
  p <- p + coord_flip()
  p <- p + theme_bw(base_size = 9) 
  p <- p + labs(x = "", y = "")
  p <- p + theme(axis.text.x = element_text(size=rel(1), angle=330, hjust=0, vjust=1), axis.text.y = element_text(size=rel(1)))
  p <- p + theme(panel.background=element_rect(fill="gray90", color='white'), legend.position = "none")
  p <- p + ggtitle(sprintf('Top %s Features for %s', n, cpd_name))
  
  return(p)
}

make_rr_betaplot <- function (params, res, cpd, n=100) {
  
  require(caret)
  require(dplyr)
  require(ggplot2)
  
  #resolve the provided compound
  if(is.numeric(cpd) & cpd <= length(params$compounds)) {
    cpd_name <- params$compounds[cpd]
  } else if ( is.character(cpd) & cpd %in% params$compounds ) {
    cpd_name <- cpd
  } else {
    stop('The cpd parameter should either be a number corresponding to the compound in params$compounds or a string matching a compound in params$compounds')
  }
  
  #work out what res is and extract importance/beta values
  betas <- process_betas (res, cpd_name)
  
  #process the betas data frame
  bn <- nrow(betas)
  colnames(betas)[1] <- 'beta_value' 
  betas$beta_name <- rownames(betas)
  betas <- betas %>% dplyr::arrange(desc(abs(beta_value))) %>% dplyr::mutate(idx=1:bn)
  
  #generte plot.data dataframe which is the top n betas
  plot.data <- betas %>% filter(idx <= n) %>% arrange(beta_value)
    
  #do the plot
  p <- ggplot ( plot.data , aes(y=beta_value, x=beta_name))
  p <- p + scale_x_discrete(limits= plot.data$beta_name)
  p <- p + coord_flip()
  p <- p + theme_bw()
  p <- p + geom_point(colour="black", size = 6)
  p <- p + geom_point(size=5, aes(colour=beta_value))
  p <- p + scale_color_gradient2(low='blue', mid='white', high='red')
  p <- p + theme(axis.text.x = element_text(size=rel(1), angle=330, hjust=0, vjust=1), axis.text.y = element_text(size=rel(1)))
  p <- p + ggtitle(sprintf('Top %s Features for %s', nrow(plot.data), cpd_name))
  #print(p)
  return(p)
  
}

getGenesOfBioProcess <- function ( proc  , mode='gene', gsc_file='/Users/pchapman/BigData/GSEA/c5.bp.v4.0.symbols.gmt' ) {
  
  #get genes for a specific gene set using piano
  
  #can either load a given filename or can use an already existing GSC object
  if (is.character(gsc_file)) {
    myGsc <- loadGSC(gsc_file)
  } else {
    myGsc <- gsc_file
  }
  
  proc <- as.character(proc)
  gs.symbols <- unlist(myGsc$gsc[proc])
  
  if (mode == 'symbol') {
    return(gs.symbols)
  } else {
    #convert to gene ids
    gs.genes <- symbolToGene(gs.symbols, as.vector=T)
    return(gs.genes)
    
  }
  
  
  
  
}

run_MANOVA_single <- function (control.df, data_mat, i) {
    
    cpd_id <- control.df [ i , 'cpd_id'  ]
    feature_id <- control.df [ i , 'feature_id'  ]
    
    resp <- data_mat$responseData[,cpd_id]
    feat <- data_mat$featureData[,feature_id]
    
    #do anova
    res <- aov (resp ~ feat)
    pval <- summary(res)[[1]][[5]][1]
    Fval <- summary(res)[[1]][[4]][1]
    
    #get diff between grps
    diffs <- tapply(resp, feat, mean, na.rm=TRUE)
    
    #calc effect size accounting for grps with one entry returned
    if (length (diffs) == 2) {
      effect_size <- diffs[1] - diffs[2]
    } else {
      effect_size <- NA
    }
    
    #size of groups
    mut_count <- length ( which ( feat == 1) )
    
    
    
    out <- data.frame(feature_name=colnames(data_mat$featureData)[feature_id] , 
                      cpd_name=colnames(data_mat$responseData)[cpd_id],
                      cpd_id=cpd_id,
                      Fval=as.numeric(Fval),
                      pval=as.numeric(pval),
                      effect_size=effect_size,
                      mut_count=mut_count,
                      stringsAsFactors=FALSE)
    return(out)

}

run_MANOVA <- function ( params, mat=NULL ) {
  
  require(parallel)
  require(doParallel)
  require(foreach)
  
  if (is.null(mat)) {
    mat <- make_matrices ( params , c('hybcap', 'cosmic') )
  }
  
  #make a data frame of all of the comparisons
  cpd_n <- ncol(mat$responseData)
  feature_n <- ncol(mat$featureData)
  control_df <- data.frame ( cpd_id = rep(1:cpd_n, times=feature_n),
                             feature_id = rep(1:feature_n, each=cpd_n))
  
  #results <- lapply (1:nrow(control.df) , run_MANOVA_single(x))
  
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  results_list <- foreach(j=1:nrow(control_df), .export='run_MANOVA_single') %dopar% run_MANOVA_single(control_df, mat, j)
  stopCluster(cl)
  
  #combine results into single df
  results <- rbind_all(results_list)
  
  #correct for multiple testing BY GENES (not by cpds)
  results <- ddply ( results ,  .(cpd_name) , transform ,  pval_fdr=p.adjust(pval , method='BH')  )
  
  #log FDR
  results$log10FDR <- -log10(results$pval_fdr)
  
  return(results)
  
}

#define reverse log10 scale
library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

make_manova_volcanoplot <- function ( params, manova.df, cpd=NULL, effect_th = 0.5, log10FDR_th = 2 ) {
  
  plot.data <- manova.df %>% filter(cpd_id == cpd)
  cpd_name <- unique(plot.data$cpd_name)
  
  sensitiser.data <- plot.data %>% filter (effect_size < -effect_th &  log10FDR > log10FDR_th)
  insensitiser.data <- plot.data %>% filter (effect_size > effect_th &  log10FDR > log10FDR_th)
  nonsig.data <- plot.data %>% filter (log10FDR <= log10FDR_th)
  
  
  
  
  xlim_param <- max(1,abs(plot.data$effect_size), na.rm=T)
  
  p <- ggplot ( plot.data , aes (x=effect_size, y=log10FDR, size=log10(mut_count)))
  p <- p + geom_point(data = nonsig.data, shape=21, fill='lightgray')
  p <- p + geom_point(data = sensitiser.data, shape=21, fill='lightgreen')
  p <- p + geom_point(data = insensitiser.data, shape=21, fill='orange')
  p <- p + theme_bw() + theme(legend.position="none")
  p <- p + ggtitle(sprintf('MANOVA for %s' , cpd_name))
  #do labels
  p <- p + geom_text(data = sensitiser.data, aes(label=feature_name), size=4, colour='darkgreen', hjust=0, vjust=0, angle=0)
  p <- p + geom_text(data = insensitiser.data, aes(label=feature_name), size=4, colour='red', hjust=0, vjust=0, angle=0)
  #do lines
  p <- p + geom_vline(xintercept = c(-effect_th,effect_th) , color='blue', alpha=0.5, linetype='dotted') + geom_hline(yintercept = log10FDR_th, color='blue', alpha=0.5, linetype='dotted')
  p <- p + xlim(c(xlim_param*(-1),xlim_param))
  
  print(p)
  
}

make_manova_heatmap <- function (params, manova.df, cpd=NULL, effect_th = 0.5, log10FDR_th = 2 ) {
  
  #sort out the manova data
  cpd.data <- subset(manova.df , manova.df$cpd_id == cpd)
  cpd_name <- as.character(unique(cpd.data$cpd_name))
  cpd.data <- subset(cpd.data , abs(cpd.data$effect_size) > effect_th &   cpd.data$log10FDR > log10FDR_th)
  cpd.data <- cpd.data [ order(cpd.data$log10FDR, decreasing=T)  ,  ]
  
  #get symbols that we'll need and modify params so that we don't have to get data for everything
  params$gene.symbols <- unique(gsub('_hybcap|_cosmic|_custom','',cpd.data$feature_name))
  params$gene.symbols.sql <- paste("'", params$gene.symbols, "'", sep='', collapse=',')
  params$gene.ids <- symbolToGene(params$gene.symbols, as.vector=T)
  
  
  #get the feature data
  hybcap_data <- get_hybcap (params)
  cosmicclp_data <- get_cosmicclp (params)
  custom_data <- get_custom (params)
  
  #combine
  feature.data <- rbind(hybcap_data$data_hm , cosmicclp_data$data_hm, custom_data$data_hm)
  
  #now filter so that we only have the ones in the beta list
  feature.data$feature_name  <- paste(feature.data$ID, feature.data$Type, sep='_')
  feature.data <- subset(feature.data , feature.data$feature_name %in% cpd.data$feature_name)
  
  #now get the response data for just this cpd
  resp_data <- get_response (params)
  resp.data.filtered  <- subset ( resp_data$data_hm , resp_data$data_hm$ID == cpd_name & !is.na(resp_data$data_hm$Value)  )
  resp.data.filtered$feature_name  <- paste(resp.data.filtered$ID, resp.data.filtered$Type, sep='_')
  
  #combine data
  plot.data <- rbind(feature.data , resp.data.filtered)
  
  #get original scaleinfo 
  scaleinfo <- rbind(resp_data$scaleinfo, hybcap_data$scaleinfo , cosmicclp_data$scaleinfo, custom_data$scaleinfo)
  
  #make new scaleinfo
  scaleinfo_new <- ddply ( plot.data , .(Type) , function (x) c(min=min(x$zscore_offset, na.rm=TRUE), mean=mean(x$zscore_offset, na.rm=TRUE), max=max(x$zscore_offset, na.rm=TRUE)) )
  scaleinfo_new <- melt(scaleinfo_new, id.vars='Type')
  colnames(scaleinfo_new) <- c('Type', 'Level', 'Value') 
  
  #get colours
  scaleinfo_new$Colour <- scaleinfo [ match( paste(scaleinfo_new$Type,scaleinfo_new$Level) , paste(scaleinfo$Type,scaleinfo$Level) ) ,'Colour'  ]
  
  #order and get rid of nas
  scaleinfo_new <- subset(scaleinfo_new, !is.na(scaleinfo_new$Colour))
  scaleinfo_new <- scaleinfo_new [ order(scaleinfo_new$Value) , ]
  
  #order feature names including resp feature name
  ordered_feature_names <- c(unique(resp.data.filtered$feature_name) , as.character(cpd.data$feature_name))
  
  #filter by cell lines that exist in response data
  respdata_cls <- unique(resp.data.filtered$CCLE_name) 
  plot.data <- subset(plot.data, plot.data$CCLE_name %in% respdata_cls)
  
  p <- ggplot(plot.data, aes(x=CCLE_name, y=feature_name)) 
  p <- p + geom_tile(aes(fill = zscore_offset), linetype=0 ) 
  p <- p + scale_fill_gradientn(colours = scaleinfo_new$Colour, values = rescale(scaleinfo_new$Value))
  p <- p + scale_y_discrete(limits=ordered_feature_names) #this orders the y axis as we want it
  p <- p + scale_x_discrete(limits=resp.data.filtered [ order(resp.data.filtered$Value, decreasing=T) , 'CCLE_name' ]) #this orders the x axis as we want it
  p <- p + coord_flip()
  p <- p + theme_bw(base_size = 9) 
  p <- p + labs(x = "", y = "")
  p <- p + theme(axis.text.x = element_text(size=rel(1), angle=330, hjust=0, vjust=1), axis.text.y = element_text(size=rel(1)))
  p <- p + theme(panel.background=element_rect(fill="gray90", color='white'), legend.position = "none")
  p <- p + ggtitle(sprintf('Top MANOVA features for %s',  cpd_name))
  
  
  return(p)
}

do_simple_analysis <- function (mat, resp, feature, feature_type='continuous', resp_type='continuous', resp_th='median') {
  
  #get data
  resp_data <- mat$responseData[,resp]
  feature_data <- mat$featureData[,feature]
  
  #do some error checking
 if ( !isTRUE(all.equal(names(feature_data), names(resp_data)) ) ) {
   print(cbind(names(feature_data), names(resp_data)))
   stop('different cell lines in response and feature data')
 }
 
 #makes feature data into a factor if required
 if (feature_type == 'discrete') {
   feature_data <- as.factor(feature_data)
 }
 
 #makes the response data into a factor if required
 if (resp_type == 'discrete') {
   
   if (resp_th == 'median') {
     resp_median <- median(resp_data, na.rm=TRUE)
     resp_data [ resp_data > resp_median  ] <- 'sensitive'
     resp_data [ resp_data <= resp_median  ] <- 'resistant'
     print(resp_median)
   } else if (is.numeric(resp_th)) {
     resp_data_l <- resp_data > resp_th
     resp_data [ resp_data_l  ] <- 'sensitive'
     resp_data [ !resp_data_l  ] <- 'resistant'
   } else {
     stop('resp_th must be median or a number')
   }
   
   resp_data <- as.factor(resp_data)
   
 }
 
 #get feature and compound names
 if (is.numeric(resp)) {
   resp_name <- colnames(mat$responseData)[resp]
 } else {
   resp_name <- resp
 }
 
 if (is.numeric(feature)) {
   feature_name <- colnames(mat$featureData)[feature]
 } else {
   feature_name <- feature
 }
 
 #make plot data
 plot_data <- data.frame(f=feature_data, r=resp_data, cell_line=names(resp_data), stringsAsFactors=F)
 
 #get rid of NA's
 plot_data <- subset (plot_data , !is.na(plot_data$f) & !is.na(plot_data$r)  )
  
 #now do a ggplot
 
 #normal correlation
 if (feature_type=='continuous' & resp_type=='continuous') {
   p <- ggplot(plot_data, aes(f,r) )
   p <- p + theme_bw() + scale_fill_gradient2(low='red', mid='yellow', high='green', midpoint=mean(range(plot_data$r)))
   p <- p + scale_y_reverse()
   p <- p + labs(x=feature_name, y=resp_name)
   p <- p + geom_point(aes(fill=r), size=7, shape=21)
   p <- p + geom_smooth(method=lm)
   
   #do an ANOVA
   stat_test <- aov(f~r, plot_data)
   stat_text <- capture.output(summary(stat_test))
   stat_text <- stat_text[1:3]
   stat_name <- 'ANOVA'
 }

 #sensitive vs resistant on x axis
 if (feature_type=='continuous' & resp_type=='discrete') {
   p <- ggplot(plot_data, aes(r,f) )
   p <- p + theme_bw() + scale_fill_brewer(palette="Set1") + scale_colour_brewer(palette="Set1")
   p <- p + labs(y=feature_name, x=resp_name)
   p <- p + geom_boxplot(aes(colour=r), outlier.size=0) #want a boxplot as well as points
   p <- p  + geom_point(aes(fill=r), size=7, shape=21, position = position_jitter(w = 0.2, h = 0))
   
   #do an ANOVA
   stat_test <- aov(f~r, plot_data)
   stat_text <- capture.output(summary(stat_test))
   stat_text <- stat_text[1:3]
   stat_name <- 'ANOVA'

 }
 
 #wild type vs mutant on x axis
 if (feature_type=='discrete' & resp_type=='continuous') {
   p <- ggplot(plot_data, aes(f,r) )
   p <- p + theme_bw() + scale_fill_brewer(palette="Set1") + scale_colour_brewer(palette="Set1")
   p <- p + scale_y_reverse()
   p <- p + labs(x=feature_name, y=resp_name)
   p <- p + geom_boxplot(aes(colour=f), outlier.size=0) #want a boxplot as well as points
   p <- p  + geom_point(aes(fill=f), size=7, shape=21, position = position_jitter(w = 0.2, h = 0))
   
   #do an ANOVA
   stat_test <- aov(r~f, plot_data)
   stat_text <- capture.output(summary(stat_test))
   stat_text <- stat_text[1:3]
   stat_name <- 'ANOVA'
   
 }
 
 #both axes discrete
 if (feature_type=='discrete' & resp_type=='discrete') {
   p <- ggplot(plot_data, aes(f,r ) )
   p <- p + theme_bw() + scale_fill_brewer(palette="Set1") + scale_colour_brewer(palette="Set1")
   p <- p + labs(x=feature_name, y=resp_name)
   p <- p  + geom_point(aes(fill=paste(r,f)), size=7, shape=21, position = position_jitter(w = 0.2, h = 0.2))
   
   #do a Fishers Exact test
   stat_test <- fisher.test(plot_data$f, plot_data$r)
   stat_text <- capture.output(print(stat_test))
   stat_text <- stat_text[c(5,7,8,10,11)]
   stat_name <- 'Fishers Exact Test'
   
 }
 
 #common settings
 p <- p + theme(axis.text = element_text(size=rel(1.5)), axis.title = element_text(size=rel(2)))
 p <- p + theme(legend.position = "none")
 
 #add stats text to plot
 annot <- tableGrob(stat_text,
                    cols=stat_name,
                    show.rownames=FALSE,
                    core.just='left',
                    padding.v = unit(2, "mm"),
                    gpar.coretext = gpar(col = "black", cex = 0.8),
                    gpar.corefill = gpar(fill = "white", col = "white"),
                    gp=gpar(fontfamily='mono'))
 
 #generate combined plot
 combo_plot <- arrangeGrob(p , annot, nrow=2, heights=c(4,1))
  
 
 return(list(plot=p,
             test=stat_test,
             annot=annot,
             data=plot_data,
             combo_plot=combo_plot))
 
  
}





