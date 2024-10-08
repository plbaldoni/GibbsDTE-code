methodsNames <- function(large.simulation = FALSE){
  if(!large.simulation){
    method <- c(c('edger-raw','edger-scaled','edger-new-raw','edger-new-scaled','sleuth-lrt','sleuth-wt','swish'),paste0('edger-scaled-n',seq(10,90,10)),paste0('edger-new-scaled-n',seq(10,90,10)))
    labels <- c('edgeR.legacy-RC','edgeR.legacy-SC','edgeR-RC','edgeR-SC','sleuth-LRT','sleuth-Wald','Swish',paste0('edgeR.legacy-SC.n',seq(10,90,10)),paste0('edgeR-SC.n',seq(10,90,10)))
    color <- c(c('black','orange','black','blue',"purple","green3","gray"),rep('black',9),rep('black',9)) # To match TranscriptDE paper
  } else{
    method <- c(c('edger-raw','edger-scaled','edger-new-raw','edger-new-scaled','sleuth-lrt','sleuth-wt','swish'),paste0('edger-scaled-n',c(2,3)),paste0('edger-new-scaled-n',c(2,3)))
    labels <- c('edgeR.legacy-RC','edgeR.legacy-SC','edgeR-RC','edgeR-SC','sleuth-LRT','sleuth-Wald','Swish',paste0('edgeR.legacy-SC.n',c(2,3)),paste0('edgeR-SC.n',c(2,3)))
    color <- c(c('black','orange','black','blue',"purple","green3","gray"),rep('black',2),rep('black',2)) # To match TranscriptDE paper
  }

  names(method) <- names(color) <- labels
  return(list(labels = labels,method = method,color = color))
}

#' @importFrom data.table data.table fread
loadResults <- function(path,genome,len,fc,read,tx.per.gene,scenario,libs.per.group,simulation,quantifier,large.simulation){

  meth <- methodsNames(large.simulation)
  path.time <- file.path(path,'time.tsv')
  path.method <- file.path(path,paste0(meth$method,'.tsv.gz'))
  names(path.method) <- meth$labels

  dt.scenario <- data.table('Genome' = genome,
                            'Length' = len,
                            'FC' = fc,
                            'Reads' = read,
                            'TxPerGene' = tx.per.gene,
                            'Scenario' = scenario,
                            'LibsPerGroup' = libs.per.group,
                            'Simulation' = simulation,
                            'Quantifier' = quantifier)

  dt.results <- lapply(seq_along(path.method),function(i){
    dt <- data.table('Method' = names(path.method)[i], fread(path.method[i],header = TRUE))
    setnames(x = dt,
             old = c('pvalue','pval','qvalue','qval','feature'),
             new = c('PValue','PValue','FDR','FDR','TranscriptID'),skip_absent = TRUE)

    dt <- dt[,c('Method','TranscriptID','PValue','FDR')]
    return(dt)
  })

  dt.results <- do.call(rbind,dt.results)
  dt.results <- cbind(dt.scenario,dt.results)

  dt.time <- fread(file = path.time,header = TRUE)
  dt.time <- cbind(dt.scenario,dt.time[,c('method','elapsed')])
  setnames(x = dt.time,old = c('method','elapsed'),new = c('Method','Time'))
  dt.time$Method <- names(meth$method)[match(dt.time$Method,meth$method)]

  out <- list('results' = dt.results,'time' = dt.time)

  return(out)
}

loadQuantResults <- function(path,genome,len,fc,read,tx.per.gene,scenario,libs.per.group,simulation,quantifier){

  path.time <- file.path(path,'time.tsv')

  dt.scenario <- data.table('Genome' = genome,
                            'Length' = len,
                            'FC' = fc,
                            'Reads' = read,
                            'TxPerGene' = tx.per.gene,
                            'Scenario' = scenario,
                            'LibsPerGroup' = libs.per.group,
                            'Simulation' = simulation,
                            'Quantifier' = quantifier)

  dt.time <- fread(file = path.time,header = TRUE)
  dt.time <- cbind(dt.scenario,dt.time[,c('method','elapsed')])
  setnames(x = dt.time,old = c('method','elapsed'),new = c('Sample','Time'))

  libsgroup <- as.numeric(gsub("libsPerGroup","",libs.per.group))
  n.libs <- rep(libsgroup,2)

  # Generating sample labels
  group <- rep(LETTERS[seq_len(length(n.libs))],times = n.libs)
  rep <- unlist(lapply(n.libs,seq_len))
  group.name <- paste0(paste(paste0('group', group), paste0('rep',rep), sep = '_'),"_R1")

  if (scenario == 'balanced') {
    lib.sizes <- rep(50e6,sum(n.libs))
  } else{
    lib.sizes <- rep(rep(c(25e6,100e6),length.out = libsgroup),2)
  }
  names(lib.sizes) <- group.name
  dt.time$LibSize <- lib.sizes[match(dt.time$Sample,names(lib.sizes))]

  # Getting time spend during bootstrapping/Gibbs resampling
  if (grepl('salmon',quantifier)){
    log.path <- file.path(dirname(path.time),dt.time$Sample,'logs/salmon_quant.log')
    log.time <- lapply(log.path,function(x){
      df.x <- readLines(x)
      if(grepl('gibbs',quantifier)){
        i <- c(grep("Starting Gibbs Sampler",df.x),grep("Finished Gibbs Sampler",df.x))
      } else{
        i <- c(grep("Optimizing over",df.x) + 1L,grep("Finished Bootstrapping",df.x))
      }
      time <- strsplit2(df.x[i],"\\[|\\]")[,2] |> as.POSIXlt(x = _)
      out <- as.numeric(difftime(time[2],time[1],units = 'sec'))
      return(out)
    })
    dt.time$ResamplingTime <- unlist(log.time)
  } else{
    # kallisto does not log bootstrapping time
    dt.time$ResamplingTime <- NA
  }

  out <- list('time' = dt.time)
  return(out)
}

loadMetadata <- function(path,genome,len,fc,read,tx.per.gene,scenario,libs.per.group,simulation){
  path.counts <- file.path(path,'counts.tsv.gz')

  dt.scenario <- data.table('Genome' = genome,
                            'Length' = len,
                            'FC' = fc,
                            'Reads' = read,
                            'TxPerGene' = tx.per.gene,
                            'Scenario' = scenario,
                            'LibsPerGroup' = libs.per.group,
                            'Simulation' = simulation)

  dt.metadata <- fread(path.counts,header = TRUE)

  dt.transcript <- dt.metadata[,c('TranscriptID')]

  dt.metadata <- cbind(dt.scenario,dt.metadata[!dt.metadata$status == 0, c('TranscriptID', 'status')])

  # If fold-change = 1, status should be 0
  if (fc == 1) dt.metadata$status <- 0

  out <- list('simulation' = dt.metadata,'transcript' = dt.transcript)

  return(out)
}

aggregateScenario <- function(path,genome,len,fc,read,tx.per.gene,scenario,libs.per.group,quantifier,nsim,large.simulation){

  subpath <- paste0('simulation-',seq_len(nsim))

  ls.quant <- lapply(seq_len(nsim),function(x){
    res.path <- file.path(path,subpath[x],paste0('quant-',quantifier))
    message("Reading simulations from:")
    message(res.path)
    loadQuantResults(res.path,genome,len,fc,read,tx.per.gene,scenario,libs.per.group,x,quantifier)
  })

  ls.results <- lapply(seq_len(nsim),function(x){
    res.path <- file.path(path,subpath[x],paste0('dte-',quantifier))
    message("Reading simulations from:")
    message(res.path)
    loadResults(res.path,genome,len,fc,read,tx.per.gene,scenario,libs.per.group,x,quantifier,large.simulation)
  })

  ls.metadata <- lapply(seq_len(nsim),function(x){
    meta.path <- file.path(path,subpath[x],'meta')
    message("Reading simulations from:")
    message(meta.path)
    loadMetadata(meta.path,genome,len,fc,read,tx.per.gene,scenario,libs.per.group,x)
  })

  dt.results <- do.call(rbind,lapply(ls.results,function(x){x[['results']]}))
  dt.time <- as.data.table(do.call(rbind,lapply(ls.results,function(x){x[['time']]})))
  dt.simulation <- do.call(rbind,lapply(ls.metadata,function(x){x[['simulation']]}))
  dt.transcript <- ls.metadata[[1]]$transcript
  dt.quanttime <- as.data.table(do.call(rbind,lapply(ls.quant,function(x){x[['time']]})))

  out <- list('results' = dt.results,'time' = dt.time,'simulation' = dt.simulation,'features' = dt.transcript,'quanttime' = dt.quanttime)
  return(out)
}

computeMetrics <- function(x,simulation,features,fdr,alpha){

  TranscriptID.DE <- x$TranscriptID[x$FDR < fdr]
  n <- length(x$TranscriptID)
  n.lt.alpha <- sum(x$PValue < alpha)

  if (length(TranscriptID.DE) > 0) {
    call.DE <- data.table(TranscriptID = TranscriptID.DE,call = 1)

    truth.DE <- simulation[Genome == x$Genome & Length == x$Length & FC == x$FC &
                             Reads == x$Reads & TxPerGene == x$TxPerGene &
                             Scenario == x$Scenario & LibsPerGroup == x$LibsPerGroup &
                             Simulation == x$Simulation, c('TranscriptID','status')]

    tb.DE <- merge(features,truth.DE,by = 'TranscriptID',all.x = TRUE)
    tb.DE <- merge(tb.DE,call.DE,by = 'TranscriptID',all.x = TRUE)

    tb.DE[is.na(status),status := 0]
    tb.DE[, status := abs(status)]
    tb.DE[is.na(call),call := 0]

    tb.DE$status <- factor(tb.DE$status,levels = c(0,1))
    tb.DE$call <- factor(tb.DE$call,levels = c(0,1))

    tb.results <- tb.DE[,table(status,call)]

    total.de <- sum(tb.results["1",])

    out <- list('N' = n,
                'N.ALPHA' = n.lt.alpha,
                'TP' = tb.results["1","1"],
                'FP' = tb.results["0","1"],
                'FPR' = tb.results["0","1"]/sum(tb.results["0",]),
                'FDR' = tb.results["0","1"]/sum(tb.results[,"1"]),
                'TPR' = ifelse(total.de == 0,NA,tb.results["1","1"]/total.de))
  } else{
    out <- list('N' = n,'N.ALPHA' = n.lt.alpha,'TP' = 0,'FP' = 0,'FPR' = 0,'FDR' = 0,'TPR' = 0)
  }
  return(lapply(out,as.double))
}

computeFDRCurve <- function(x,simulation,features,fdr,seq.n){

  truth.DE <- simulation[Genome == x$Genome & Length == x$Length & FC == x$FC &
                           Reads == x$Reads & TxPerGene == x$TxPerGene &
                           Scenario == x$Scenario & LibsPerGroup == x$LibsPerGroup &
                           Simulation == x$Simulation, c('TranscriptID','status')]

  tb.DE <- merge(features,truth.DE,by = 'TranscriptID',all.x = TRUE)
  tb.DE[is.na(status),status := 0]
  tb.DE[, status := abs(status)]

  feature.DE <- data.table(TranscriptID = x$TranscriptID,FDR = x$FDR,call = 1)

  tb.DE <- merge(tb.DE,feature.DE,by = 'TranscriptID',all.x = TRUE)
  tb.DE[is.na(call),call := 0]
  tb.DE$status <- factor(tb.DE$status,levels = c(0,1))
  tb.DE$call <- factor(tb.DE$call,levels = c(0,1))
  tb.DE <- tb.DE[order(FDR),]

  out <- lapply(seq.n,function(w){
    tb.results <- tb.DE[seq(1,w),][,table(status,call)]
    return(tb.results["0","1"])
  })

  names(out) <- paste0('n.',seq.n)

  return(out)
}

summarizePValue <- function(x,byvar,step = 0.05){
  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  cut.sq <- seq(step,1 - step,by = step)
  cut.match <- c(0,cut.sq) + step/2
  names(cut.match) <- paste0('(',roundPretty(c(0,cut.sq),2),'-',roundPretty(c(cut.sq,1),2),']')
  n.groups <- length(cut.match)

  x.sub <- copy(x)
  x.sub[,PValue := cut(x = PValue,breaks = c(-Inf,cut.sq,Inf),labels = names(cut.match))]

  x.sub.method <- x.sub[,list(N = .N),by = byvar]

  table <- x.sub[,list(N.cat = .N),by = c(byvar,'PValue')]
  table <- merge(table,x.sub.method,by = byvar,all.x = TRUE)
  table$PValue <- factor(table$PValue,levels = names(cut.match))

  table <- table[,list(Density.Avg = mean(n.groups*N.cat/N)),by = c(sub.byvar,'PValue')]

  table$PValue.Midpoint <- cut.match[match(table$PValue,names(cut.match))]

  return(table)
}

#' @importFrom data.table melt
summarizeFDRCurve <- function(x,byvar){

  cnames <- colnames(x)

  sub.cnames <- cnames[grepl('n\\.',cnames)]
  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  x.mean <- x[,lapply(.SD,mean),by = sub.byvar,.SDcols = sub.cnames]

  table <- melt(data = x.mean,id.vars = sub.byvar,variable.name = 'N',value.name = 'FDR')

  table[,N := as.numeric(gsub('n\\.','',N))]

  return(table)
}

summarizeMetrics <- function(x,byvar){

  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  table <- x[,.(P.SIG = mean(N.ALPHA/N),TP = mean(TP),FP = mean(FP)),sub.byvar]

  return(table)
}

summarizeTime <- function(x,byvar){

  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  table <- x[,.(Time = mean(Time/60)),sub.byvar]

  return(table)
}

summarizeQQ <- function(x,byvar,step = 0.001){

  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  x.quant <- x[, list(q.sample = quantile(PValue,probs = seq(0,1,length.out = .N)),
                      q.theory = seq(0,1,length.out = .N)),by = byvar]

  quant.sq <- seq(step,1 - step,by = step)
  quant.match <- c(0,quant.sq) + step/2
  names(quant.match) <- paste0('(',c(0,quant.sq),'-',c(quant.sq,1),']')

  x.quant[,Q.Theory.Cat := cut(x = q.theory,breaks = c(-Inf,quant.sq,Inf),labels = names(quant.match))]

  table <- x.quant[,.(Q.Sample.Avg = mean(q.sample)),by = c(sub.byvar,'Q.Theory.Cat')]

  table$Q.Theory.Midpoint <- quant.match[match(table$Q.Theory.Cat,names(quant.match))]

  return(table)
}

summarizeOverdispersion <- function(path, genome, len, fc, read, tx.per.gene,
                                    scenario, libs.per.group, quantifier, nsim){

  subpath <- paste0('simulation-',seq_len(nsim))
  catchFunction <- get(ifelse(grepl("salmon",quantifier),'catchSalmon2','catchKallisto2'))

  ls.results <- lapply(seq_len(nsim),function(x){
    res.path <- file.path(path,subpath[x],paste0('quant-',quantifier))
    meta.path <- file.path(path,subpath[x],'meta/counts.tsv.gz')

    meta <- fread(meta.path)

    catch <- catchFunction(list.dirs(res.path,recursive = FALSE))
    rownames(catch$annotation) <- strsplit2(rownames(catch$annotation),"\\|")[,1]

    keep <- rownames(catch$annotation) %in% meta$TranscriptID[!is.na(meta$status)]

    out <- catch$annotation$Overdispersion[keep]
    return(out)
  })

  res <- log10(unlist(ls.results))

  out <- data.table('Genome' = genome,
                    'Length' = len,
                    'FC' = fc,
                    'Reads' = read,
                    'TxPerGene' = tx.per.gene,
                    'Scenario' = scenario,
                    'LibsPerGroup' = libs.per.group,
                    'Quantifier' = quantifier,
                    'Mean' = mean(res),
                    'SD' = sd(res),
                    '2.5Pct' = quantile(res,0.025),
                    '25Pct' = quantile(res,0.25),
                    '50Pct' = quantile(res,0.5),
                    '75Pct' = quantile(res,0.75),
                    '97.5Pct' = quantile(res,0.975))

  return(out)
}

computeROCCurve <- function(x,simulation,features,seq.fdr){
  truth.DE <- simulation[Genome == x$Genome & Length == x$Length & FC == x$FC &
                           Reads == x$Reads &
                           Scenario == x$Scenario & LibsPerGroup == x$LibsPerGroup &
                           Simulation == x$Simulation, c('TranscriptID','status')]

  tb.DE <- merge(features,truth.DE,by = 'TranscriptID',all.x = TRUE)
  tb.DE[is.na(status),status := 0]
  tb.DE[, status := abs(status)]

  feature.DE <- data.table(TranscriptID = x$TranscriptID,FDR = x$FDR)

  tb.DE <- merge(tb.DE,feature.DE,by = 'TranscriptID',all.x = TRUE)
  tb.DE[,call := 0]
  tb.DE$status <- factor(tb.DE$status,levels = c(0,1))
  tb.DE$call <- factor(tb.DE$call,levels = c(0,1))
  tb.DE[is.na(FDR),FDR := 1]

  out <- lapply(seq.fdr,function(w){
    tb.results <- copy(tb.DE)
    tb.results[FDR < w/100, call := "1"]
    tb.results <- tb.results[,table(status,call)]
    return(c('TPR' = tb.results[2,2]/sum(tb.results[2,]),'FDR' = tb.results[1,2]/sum(tb.results[,2])))
  })

  out.tpr <- lapply(out,function(x){x['TPR']})
  out.fdr <- lapply(out,function(x){x['FDR']})

  names(out.tpr) <- paste0('tpr.',seq.fdr)
  names(out.fdr) <- paste0('fdr.',seq.fdr)

  return(c(out.tpr,out.fdr))
}

summarizeROCCurve <- function(x,byvar){

  cnames <- colnames(x)

  sub.cnames1 <- cnames[grepl('tpr\\.',cnames)]
  sub.cnames2 <- cnames[grepl('fdr\\.',cnames)]

  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  x.mean <- x[,lapply(.SD,mean),by = sub.byvar,.SDcols = c(sub.cnames1,sub.cnames2)]

  table1 <- melt(data = x.mean,id.vars = sub.byvar,measure.vars = sub.cnames1,variable.name = 'nFDR',value.name = 'oTPR')
  table2 <- melt(data = x.mean,id.vars = sub.byvar,measure.vars = sub.cnames2,variable.name = 'nFDR',value.name = 'oFDR')

  table1$nFDR <- as.numeric(gsub("tpr.","",table1$nFDR))/100
  table2$nFDR <- as.numeric(gsub("fdr.","",table2$nFDR))/100

  table <- merge(table1,table2,by = c(sub.byvar,'nFDR'),all.x = TRUE)

  return(table)
}

#' @importFrom data.table fwrite
summarizeQuantification <- function(path,dest,genome,fc,read,len,large.simulation,
                                    tx.per.gene,scenario,libs.per.group,quantifier,
                                    nsim = 20, fdr = 0.05, seq.n = seq(100,3000,100),seq.fdr = c(1,5,10,15,20),alpha = fdr){

  byvar <- c('Genome','Length','FC','Reads','TxPerGene','Scenario','LibsPerGroup','Quantifier','Method','Simulation')

  res <- aggregateScenario(path, genome, len, fc, read, tx.per.gene, scenario, libs.per.group, quantifier, nsim,large.simulation)

  table.metrics <- res$results[,computeMetrics(c(.BY,.SD),simulation = res$simulation,features = res$features,fdr = fdr,alpha = alpha),by = byvar]
  table.metrics.1pct <- res$results[,computeMetrics(c(.BY,.SD),simulation = res$simulation,features = res$features,fdr = 0.01,alpha = 0.01),by = byvar]
  table.fdr <- res$results[,computeFDRCurve(c(.BY,.SD),simulation = res$simulation,features = res$features,fdr = fdr,seq.n = seq.n),by = byvar]
  table.overdispersion <- summarizeOverdispersion(path, genome, len, fc, read, tx.per.gene, scenario, libs.per.group, quantifier, nsim)
  table.roc <- res$results[,computeROCCurve(c(.BY,.SD),simulation = res$simulation,features = res$features,seq.fdr = seq.fdr),by = byvar]

  out <- list('fdr' = summarizeFDRCurve(table.fdr,byvar),
              'metrics' = summarizeMetrics(table.metrics,byvar),
              'metrics.1pct' = summarizeMetrics(table.metrics.1pct,byvar),
              'time' = summarizeTime(res$time,byvar),
              'quantile' = summarizeQQ(res$results,byvar),
              'pvalue' = summarizePValue(res$results,byvar),
              'overdispersion' = table.overdispersion,
              'quanttime' = res$quanttime,
              'roc' = summarizeROCCurve(table.roc,byvar))

  # plotFDRCurve(out$fdr,max.n = max(seq.n))
  # plotPowerBars(out$metrics,fdr,max.n)
  # plotType1Error(out$metrics,alpha)
  # plotTime(out$time)
  # plotQQPlot(out$quantile)
  # plotPValues(out$pvalue)

  dir.create(dest,showWarnings = FALSE,recursive = TRUE)
  lapply(names(out),function(x){fwrite(x = out[[x]], file = file.path(dest,paste0(x,'.tsv.gz')),quote = FALSE,sep = '\t')})

  return(invisible())
}

summarizeScenario <- function(x,table,path,dest,large.simulation){
  dt <- as.character(table[x,])
  names(dt) <- colnames(table)
  table.names = c('fdr','metrics','time','quantile','pvalue','overdispersion','quanttime','roc')

  in.path <- file.path(path,do.call(file.path,as.list(dt)))
  out.path <- file.path(dest,do.call(file.path,as.list(dt)))

  # Check is simulation directory exists
  if (!dir.exists(in.path)) return(invisible())
  message("Reading simulations from:")
  message(in.path)

  # Summarizing results with Salmon (bootstrap)
  if (!all(file.exists(file.path(out.path,'dte-salmon',paste0(table.names,'.tsv.gz'))))){
    # Verbose
    summarizeQuantification(path = in.path,quantifier = 'salmon',
                            large.simulation = large.simulation,
                            dest = file.path(out.path,'dte-salmon'),
                            genome = dt['genome'],fc = dt['fc'],read = dt['read'],
                            tx.per.gene = dt['tx.per.gene'],scenario = dt['scenario'],
                            libs.per.group = dt['libs.per.group'],len = dt['len'])
  }

  # Summarizing results with Salmon (Gibbs)
  if (!all(file.exists(file.path(out.path,'dte-salmon-gibbs',paste0(table.names,'.tsv.gz'))))){
    # Verbose
    summarizeQuantification(path = in.path,quantifier = 'salmon-gibbs',
                            large.simulation = large.simulation,
                            dest = file.path(out.path,'dte-salmon-gibbs'),
                            genome = dt['genome'],fc = dt['fc'],read = dt['read'],
                            tx.per.gene = dt['tx.per.gene'],scenario = dt['scenario'],
                            libs.per.group = dt['libs.per.group'],len = dt['len'])
  }

  # Summarizing results with kallisto
  if (!all(file.exists(file.path(out.path,'dte-kallisto',paste0(table.names,'.tsv.gz'))))){
    # Verbose
    summarizeQuantification(path = in.path,quantifier = 'kallisto',
                            large.simulation = large.simulation,
                            dest = file.path(out.path,'dte-kallisto'),
                            genome = dt['genome'],fc = dt['fc'],read = dt['read'],
                            tx.per.gene = dt['tx.per.gene'],scenario = dt['scenario'],
                            libs.per.group = dt['libs.per.group'],len = dt['len'])
  }
}

summarizeSimulation <- function(path,
                                dest,
                                genome = c('mm39'),
                                len = c(50,75,100,125,150),
                                fc = c(1,2),
                                read = c('single-end','paired-end'),
                                tx.per.gene = c(2,3,4,5,9999),
                                scenario = c('balanced','unbalanced'),
                                libs.per.group = c(3,5,10),
                                large.simulation = FALSE,...){

  path <- normalizePath(path)
  dir.create(dest,showWarnings = FALSE,recursive = TRUE)
  dest <- normalizePath(dest)

  dt.scenario <- expand.grid('genome' = genome,
                             'len' = paste0('readlen-',len),
                             'fc' = paste0('fc',fc),
                             'read' = read,
                             'tx.per.gene' = paste0(tx.per.gene,'TxPerGene'),
                             'scenario' = scenario,
                             'libs.per.group' = paste0(libs.per.group,'libsPerGroup'),
                             stringsAsFactors = FALSE)

  bplapply(seq_len(nrow(dt.scenario)),summarizeScenario,table = dt.scenario,dest = dest,path = path,large.simulation = large.simulation,...)

  return(invisible())
}

