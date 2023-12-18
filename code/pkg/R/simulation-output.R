roundPretty <- function(x,digits = 1){
  formatC(round(x,digits),digits = digits,format = 'f')
}

#' @importFrom kableExtra kbl add_header_above kable_styling pack_rows cell_spec
#' @importFrom kableExtra landscape collapse_rows
tabulateMetrics <- function(x,cap,
                            seq.len = seq(50,150,25),
                            lib.group = c(3,5,10),
                            lib.size = c('50M','25/100M'),
                            color = TRUE,
                            font_size = NULL,
                            color.fdr = 0.05,
                            format = 'latex',...){

  dt <- copy(x)

  methods <- methodsNames()$labels

  dt$Length %<>% factor(levels = paste0(seq.len,'bp'))
  dt$LibsPerGroup %<>% mapvalues(from = paste0("#Lib/Group = ",lib.group),to = lib.group)
  dt$Scenario %<>% mapvalues(from = c('balanced','unbalanced'),to = lib.size) %<>% factor(levels = lib.size)

  dt[,A.Power := TP/3000]
  dt[,B.FDR := ifelse((FP+TP) == 0,NA,FP/(FP+TP))]

  dt.dcast <- dcast(dt, Reads + LibsPerGroup + Scenario + Length ~ Method,value.var = c('A.Power','B.FDR'))

  dt.dcast <- dt.dcast[order(Reads,LibsPerGroup,Scenario,Length),]
  dt.dcast <- dt.dcast[,c('Reads','LibsPerGroup','Scenario','Length',paste0('A.Power_',methods),paste0('B.FDR_',methods)),with = FALSE]

  dt.dcast[,c(paste0('A.Power_',methods),paste0('B.FDR_',methods)) := lapply(.SD,roundPretty,digits = 3),.SDcols = c(paste0('A.Power_',methods),paste0('B.FDR_',methods))]

  dt.dcast.color <- copy(dt.dcast)

  if(color == TRUE){
    mat.power <- as.matrix(dt.dcast[,paste0('A.Power_',methods),with = FALSE])
    class(mat.power) <- 'numeric'
    mat.fdr <- as.matrix(dt.dcast[,paste0('B.FDR_',methods),with = FALSE])
    class(mat.fdr) <- 'numeric'

    col.power <- t(sapply(1:nrow(mat.power),FUN = function(i){
      ifelse(mat.power[i,] == max(mat.power[i,which(mat.fdr[i,] < color.fdr)]) &
               mat.fdr[i,] < color.fdr,'blue','black')
    }))
    col.power[is.na(col.power)] <- 'black'
    col.fdr <- t(apply(mat.fdr,1,function(x){ifelse(x > color.fdr,'red','black')}))
    col.fdr[is.na(col.fdr)] <- 'black'

    for(imethod in methods){
      dt.dcast.color[[paste0('A.Power_',imethod)]] <-
        cell_spec(x = dt.dcast.color[[paste0('A.Power_',imethod)]],color = col.power[,paste0('A.Power_',imethod)],format = format)
      dt.dcast.color[[paste0('A.Power_',imethod)]] <- gsub("NA","-",dt.dcast.color[[paste0('A.Power_',imethod)]])

      dt.dcast.color[[paste0('B.FDR_',imethod)]] <-
        cell_spec(x = dt.dcast.color[[paste0('B.FDR_',imethod)]],color = col.fdr[,paste0('B.FDR_',imethod)],format = format)
      dt.dcast.color[[paste0('B.FDR_',imethod)]] <- gsub("NA","-",dt.dcast.color[[paste0('B.FDR_',imethod)]])
    }
  }

  kb <- kbl(dt.dcast.color,
            escape = FALSE,
            format = format,
            booktabs = TRUE,
            align = c('l',rep('r',13)),
            caption = cap,
            col.names = c('Read','Samples/Group','Library Size','Read Length',rep(methods,2)),...) %>%
    add_header_above(c(" " = 4, "Power" = length(methods), "False Discovery Rate" = length(methods))) %>%
    { if(format == 'latex'){
      dotDotDot <- match.call(expand.dots = FALSE)
      if('longtable' %in% names(dotDotDot$...)){
        getLongTable <- get('longtable')
      } else{
        getLongTable <- FALSE
      }
      if(isTRUE(getLongTable)){
        opts <- c("repeat_header")
      } else{
        opts <- c("scale_down")
      }
      kable_styling(kable_input = .,latex_options = opts,font_size = font_size) %>%
        collapse_rows(kable_input = .,columns = 1, latex_hline = "major", valign = "top")
    } else{
      kable_styling(kable_input = .,bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
        collapse_rows(kable_input = .,columns = 1, valign = "top") %>%
        landscape()
    }
    }

  return(kb)
}

plotQQPlot <- function(x,base_size = 8){
  meth <- methodsNames()

  plot <- ggplot(x,
                 aes(x = Q.Theory.Midpoint,y = Q.Sample.Avg,color = Method,group = Method)) +
    facet_grid(rows = vars(LibsPerGroup),cols = vars(TxPerGene)) +
    # geom_abline(intercept = 0,slope = 1,colour = 'black',linetype = 'dashed') +
    geom_line(inherit.aes = FALSE,aes(x=0,y = 0,color = Method,group = Method)) +
    geom_point(pch = '.',size = 2) +
    theme_bw(base_size = base_size,base_family = 'sans') +
    scale_color_manual(values = meth$color) +
    scale_x_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1)) +
    theme(strip.text = element_text(colour = 'black',size = base_size),
          legend.background = element_rect(fill = alpha('white', 0)),
          legend.text = element_text(size = base_size),
          legend.title = element_blank(),
          # legend.key.size = unit(2,"line"),
          legend.position = 'top',
          panel.grid = element_blank(),
          axis.text = element_text(colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size)) +
    labs(y = 'Sample Quantiles',x = 'Theoretical Quantiles')

  return(plot)
}


#' @importFrom ggplot2 facet_wrap geom_histogram geom_hline scale_x_discrete rel
plotPValues <- function(x,base_size = 8){
  plot <- ggplot(data = x,aes(x = PValue,y = Density.Avg)) +
    facet_grid(rows = vars(LibsPerGroup),cols = vars(Method)) +
    geom_col(col = 'black',fill = 'grey',position = position_dodge(0.9),width = 0.8) +
    geom_hline(yintercept = 1,col = 'red',linetype = 'dashed') +
    theme_bw(base_size = base_size,base_family = 'sans') +
    theme(strip.text = element_text(colour = 'black',size = base_size),
          legend.background = element_rect(fill = alpha('white', 0)),
          legend.text = element_text(size = base_size),
          legend.title = element_blank(),
          legend.key.size = unit(0.75,"line"),
          legend.position = 'top',
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90),
          axis.text = element_text(colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size)) +
    labs(x = 'P-values',y = 'Density')
  return(plot)
}

#' @importFrom ggplot2 geom_bar position_dodge element_text
plotTime <- function(x,base_size = 8){
  plot <- ggplot(data = x,aes(x = Method,y = Time)) +
    facet_grid(rows = vars(LibsPerGroup),cols = vars(TxPerGene)) +
    geom_bar(stat = 'identity',position = position_dodge()) +
    theme_bw(base_size = base_size,base_family = 'sans') +
    scale_y_continuous(limits = c(0,max(5,max(x$Time)))) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = 'top',
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 90,colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size),
          strip.text = element_text(colour = 'black',size = base_size)) +
    labs(y = 'Time (min)',x = NULL)
  return(plot)
}


#' @importFrom thematic okabe_ito
#' @importFrom ggplot2 ggplot geom_line theme_bw scale_color_manual
#' @importFrom ggplot2 scale_x_continuous theme element_blank labs aes
#' @importFrom ggplot2 scale_y_continuous geom_abline facet_grid vars
plotFDRCurve <- function(x,max.n,base_size = 8){

  meth <- methodsNames()

  plot <- ggplot(x,aes(x = N,y = FDR,color = Method,group = Method)) +
    facet_grid(rows = vars(LibsPerGroup),cols = vars(TxPerGene),scales = 'free_y') +
    geom_line(size = 0.75) +
    theme_bw(base_size = base_size,base_family = 'sans') +
    scale_color_manual(values = meth$color) +
    # scale_y_continuous(limits = c(0,1300)) +
    theme(strip.text = element_text(colour = 'black',size = base_size),
          legend.background = element_rect(fill = alpha('white', 0)),
          legend.text = element_text(size = base_size),
          legend.title = element_blank(),
          legend.key.size = unit(0.75,"line"),
          legend.position = 'top',
          panel.grid = element_blank(),
          axis.text = element_text(colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size)) +
    labs(y = 'False discoveries',x = 'Transcripts chosen')

  return(plot)
}

#' @importFrom ggplot2 geom_col geom_text scale_fill_brewer
plotPowerBars <- function(x,fdr,max.n,base_size = 8){

  sub.byvar <- colnames(x)[-which(colnames(x) %in% c('P.SIG','TP','FP'))]

  gap <- 0.05*max(x$TP + x$FP)

  x.melt <- melt(x,id.vars = sub.byvar,
                 measure.vars = c('TP','FP'),
                 variable.name = 'Type',
                 value.name = 'Value')
  x.melt$Type <- factor(x.melt$Type,levels = c('FP','TP'),labels = c('False Positive','True Positive'))

  plot <- ggplot(x.melt,aes(x = Method,y = Value,fill = Type)) +
    facet_grid(rows = vars(LibsPerGroup),cols = vars(TxPerGene)) +
    geom_col() +
    geom_text(aes(x = Method,y = (TP + FP) + gap,
                  label = roundPretty(ifelse((FP+TP) == 0,NA,100*FP/(FP+TP)),1)),
              inherit.aes = FALSE,data = x,vjust = 0,size = base_size/.pt) +
    theme_bw(base_size = base_size,base_family = 'sans') +
    scale_fill_brewer(palette = 'Set1') +
    labs(x = NULL,y = paste('# DE Transcripts at FDR <',roundPretty(fdr,2))) +
    scale_y_continuous(limits = c(0,max.n)) +
    theme(strip.text = element_text(colour = 'black',size = base_size),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = 'top',
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 90,colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size),
          legend.background = element_rect(fill = alpha('white', 0)),
          legend.text = element_text(size = base_size),
          legend.key.size = unit(0.75,"line"))

  return(plot)
}

plotType1Error <- function(x,alpha,base_size = 8){

  sub.byvar <- colnames(x)[-which(colnames(x) %in% c('P.SIG','TP','FP'))]

  x.melt <- melt(x,id.vars = sub.byvar,
                 measure.vars = c('P.SIG'),variable.name = 'Type',value.name = 'Value')

  plot <- ggplot(x.melt,aes(x = Method,y = Value)) +
    facet_grid(rows = vars(LibsPerGroup),cols = vars(TxPerGene)) +
    geom_col() +
    theme_bw(base_size = base_size,base_family = 'sans') +
    geom_hline(yintercept = alpha,color = 'red',linetype = 'dashed') +
    labs(x = NULL,y = paste('Proportion of transcripts with p-value <',roundPretty(alpha,2))) +
    theme(strip.text = element_text(colour = 'black',size = base_size),
          legend.background = element_rect(fill = alpha('white', 0)),
          legend.text = element_text(size = base_size),
          legend.title = element_blank(),
          legend.key.size = unit(0.75,"line"),
          legend.position = 'top',
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 90),
          axis.text = element_text(colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size))

  return(plot)
}
