# ---- plots ----
library(ggplot2)

#' plot gene name annotations
#' @param dat a matrix of gene names with 'start' and 'end' base-pair position
#' @param xrange range of x axis, base-pair position
plot_geneName = function(dat, xrange){
  ngene = 2:nrow(dat)
  line = 1
  dat$lines = NA
  dat$lines[1] = 1
  gene.end = dat[1, 'end']
  while(length(ngene) != 0){
    id = which(dat[ngene, 'start'] > gene.end + 0.02)[1]
    if(!is.na(id)){
      dat$lines[ngene[id]] = line
      gene.end = dat[ngene[id],'end']
      ngene = ngene[-id]
    }else{
      line = line + 1
      dat$lines[ngene[1]] = line
      gene.end = dat[ngene[1],'end']
      ngene = ngene[-1]
    }
  }
  
  dat$mean = rowMeans(dat[,c('start', 'end')])
  
  pl = ggplot(dat, aes(xmin = xrange[1], xmax = xrange[2], ymin = -lines-0.4)) + 
    geom_rect(aes(xmin = start, xmax = end, ymin = -lines-0.05, ymax = -lines+0.05), fill='blue') +
    geom_text(aes(x = mean, y=-lines-0.2, label=geneName), size=2) + xlab('base-pair position') + ylab('Gene') + 
    theme_bw() + theme(axis.text.x=element_blank(),
                       axis.ticks = element_blank(),
                       axis.text.y=element_blank(),
                       plot.title=element_text(size=11),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
  pl
}

discrete_gradient_pal <- function(colours, bins = 5) {
  ramp <- scales::colour_ramp(colours)
  
  function(x) {
    if (length(x) == 0) return(character())
    
    i <- floor(x * bins)
    i <- ifelse(i > bins-1, bins-1, i)
    ramp(i/(bins-1))
  }
}

scale_colour_discrete_gradient <- function(..., colours, bins = 5, na.value = "grey50", guide = "colourbar", aesthetics = "colour", colors)  {
  colours <- if (missing(colours)) 
    colors
  else colours
  continuous_scale(
    aesthetics,
    "discrete_gradient",
    discrete_gradient_pal(colours, bins),
    na.value = na.value,
    guide = guide,
    ...
  )
}

#' Locuszoom plot
#' @param z a vector of z scores with SNP names
#' @param pos base-pair positions
#' @param gene.pos.map a matrix of gene names with 'start' and 'end' base-pair position
#' @param z.ref.name the reference SNP
#' @param ld correlations between teh reference SNP and the rests
#' @param title title of the plot
#' @param title.size the size of the title
#' @param true the true value
#' @param y.height height of -log10(p) plot and height of the gene name annotation plot
#' @param y.lim range of y axis
locus.zoom = function(z, pos, gene.pos.map=NULL, z.ref.name=NULL, ld=NULL, 
                      title = NULL, title.size = 10, true = NULL, 
                      y.height=c(5,1.5), y.lim=NULL){
  tmp = data.frame(POS = pos, p = -(pnorm(-abs(z), log.p = T) + log(2))/log(10))
  pl_zoom = ggplot(tmp, aes(x = POS, y = p)) + geom_point(color = 'darkblue') + 
    ylab("-log10(p value)") + ggtitle(title) + 
    theme_bw() + theme(axis.title.x=element_blank(),
                       plot.title = element_text(size=title.size))
  if(!is.null(ld) && !is.null(z.ref.name)){
    tmp$ref = names(z) == z.ref.name
    tmp$r2 = ld^2
    pl_zoom = ggplot(tmp, aes(x = POS, y = p, shape = ref, size=ref, color=r2)) + geom_point() + 
      ylab("-log10(p value)") + ggtitle(title) + 
      scale_colour_discrete_gradient(
        colours = c("darkblue", "deepskyblue", "lightgreen", "orange", "red"),
        limits = c(0, 1.01),
        breaks = c(0,0.2,0.4,0.6,0.8,1),
        guide = guide_colourbar(nbin = 100, raster = FALSE, frame.colour = "black", ticks.colour = NA)
      ) + 
      scale_shape_manual(values=c(20, 18), guide=FALSE) + scale_size_manual(values=c(2,5), guide=FALSE) + 
      theme_bw() + theme(axis.title.x=element_blank(),
                         plot.title = element_text(size=title.size))
  }
  
  if(!is.null(y.lim)){
    pl_zoom = pl_zoom + ylim(y.lim[1], y.lim[2])
  }
  pl_zoom = pl_zoom + geom_hline(yintercept=-log10(5e-08), linetype='dashed', color = 'red')
  if(!is.null(true)){
    tmp.true = data.frame(POS = which(true!=0), p = tmp$p[which(true!=0)], 
                          ref = (names(z) == z.ref.name)[which(true!=0)],
                          label = paste0('SNP',1:length(which(true!=0))))
    pl_zoom = pl_zoom + geom_point(data=tmp.true, aes(x=POS, y=p), 
                                   color='red', show.legend = FALSE, shape=1, stroke = 1) + 
      geom_text(data=tmp.true, aes(x = POS-30, y=p+1, label=label), size=3, color='red')
  }
  if(!is.null(gene.pos.map)){
    pl_gene = plot_geneName(gene.pos.map, xrange = c(min(pos), max(pos)))
    g = egg::ggarrange(pl_zoom, pl_gene, nrow=2, heights = y.height, draw=FALSE)
  }else{
    g = pl_zoom
  }
  g
}

#' SuSiE plot with Locuszoom plot
#' @param z a vector of z scores with SNP names
#' @param model the fitted SuSiE model
#' @param pos base-pair positions
#' @param gene.pos.map a matrix of gene names with 'start' and 'end' base-pair position
#' @param z.ref.name the reference SNP
#' @param ld correlations between teh reference SNP and the rests
#' @param title title of the plot
#' @param title.size the size of the title
#' @param true the true value
#' @param plotz whether to plot locuszoom plot
#' @param y.lim range of y axis
#' @param y.susie the y axis of the SuSiE plot, 'PIP' or 'p', 'p' refers to -log10(p)
susie_plot_locuszoom = function(z, model, pos, gene.pos.map = NULL, z.ref.name, ld, 
                                title = NULL, title.size = 10, true = NULL, 
                                plotz = TRUE, y.lim=NULL, y.susie='PIP'){
  if(plotz){
    pl_zoom = locus.zoom(z, pos = pos, ld=ld, z.ref.name = z.ref.name, title = title, title.size = title.size, y.lim=y.lim)
  }
  pip = model$pip
  tmp = data.frame(POS = pos, PIP = pip, p = -(pnorm(-abs(z), log.p = T) + log(2))/log(10))
  if(y.susie == 'PIP'){
    if(plotz){
      pl_susie = ggplot(tmp, aes(x = POS, y = PIP)) + geom_point(show.legend = FALSE, size=3) + 
        theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
    }else{
      pl_susie = ggplot(tmp, aes(x = POS, y = PIP)) + geom_point(show.legend = FALSE, size=3) + 
        theme_bw() + theme(axis.title.x=element_blank())
    }
  }else if(y.susie == 'p'){
    if(plotz){
      pl_susie = ggplot(tmp, aes(x = POS, y = p)) + geom_point(show.legend = FALSE, size=3) + 
        ylab("-log10(p value)") + 
        theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
      pl_susie = pl_susie + geom_hline(yintercept=-log10(5e-08), linetype='dashed', color = 'red')
    }else{
      pl_susie = ggplot(tmp, aes(x = POS, y = p)) + geom_point(show.legend = FALSE, size=3) +
        ylab("-log10(p value)") + 
        theme_bw() + theme(axis.title.x=element_blank())
    }
    
  }
  
  if(!is.null(true)){
    tmp.true = data.frame(POS = pos[which(true!=0)], PIP = pip[which(true!=0)])
    pl_susie = pl_susie + geom_point(data=tmp.true, aes(x=POS, y=PIP), 
                                     color='red', size=3, show.legend = FALSE)
  }
  
  model.cs = model$sets$cs
  if(!is.null(model.cs)){
    colors = c('red', 'cyan', 'green', 'orange', 'dodgerblue', 'violet', 'gold',
               '#FF00FF', 'forestgreen', '#7A68A6')
    for(i in 1:length(model.cs)){
      tmp.cs = tmp[model.cs[[i]],]
      if(y.susie == 'PIP'){
        pl_susie = pl_susie + geom_point(data=tmp.cs, aes(x=POS, y=PIP), 
                                       color = colors[i], size=3, show.legend = FALSE, shape=1, stroke = 2)
      }else if(y.susie == 'p'){
        pl_susie = pl_susie + geom_point(data=tmp.cs, aes(x=POS, y=p), 
                                       color = colors[i], size=3, show.legend = FALSE, shape=1, stroke = 2)
      }
    }
  }
  
  if(!is.null(gene.pos.map)){
    pl_gene = plot_geneName(gene.pos.map, xrange = c(min(pos), max(pos)))
    if(plotz){
      g = egg::ggarrange(pl_zoom, pl_susie, pl_gene, nrow=3, heights = c(4,4,1.5), draw=FALSE)
    }else{
      g = egg::ggarrange(pl_susie, pl_gene, nrow=2, heights = c(4,2), draw=FALSE)
    }
    
  }else{
    if(plotz){
      g = egg::ggarrange(pl_zoom, pl_susie, nrow=2, heights = c(4,4), draw=FALSE)
    }else{
      g = pl_susie
    }
  }
  g
}
