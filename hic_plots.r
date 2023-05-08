
##### shaman plots #####
plot.whole.chrom <- function(chrom.cov, chrom.ab, chrom.scores, all.tracks.conts, chrom.index, plot.path, coord.range=NULL, gene.intervs=NULL, exp.per.cluster=NULL, rotate=F, mark.coord=NULL, smooth.ab=F, smooth.cov=F, cls.col=c(1, 2), point_size=1, sel.attr.plots=NULL, width=NULL, height.vec=NULL, show.coord.text=F, coord.on.bottom=NULL) {

  chrom.end = get.chrom.sizes()[chrom.index]
  if (is.null(coord.range)) {
    whole.chrom = T
    coord.range = c(1, chrom.end)
    breaks_seq = seq(5e6, chrom.end, 5e6)
  } else {
    whole.chrom = F
    breaks_seq = seq(1, chrom.end, 2e5)
    breaks_seq = breaks_seq[between(breaks_seq, coord.range[1], coord.range[2])]
  }
  breaks_text = breaks_seq - 1
  breaks_text[as.numeric(breaks_text) %% 1e6 != 0] = ''
  #if (!show.coord.text & length(height.vec) > 2) {
  #  breaks_text[T] = ''
  #}
  shaman_breaks_text = breaks_text
  lines_breaks_text = breaks_text
  if (!show.coord.text) {
    shaman_breaks_text[T] = ''
  }

  chr.name = paste0('chr', chrom.index)
  chrom.cov = as.data.frame(chrom.cov[[1]])
  chrom.ab = as.data.frame(chrom.ab[[1]])
  chrom.scores = as.data.frame(chrom.scores[[1]])

  chrom.cov.bin.size = min(diff(sort(unique(sapply(strsplit(rownames(chrom.cov), '_'), function(x) as.numeric(x[2]))))))
  chrom.ab.bin.size = min(diff(sort(unique(sapply(strsplit(rownames(chrom.ab), '_'), function(x) as.numeric(x[2]))))))

  chr.name_ = paste0(chr.name, '_')
  chrom.cov = chrom.cov[grepl(chr.name_, rownames(chrom.cov)),]
  chrom.ab = chrom.ab[grepl(chr.name_, rownames(chrom.ab)),]
  chrom.scores = chrom.scores[grepl(chr.name_, rownames(chrom.scores)),]

  chrom.cov = chrom.cov[between(sapply(strsplit(rownames(chrom.cov), '_'), function(x) as.numeric(x[2])), coord.range[1], coord.range[2]),]
  chrom.ab = chrom.ab[between(sapply(strsplit(rownames(chrom.ab), '_'), function(x) as.numeric(x[2])), coord.range[1], coord.range[2]),]
  chrom.scores = chrom.scores[between(sapply(strsplit(rownames(chrom.scores), '_'), function(x) as.numeric(x[2])), coord.range[1], coord.range[2]),]

  all.cov.bins = get.all.bins(chr.name, coord.range[1], coord.range[2], bin.size=chrom.cov.bin.size)[[1]]
  all.ab.bins = get.all.bins(chr.name, coord.range[1], coord.range[2], bin.size=chrom.ab.bin.size)[[1]]
  chrom.cov[setdiff(all.cov.bins, rownames(chrom.cov)),] = NA
  chrom.ab[setdiff(all.ab.bins, rownames(chrom.ab)),] = NA

  if (smooth.ab) {
    coord = as.numeric(sapply(strsplit(rownames(chrom.ab), '_'), function(x) x[[2]]))
    chrom.ab = chrom.ab[order(coord),]
    for (j in 1:ncol(chrom.ab)) {
      chrom.ab[,j] = rollmean(chrom.ab[,j], 5, na.rm=T, fill=NA)
      #chrom.ab[,j] = rollmean(chrom.ab[,j], 20, na.rm=T, fill=NA)
    }
  }
  if (smooth.cov) {
    coord = as.numeric(sapply(strsplit(rownames(chrom.cov), '_'), function(x) x[[2]]))
    chrom.cov = chrom.cov[order(coord),]
    for (j in 1:ncol(chrom.cov)) {
      chrom.cov[,j] = rollmean(chrom.cov[,j], 10, na.rm=T, fill=NA)
    }
  }

  chrom.conts1 = all.tracks.conts[[1]][all.tracks.conts[[1]]$chrom1 == chr.name, ]
  chrom.conts2 = all.tracks.conts[[2]][all.tracks.conts[[2]]$chrom1 == chr.name, ]
  chrom.conts = rbind(chrom.conts1[chrom.conts1$start1 < chrom.conts1$start2,], 
                      chrom.conts2[chrom.conts2$start1 > chrom.conts2$start2,])

  shaman.gplot1 = my.shaman.plot(chrom.conts1, list(start=coord.range[1], end=coord.range[2]), rotate=T, breaks_seq=breaks_seq, mark.coord=mark.coord, point_size=point_size, mark.on.top=T, mark.width=1.5, breaks_text=shaman_breaks_text)
  shaman.gplot2 = my.shaman.plot(chrom.conts2, list(start=coord.range[1], end=coord.range[2]), rotate=T, breaks_seq=breaks_seq, mark.coord=mark.coord, point_size=point_size, mark.on.top=T, mark.width=1.5, breaks_text=shaman_breaks_text)

  
  attr.names = c('cov', 'ab', 'scores')
  display.attr.names = c('Coverage', 'A score', 'Scores')
  attr.plots = list()
  for (i in seq_along(attr.names)) {
    attr.name = attr.names[i]
    attr.df = get(paste0('chrom.', attr.name))
    nlines = ncol(attr.df)
    colnames(attr.df) = paste0('attr', 1:nlines)
    attr.df$coord = as.numeric(sapply(strsplit(rownames(attr.df), '_'), function(x) x[[2]]))
    attr.plot = ggplot(data=attr.df)
    for (j in 1:nlines) {
      # workaround to avoid environment problems. Don't know how to solve it elegantly.
      if (j == 1) {
       attr.plot = attr.plot + geom_line(aes(x=coord, y=attr1, color='col1'), size=3)
      } else if (j == 2) {
       attr.plot = attr.plot + geom_line(aes(x=coord, y=attr2, color='col2'), size=3)
       #attr.plot = attr.plot + geom_smooth(aes(x=coord, y=attr2, color='col2'))
      } 
    }
    if (!is.null(coord.on.bottom)) {
      #attr.plot = attr.plot + geom_vline(xintercept = coord.on.bottom, size=1)
      attr.plot = attr.plot + geom_vline(xintercept = coord.on.bottom, size=1, colour=rgb(0, 0, 0, alpha=0.3))

    }
    attr.plot = attr.plot + scale_color_manual(values=c(col1=cls.col[1], col2=cls.col[2])) + 
                xlim(coord.range[1], coord.range[2]) + ylab(NULL) + theme(legend.position='none') +
                #xlim(min(all.bins.in.range1[[2]]), max(all.bins.in.range1[[2]])) + ylab(display.attr.names[i]) + theme(legend.position='none') +
                #xlim(min(attr.df$coord), max(attr.df$coord)) + ylab(display.attr.names[i]) + theme(legend.position='none') +
                #expand_limits(x = c(1, chrom.end), y = 0) #+ scale_x_continuous(expand = c(0, 0))
                theme(axis.title.x=element_blank(), panel.grid.minor.x = element_line(colour="grey", size=0.5), panel.grid.major.x = element_line(colour="grey", size=0.5)) + scale_x_continuous(expand = c(0, 0), minor_breaks=breaks_seq, breaks=breaks_seq, labels=lines_breaks_text, limits=coord.range) + expand_limits(x = coord.range, y = 0) 

    attr.plots[[i]] = attr.plot
  }

  if (!is.null(gene.intervs)) {
    gene.intervs = gene.intervs[gene.intervs$chrom == chr.name & between(gene.intervs$tss, coord.range[1], coord.range[2]), ]
    for (i in 1:ncol(exp.per.cluster)) {
      #gene.intervs[, paste0('exp', i)] = exp.per.cluster[match.gene.names(gene.intervs$geneSymbol, rownames(exp.per.cluster)), i]
      gene.intervs[, paste0('exp', i)] = exp.per.cluster[gene.intervs$geneSymbol, i]
    }
    
    exp.plot = ggplot(data=gene.intervs) + geom_point(aes(x=tss, y=exp1, color='col1')) + geom_point(aes(x=tss, y=exp2, color='col2')) + scale_color_manual(values=c(col1=cls.col[1], col2=cls.col[2])) + 
                xlim(coord.range[1], coord.range[2]) + ylab(NULL) + theme(legend.position='none') +
                theme(panel.grid.minor.x = element_line(colour="grey", size=0.5), panel.grid.major.x = element_line(colour="grey", size=0.5)) + expand_limits(x = coord.range, y = 0) + scale_x_continuous(expand = c(0, 0), minor_breaks=breaks_seq, breaks=breaks_seq, labels=lines_breaks_text)

    attr.plots = c(list(exp.plot), attr.plots)
  }

  if (!is.null(sel.attr.plots)) {
    attr.plots = attr.plots[sel.attr.plots]
  }

  #all.plots = do.call(align_plots, c(list(shaman.gplots), attr.plots, list(align='v')))
  all.plots = do.call(align_plots, c(list(shaman.gplot1, shaman.gplot2), attr.plots, list(align='v')))
  if (is.null(height.vec)) {
    if (whole.chrom) {
      width=40
      height.vec = c(20, 20, 2, 2, 2)
    } else {
      width=6
      if (!is.null(gene.intervs)) {
        height.vec = c(3, 3, 1.5, 1.5, 1.5, 1.5)
      } else {
        height.vec = c(3, 3, 1.5, 1.5, 1.5)
      }
    }
  } else {
    width = 6
  }
  all.plots = do.call(grid.arrange, c(all.plots, list(ncol=1, nrow=length(height.vec), heights=height.vec, widths=6)))
  if (!is.null(plot.path)) {
    for (cur.plot.path in plot.path) {
      dev = ifelse(whole.chrom, 'jpeg', 'png')
      ggsave(cur.plot.path, all.plots, device=dev, height=sum(height.vec), width=width, limitsize=F)
    }
  }
}

plot.shaman.gene <- function(track.conts, interv, fig.path, window.size=1e6, point_size=2.5, rotate=F) {
  chr.name = interv$chrom
  coord.range = c(interv$tss - window.size / 2, interv$tss + window.size / 2)
  chrom.index = as.numeric(substr(chr.name, 4, 6))
  chrom.end = get.chrom.sizes()[chrom.index]
  breaks_seq = seq(1, chrom.end, 2e5)
  breaks_seq = breaks_seq[between(breaks_seq, coord.range[1], coord.range[2])]
  chrom.conts = track.conts[track.conts$chrom1 == chr.name, ]
  chrom.conts = chrom.conts[chrom.conts$start1 < chrom.conts$start2,]

  shaman.gplots = my.shaman.plot(chrom.conts, list(start=coord.range[1], end=coord.range[2]), rotate=rotate, 
                                 breaks_seq=breaks_seq, mark.coord=interv$tss, text.interv=2e5, point_size=point_size, mark.on.top=T, mark.width=1.5)
  ggsave(fig.path, plot=shaman.gplots, device='png', width=12, height=6)
}

my.shaman.plot <- function (points_score, interval_range = NA, rotate = TRUE, point_size = 1,
    add_axis = TRUE, breaks_seq=NULL, mark.coord=NULL, text.interv=1e6, mark.on.top=F, mark.width=0.5, breaks_text=NULL) {
    
    if (!all(c("start1", "start2", "score") %in% colnames(points_score))) {
        stop("points_score data frame must contain the following columns: start1, start2, score")
    }
    if (is.na(interval_range)[1]) {
        interval_range = data.frame(start = min(points_score$start1),
            end = max(points_score$start1))
    }
    col.scores = shaman_score_pal()
    ylim_max = (interval_range$end - interval_range$start) / 2
    points_score$x_val = (points_score$start1 + points_score$start2) / 2
    points_score$y_val = (points_score$start2 - points_score$start1) / 2
    points_score = points_score[between(points_score$x_val, interval_range$start, interval_range$end) &
                                between(points_score$y_val, 0, ylim_max),]
    map_gplot <- ggplot2::ggplot(points_score[order(points_score$score),
        ], ggplot2::aes(x = x_val, y = y_val, color = factor(floor(score)))) + ggplot2::coord_cartesian(xlim = c(interval_range$start,
        interval_range$end),ylim = c(0, ylim_max)) + 
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank())
    if (!is.null(mark.coord)) {
      map_gplot <- map_gplot + ggplot2::geom_abline(slope=1, intercept=-mark.coord, color='black', size=0.5) +
                              ggplot2::geom_abline(slope=-1, intercept=mark.coord, color='black', size=0.5)
    }
    if (is.null(breaks_seq)) {
      breaks_seq = seq(5e6, interval_range$end, 5e6)
    }
    if (is.null(breaks_text)) {
      breaks_text = breaks_seq - 1
      breaks_text[as.numeric(breaks_text) %% text.interv != 0] = ''
    }
    #breaks_text[T] = ''
    map_gplot <- map_gplot + ggplot2::geom_point(size = point_size) + xlim(interval_range$start, interval_range$end) + 
        #theme(panel.grid.minor = element_line(colour="grey", size=1), panel.grid.major = element_line(colour="grey", size=1)) +
        ggplot2::scale_colour_manual(values = col.scores[101 +
            sort(unique(floor(points_score$score)))]) + ggplot2::scale_x_continuous(expand = c(0,
        0), minor_breaks=breaks_seq, breaks=breaks_seq, labels=breaks_text) + ggplot2::scale_y_continuous(expand = c(0, 0), minor_breaks=breaks_seq, breaks=breaks_seq, labels=breaks_text) +
        ggplot2::theme_bw() + ggplot2::theme(legend.position = "none",
        panel.border = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(), axis.line.x = ggplot2::element_line(size = 0.2)) + expand_limits(x = c(interval_range$start, interval_range$end))
    map_gplot = map_gplot + theme(panel.grid.minor = element_line(colour="grey", size=0.4), panel.grid.major = element_line(colour="grey", size=0.4))
    if (add_axis == FALSE) {
        map_gplot <- map_gplot + ggplot2::theme(axis.ticks = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_blank(), axis.line.y = ggplot2::element_blank(),
            axis.ticks.length = ggplot2::unit(0, "null"))
    }

    if (rotate & !is.null(mark.coord) & mark.on.top) {
      map_gplot <- map_gplot + ggplot2::geom_abline(slope=1, intercept=-mark.coord, color='black', size=mark.width) +
                              ggplot2::geom_abline(slope=-1, intercept=mark.coord, color='black', size=mark.width)
    }
    return(map_gplot)
}
