get.de.ab.df <- function(de.ab.bins, chrom.ab) {
  de.ab.bins[[1]]$start = as.numeric(de.ab.bins[[1]]$start)
  de.ab.bins[[2]]$start = as.numeric(de.ab.bins[[2]]$start)
  de.ab.bins[[1]]$end = de.ab.bins[[1]]$start + 2e5 * as.numeric(de.ab.bins[[1]]$seq_len) - 1
  de.ab.bins[[2]]$end = de.ab.bins[[2]]$start + 2e5 * as.numeric(de.ab.bins[[2]]$seq_len) - 1
  de.ab.bins[[1]]$is_ecto = T
  de.ab.bins[[1]]$is_meso = F
  de.ab.bins[[2]]$is_ecto = F
  de.ab.bins[[2]]$is_meso = T

  ecto.a.score.data = do.call(rbind, lapply(1:nrow(de.ab.bins[[1]]), function(i) {
    cur.interv = de.ab.bins[[1]][i,,drop=F]
    bin.names = get.all.bins(cur.interv$chrom, cur.interv$start, cur.interv$end - 1)[[1]]
    all.diffs = chrom.ab[[1]][bin.names, 1] - chrom.ab[[1]][bin.names, 2]
    sel.bin = bin.names[which.max(all.diffs)]
    ret = c(chrom.ab[[1]][sel.bin,1], chrom.ab[[1]][sel.bin,2], max(all.diffs))
    names(ret) = c('ecto_a', 'meso_a', 'diff')
    return(ret)
  }))
  meso.a.score.data = do.call(rbind, lapply(1:nrow(de.ab.bins[[2]]), function(i) {
    cur.interv = de.ab.bins[[2]][i,,drop=F]
    bin.names = get.all.bins(cur.interv$chrom, cur.interv$start, cur.interv$end - 1)[[1]]
    all.diffs = chrom.ab[[1]][bin.names, 2] - chrom.ab[[1]][bin.names, 1]
    sel.bin = bin.names[which.max(all.diffs)]
    ret = c(chrom.ab[[1]][sel.bin,1], chrom.ab[[1]][sel.bin,2], -max(all.diffs))
    names(ret) = c('ecto_a', 'meso_a', 'diff')
    return(ret)
  }))
  ab.df = rbind(cbind(de.ab.bins[[1]], ecto.a.score.data),
                cbind(de.ab.bins[[2]], meso.a.score.data))
  return(ab.df)
}

ecto.meso.a.analysis <- function() {
  set.misha(SCHIC.MISHA.PATH)
  repl.data.dir = get.data.file.dir()
  load(MODEL.PATH)
  load(CLUSTERING.PATH)
  decay.metrics = get.common.decay.metrics()
  repl.erys = get.repl.erys(decay.metrics)
  erys = get.erys()
  repl.esc = get.repl.esc(esc_sch_decay_metrics)

  full.cluster.assignment = all.cluster.assignment
  full.cluster.assignment[repl.erys] = 3
  full.cluster.assignment[repl.esc] = 4
  rownames(ab.model.params$bin.clusters) = rownames(ab.model.params$total.cov)

  ab.form.path = file.path(repl.data.dir, 'ab_form_all_cells')
  ab.form = get.or.create(ab.form.path, function() total.ab.format(unique(bin_cell_ab$cell), bin_cell_ab))
  ab.form.orig.cls = ab.form
  ab.form$atab = ab.form$atab[names(full.cluster.assignment),]
  ab.form$btab = ab.form$btab[names(full.cluster.assignment),]

  unnorm.total.cov = bin_cell_cov
  deeply.covered.bins = rownames(unnorm.total.cov)[rowSums(unnorm.total.cov[, names(all.cluster.assignment)]) > 5e3]
  for (chrom in c('chrM', 'chrX', 'chrY')) {
    deeply.covered.bins = deeply.covered.bins[!grepl(chrom, deeply.covered.bins)]
  }
  stopifnot(all(rownames(ab.model.params$total.cov) %in% deeply.covered.bins))
  cls.col = MODEL.CLUSTER.COLORS[-2]

  chrom.cov = compute.chrom.cov(full.cluster.assignment, unnorm.total.cov[deeply.covered.bins,], decay.metrics)
  chrom.ab = compute.chrom.ab(full.cluster.assignment, ab.form, min.cov.discard=300, min.cov.per.clust=100, min.cov.cells=names(all.cluster.assignment))

  de.ab.bins = lapply(1:2, function(i) get.ab.up.bins(chrom.ab, i, 3 - i, 0.098))
  de.ab.bins.df.list = lapply(de.ab.bins, group.bins)
  print(lapply(de.ab.bins.df.list, nrow))

  print('writing table S2')
  de.ab.df = get.de.ab.df(de.ab.bins.df.list, chrom.ab)
  genes.per.bin = sapply(1:nrow(de.ab.df), function(i) {
    cur.chrom = de.ab.df[i, 'chrom']
    cur.start = de.ab.df[i, 'start']
    cur.end = de.ab.df[i, 'end']
    genes.in.bin = intersect(unique(filter(gene.md, chrom == cur.chrom & between(start, cur.start, cur.end))$geneSymbol), rownames(e9.exp@e_gc))
    genes.in.bin.str = paste(genes.in.bin, collapse=', ')
    return(genes.in.bin.str)
  })
  de.ab.df$genes = genes.per.bin
  de.ab.df$start = trimws(format(de.ab.df$start, scientific=F))
  de.ab.df$end = trimws(format(de.ab.df$end, scientific=F))
  write.csv(de.ab.df, file=file.path(repl.data.dir, 'ecto_meso_ab_diff'))

  gene.md = get.gene.metadata.jumps()
  exp.ret = get.or.create(file.path(repl.data.dir, 'atlases'), function() get.atlases(gene.md))
  atlas = exp.ret$atlas
  e9.exp = exp.ret$e9.exp

  ab.matching.ret = metacell.matching(log2(e9.exp@mc_fp), chrom.ab, 1:2, de.genes=NULL)
  exp.per.cluster = do.call(cbind, lapply(1:2, function(i) {
    cluster.exp = apply(log2(e9.exp@e_gc + 1e-5)[, ab.matching.ret$top.mc.indices[[i]]], 1, max)
    return(cluster.exp)
  }))
  write.csv(exp.per.cluster, file=MESO.ECTO.EXP)
}

get.enhancers <- function(mod.name='k4mono') {
  tname = list(k4mono='h3k4me1_rep1', k27='h3k27me3_sum')[[mod.name]]
  track.names = sprintf("ENCODE.bing_ren.e10_5.%s.%s", c('forebrain', 'heart', 'hindbrain', 'limb', 'midbrain'), tname)
  vtrack1 = gvtrack.create('vtrack1', track.names[1], 'global.percentile')
  vtrack2 = gvtrack.create('vtrack2', track.names[2], 'global.percentile')
  vtrack3 = gvtrack.create('vtrack3', track.names[3], 'global.percentile')
  vtrack4 = gvtrack.create('vtrack4', track.names[4], 'global.percentile')
  vtrack5 = gvtrack.create('vtrack5', track.names[5], 'global.percentile')
  
  intervs = gscreen('log2(1-vtrack1) < -7 | log2(1-vtrack2) < -7 | log2(1-vtrack3) < -7 | log2(1-vtrack4) < -7 | log2(1-vtrack5) < -7')
  intervs$mid = (intervs$start + intervs$end) / 2
  
  return(intervs)
}

get.enhancers.with.scores <- function(enhancer.intervs=NULL, mod.name='k4mono', use.max=T) {
  if (is.null(enhancer.intervs)) {
    enhancer.intervs = get.enhancers()
  }

  tname = list(k4mono='h3k4me1_rep1', k27='h3k27me3_sum')[[mod.name]]
  track.names = sprintf("ENCODE.bing_ren.e10_5.%s.%s", c('forebrain', 'heart', 'hindbrain', 'limb', 'midbrain'), tname)
  vtrack1 = gvtrack.create('vtrack1', track.names[1], 'global.percentile')
  vtrack2 = gvtrack.create('vtrack2', track.names[2], 'global.percentile')
  vtrack3 = gvtrack.create('vtrack3', track.names[3], 'global.percentile')
  vtrack4 = gvtrack.create('vtrack4', track.names[4], 'global.percentile')
  vtrack5 = gvtrack.create('vtrack5', track.names[5], 'global.percentile')

  all.enh.scores = NULL
  for (i in 1:5) {
    print(paste('calulating enhancer scores', i))
    cur.vtrack = paste0('vtrack', i)
    enh.scores = gextract(cur.vtrack, enhancer.intervs)
    if (is.null(all.enh.scores)) {
      all.enh.scores = enh.scores
    }
    all.enh.scores[,cur.vtrack] = -log2(1 - enh.scores[,cur.vtrack])
  }

  if (use.max) {
    score.per.interv = do.call(rbind, tapply(1:nrow(all.enh.scores), 
           all.enh.scores$intervalID, 
  	 function(x) apply(all.enh.scores[x, paste0('vtrack', 1:5)], 2, max)))
  } else {
    score.per.interv = do.call(rbind, tapply(1:nrow(all.enh.scores), 
           all.enh.scores$intervalID, 
  	 function(x) apply(all.enh.scores[x, paste0('vtrack', 1:5)], 2, mean)))
  }
  ret = cbind(enhancer.intervs, score.per.interv)
  return(ret)
}


plot.tbx.exp <- function(fig.dir=get.fig.dir()) {
  scdb_init(EXPRESSION.MC.DIR, force_reinit=T)
  mc = scdb_mc("e9_recolor2")
  lfp = log2(mc@mc_fp)
  png(file.path(fig.dir, 'tbx_exp.png'))
  plot(lfp['Tbx3',], lfp['Tbx5',], bg=mc@colors, pch=21, cex=3)
  grid(col='grey', lwd=0.5)
  dev.off()
}


plot.histone.tracks <- function(interv, k4.scores, k27.scores, all.tracks.conts, fig.path) {
  cols = MODEL.CLUSTER.COLORS[-2]

  interv.chrom = interv$chrom
  all.scores = list(k4.scores, k27.scores)
  all.scores = lapply(all.scores, function(hist.scores) {
    hist.scores$ecto = apply(hist.scores[,paste0('vtrack', c(1, 3, 5))], 1, max)
    hist.scores$meso = apply(hist.scores[,paste0('vtrack', c(2, 4))], 1, max)
    hist.scores.fil = filter(hist.scores, chrom == interv.chrom)
    hist.scores.fil = filter(hist.scores.fil, between(hist.scores.fil$mid, interv$start, interv$end))
    return(hist.scores.fil)
  }) 
  num.panels = 5
  png(fig.path, width=1600, height=num.panels * 400)
  par(mfrow=c(num.panels, 1), mar=c(1, 3, 1, 1))

  # plot 4c
  gene.4c.fil = get.4c.trace.mult.dist(all.tracks.conts, interv, use.max.dist=T)
  gene.4c.fil$coord = sapply(strsplit(rownames(gene.4c.fil), '_'), function(x) as.numeric(x[2]))
  gene.4c.fil = gene.4c.fil[order(gene.4c.fil$coord),]
  xlim = c(interv$start, interv$end)
  # Patch to mark the coordinate
  plot(gene.4c.fil$coord, gene.4c.fil[,1], type='l', col='white', lwd=4, xlim=xlim, ylim=c(-100, 100), ylab='', xlab='', xaxs='i', main='', cex.lab=4, yaxs='i', labels=F)
    if (!is.null(interv$tss1)) {
      abline(v=interv$tss1, lwd=15, col=rgb(0, 0, 0, alpha=0.3))
      abline(v=interv$tss2, lwd=15, col=rgb(0, 0, 0, alpha=0.3))
    }

  abline(v=seq(- 1e6, 500e6, 2e5), col='grey', lwd=4)
  abline(h=seq(-50, 50, 50), col='grey', lwd=4)
  for (i in 1:length(all.tracks.conts)) {
    lines(gene.4c.fil$coord, gene.4c.fil[,i], col=cols[i], lwd=20, xlim=xlim, ylim=c(-100, 100))
  }

  ylabs=c('k4mono', 'k27')
  num.tracks = 2
  for (i in 1:num.tracks) {
    hist.scores = all.scores[[i]]
    ecto.scores = hist.scores$ecto[order(hist.scores$mid)]
    meso.scores = hist.scores$meso[order(hist.scores$mid)]
    coords = hist.scores$mid[order(hist.scores$mid)]

    ecto.scores.smoothed = ecto.scores
    meso.scores.smoothed = meso.scores

    ylim.ecto = c(5, 15)
    plot(coords, ecto.scores.smoothed, ylim=ylim.ecto, xlim=xlim, type='h', col=cols[1], ylab=paste(ylabs[i], 'ecto'), cex.lab=3, labels=F, lwd=5, xaxs='i')
    abline(h=7)
    abline(h=10)
    abline(h=13)
    abline(v=seq(- 1e6, 500e6, 2e5), col='grey', lwd=4)
    if (!is.null(interv$tss1)) {
      abline(v=interv$tss1, lwd=15, col=rgb(0, 0, 0, alpha=0.3))
      abline(v=interv$tss2, lwd=15, col=rgb(0, 0, 0, alpha=0.3))
    }

    ylim.meso = c(5, 15)
    plot(coords, meso.scores.smoothed, ylim=ylim.meso, xlim=xlim, type='h', col=cols[2], ylab=paste(ylabs[i], 'meso'), cex.lab=3, labels=F, lwd=5, xaxs='i')
    abline(h=7)
    abline(h=10)
    abline(h=13)
    abline(v=seq(- 1e6, 500e6, 2e5), col='grey', lwd=4)
    if (!is.null(interv$tss1)) {
      abline(v=interv$tss1, lwd=15, col=rgb(0, 0, 0, alpha=0.3))
      abline(v=interv$tss2, lwd=15, col=rgb(0, 0, 0, alpha=0.3))
    }
  }
  dev.off()
}

epigenetic.analysis <- function() {
  cls.col=MODEL.CLUSTER.COLORS[-2]
  repl.data.dir = get.data.file.dir()
  fig.dir = file.path(get.fig.dir(), 'epigenetic_analysis')
  dir.create(fig.dir, showWarnings=F)

  exp.per.cluster.orig = read.csv(MESO.ECTO.EXP, stringsAsFactors=F)
  exp.per.cluster = exp.per.cluster.orig[!duplicated(exp.per.cluster.orig[,1]),]
  rownames(exp.per.cluster) = exp.per.cluster[,1]
  exp.per.cluster = exp.per.cluster[,-1]


  set.misha(SHAMAN.MISHA.PATH)
  track.names = paste0(EMB.CLS.ANALYSIS.TRACK.NAMES, '_equal_size')
  score.track.names = paste0(track.names[1:2], '_matshuff_scores')
  all.tracks.conts = lapply(score.track.names, function(track.name) gextract(track.name, gintervals.2d.all(), colnames='score'))

  set.misha(SCHIC.MISHA.PATH)
  gene.md = get.gene.metadata.jumps()
  gene.md$tss = gene.md$start
  gene.md = gene.md[!duplicated(gene.md$geneSymbol),]
  rownames(gene.md) = gene.md$geneSymbol

  enhancer.data.dir = file.path(repl.data.dir, 'enhancer_dir')
  dir.create(enhancer.data.dir, showWarnings=F)
  enhancer.path = file.path(enhancer.data.dir, 'enhancers')
  k27.path = file.path(enhancer.data.dir, 'k27')
  k4.int = get.or.create(enhancer.path, get.enhancers)
  k27.int = get.or.create(k27.path, function() get.enhancers('k27'))

  k4.scores = get.or.create(file.path(enhancer.data.dir, 'k4mono_scores'), function() get.enhancers.with.scores(k4.int, 'k4mono'))
  k27.scores = get.or.create(file.path(enhancer.data.dir, 'k27_scores'), function() get.enhancers.with.scores(k27.int, 'k27'))

  k4.shaman.scores = get.or.create(file.path(enhancer.data.dir, 'k4_shaman_scores_matshuff'), function() enhancer.meso.ecto.analysis2(k4.scores, exp.per.cluster, gene.md, all.tracks.conts, enh.thresh=9))
  k27.shaman.scores = get.or.create(file.path(enhancer.data.dir, 'k27_shaman_scores_matshuff'), function() enhancer.meso.ecto.analysis2(k27.scores, exp.per.cluster, gene.md, all.tracks.conts, enh.thresh=7, is.enh=F))

  print('printing k4 statistics')
  unfiltered.enh.scores = k4.shaman.scores$unfiltered.enh.scores
  print('number of ecto enhancers')
  print(sum(unfiltered.enh.scores$is_selected & unfiltered.enh.scores$is_ecto))
  print('number of meso enhancers')
  print(sum(unfiltered.enh.scores$is_selected & unfiltered.enh.scores$is_meso))

  ecto.up.genes = rownames(exp.per.cluster)[exp.per.cluster[,1] - exp.per.cluster[,2] > 1]
  meso.up.genes = rownames(exp.per.cluster)[exp.per.cluster[,1] - exp.per.cluster[,2] < -1]
  print('number of ecto genes')
  print(length(ecto.up.genes))
  print('number of meso genes')
  print(length(meso.up.genes))
  print('number of either genes')
  print(length(c(ecto.up.genes, meso.up.genes)))

  plot.enh.dist.distrib(k4.shaman.scores$unfiltered.enh.scores, k4.shaman.scores$gene.md, data.dir=enhancer.data.dir, fig.dir=fig.dir)

  plot.selected.enhancers(k4.shaman.scores, fig.dir, 'k4_selection.png', cols=cls.col)
  plot.selected.enhancers(k27.shaman.scores, fig.dir, 'k27_selection.png', cols=cls.col)
  print('plotting diff k4')
  plot.diff.enhancers(k4.shaman.scores, fig.dir, 'k4_shaman_scores.png', cols=cls.col)
  print('plotting diff k27')
  plot.diff.enhancers(k27.shaman.scores, fig.dir, 'k27_shaman_scores.png', cols=cls.col)

  loc.center = 120120677
  epigenetic.hot.locus = list(chrom='chr5', tss=loc.center, tss1=loc.center, tss2=120284671, start=loc.center - 5e5, end=loc.center + 5e5)
  plot.shaman.gene(all.tracks.conts[[1]], epigenetic.hot.locus, file.path(fig.dir, 'tbx3_zoom1.png'), rotate=T)
  plot.shaman.gene(all.tracks.conts[[2]], epigenetic.hot.locus, file.path(fig.dir, 'tbx3_zoom2.png'), rotate=T)

  shaman.diff.thresh = 15
  print('k4 statisitics')
  print('ecto 3-way pairings')
  is.ecto.support = k4.shaman.scores[[2]]$is_ecto & ((k4.shaman.scores[[2]]$V1 > 0 & k4.shaman.scores[[2]]$V1 > pmax(k4.shaman.scores[[2]]$V2, 0) + shaman.diff.thresh) |
                                                     (k4.shaman.scores[[2]]$V1 < 0 & k4.shaman.scores[[2]]$V1 > k4.shaman.scores[[2]]$V2 + shaman.diff.thresh))
  is.shaman.ecto.support =((k4.shaman.scores[[2]]$V1 > 0 & k4.shaman.scores[[2]]$V1 > pmax(k4.shaman.scores[[2]]$V2, 0) + shaman.diff.thresh) |
                                                     (k4.shaman.scores[[2]]$V1 < 0 & k4.shaman.scores[[2]]$V1 > k4.shaman.scores[[2]]$V2 + shaman.diff.thresh))
  print(mean(k4.shaman.scores[[2]][is.shaman.ecto.support, 'is_ecto']))
  print(sum(is.ecto.support))

  print('meso 3-way pairings')
  is.meso.support = k4.shaman.scores[[2]]$is_meso & ((k4.shaman.scores[[2]]$V2 > 0 & k4.shaman.scores[[2]]$V2 > pmax(k4.shaman.scores[[2]]$V1, 0) + shaman.diff.thresh) | 
                                             (k4.shaman.scores[[2]]$V2 < 0 & k4.shaman.scores[[2]]$V2 > k4.shaman.scores[[2]]$V1 + shaman.diff.thresh))
  is.shaman.meso.support = ((k4.shaman.scores[[2]]$V2 > 0 & k4.shaman.scores[[2]]$V2 > pmax(k4.shaman.scores[[2]]$V1, 0) + shaman.diff.thresh) | 
                                             (k4.shaman.scores[[2]]$V2 < 0 & k4.shaman.scores[[2]]$V2 > k4.shaman.scores[[2]]$V1 + shaman.diff.thresh))
  print(mean(k4.shaman.scores[[2]][is.shaman.meso.support, 'is_meso']))
  print(sum(is.meso.support))

  k4.table = k4.shaman.scores[[2]]
  is.three.way.support = is.ecto.support | is.meso.support
  k4.is.three.way.support = is.three.way.support
  three.way.support = k4.table[is.three.way.support, c('geneSymbol', 'chrom', 'tss', 'mid', 'is_ecto', 'is_meso', 'V1', 'V2')]
  colnames(three.way.support) = c('geneSymbol', 'chrom', 'tss', 'enhancer_coord', 'is_ecto', 'is_meso', 'ecto_shaman', 'meso_shaman')
  write.csv(three.way.support, file=file.path(enhancer.data.dir, 'k4_epi_hot'))


  k27.table = k27.shaman.scores[[2]]
  is.k27.ecto.support = k27.table$is_ecto & ((k27.table$V2 > 0 & k27.table$V2 > pmax(k27.table$V1, 0) + shaman.diff.thresh) |
                                              (k27.table$V2 < 0 & k27.table$V2 > k27.table$V1 + shaman.diff.thresh))
  is.k27.meso.support = k27.table$is_meso & ((k27.table$V1 > 0 & k27.table$V1 > pmax(k27.table$V2, 0) + shaman.diff.thresh) |
                                              (k27.table$V1 < 0 & k27.table$V1 > k27.table$V2 + shaman.diff.thresh))
  is.k27.three.way.support = is.k27.ecto.support | is.k27.meso.support

  k27.three.way.support = k27.table[is.k27.three.way.support, c('geneSymbol', 'chrom', 'tss', 'mid', 'is_ecto', 'is_meso', 'V1', 'V2')]
  colnames(k27.three.way.support) = c('geneSymbol', 'chrom', 'tss', 'enhancer_coord', 'is_k27_site_ecto', 'is_k27_site_meso', 'ecto_shaman', 'meso_shaman')
  write.csv(k27.three.way.support, file=file.path(enhancer.data.dir, 'k27_epi_hot'))


  set.misha(SCHIC.MISHA.PATH)
  k4.resolution = 1000
  all.intervs = data.frame(chrom=epigenetic.hot.locus$chrom, 
                           start=seq(loc.center - 2e6, loc.center + 2e6 - 1, by=k4.resolution), 
                           end=seq(loc.center - 2e6 + k4.resolution, loc.center + 2e6, by=k4.resolution))
  all.intervs$mid = (all.intervs$start + all.intervs$end) / 2
  k4.for.epi.locus = get.or.create(file.path(enhancer.data.dir, 'k4_for_epi_locus'), function() get.enhancers.with.scores(all.intervs, 'k4mono', use.max=F))
  k27.for.epi.locus = get.or.create(file.path(enhancer.data.dir, 'k27_for_epi_locus'), function() get.enhancers.with.scores(all.intervs, 'k27', use.max=F))
  plot.histone.tracks(epigenetic.hot.locus, k4.for.epi.locus, k27.for.epi.locus, all.tracks.conts, fig.path=file.path(fig.dir, 'tbx_histone_track.png'))
  plot.tbx.exp(fig.dir)

  # plot 4C for genes
  set.misha(SHAMAN.MISHA.PATH)
  genes.to.plot = c('Sox2', 'Twist1', 'Foxd2', 'Msx2', 'Foxb1', 'Rfx4', 'Gata6', 'Foxc1', 'Otx2', 'Pou3f2')
  epigenetic.loci = lapply(genes.to.plot, function(gene) {
    gene.tss = gene.md[gene, 'tss']
    enh.options = filter(k4.shaman.scores[[2]], geneSymbol == gene)
    enh.tss = enh.options$mid[which.max(abs(enh.options$V1 - enh.options$V2))]
    #return(list(chrom=gene.md[gene, 'chrom'], tss=gene.tss, tss1=gene.tss, tss2=NULL, start=gene.tss- 1e6, end=gene.tss + 1e6, gene_name=gene))
    return(list(chrom=gene.md[gene, 'chrom'], tss=gene.tss, tss1=gene.tss, tss2=enh.tss, start=gene.tss- 1e6, end=gene.tss + 1e6, gene_name=gene))
  })
  for (i in seq_along(epigenetic.loci)) {
    cur.locus = epigenetic.loci[[i]]
    print(paste(i, cur.locus$gene_name))
    plot.histone.tracks(cur.locus, k4.scores, k27.scores, all.tracks.conts, 
                fig.path=file.path(fig.dir, sprintf('%s_histone_track.png', cur.locus$gene_name)))
  }

  # p-value window analysis
  all.shaman.scores2 = k4.shaman.scores$all.shaman.scores.with.unselected
  compute.shaman.window.pvals(all.shaman.scores2, enhancer.data.dir, fig.dir=fig.dir)


  set.misha(SHAMAN.MISHA.PATH)
  diff.regions = get.or.create(file.path(enhancer.data.dir, 'diff_shaman_regions'), function() get.diff.scored.regions(score.track.names, threshold=40))
  is.real.diff = ((diff.regions$orig.scores1 * diff.regions$orig.scores2 >= 0) & abs(diff.regions$orig.scores1 - diff.regions$orig.scores2) > 40) | 
                 ((diff.regions$orig.scores1 * diff.regions$orig.scores2 <= 0) & pmax(abs(diff.regions$orig.scores1), abs(diff.regions$orig.scores2)) > 40) 
  diff.regions = diff.regions[is.real.diff,]
  k4.shaman.intervs = k4.shaman.scores[[2]][k4.is.three.way.support,]
  #k4.shaman.intervs = k4.shaman.scores[[2]]
  k4.shaman.intervs$chrom1 = k4.shaman.intervs$chrom
  k4.shaman.intervs$chrom2 = k4.shaman.intervs$chrom
  k4.shaman.intervs$start1 = k4.shaman.intervs$start
  k4.shaman.intervs$start2 = k4.shaman.intervs$mid
  k4.shaman.intervs$end1 = k4.shaman.intervs$start + 1
  k4.shaman.intervs$end2 = k4.shaman.intervs$mid + 1
  k4.shaman.intervs = k4.shaman.intervs[,c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2')]
  diff.regions.intervs = diff.regions[,c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2')]
  #neighbors = gintervals.neighbors(diff.regions.intervs, k4.shaman.intervs, maxdist1=1e5, maxdist2=1e5)
  neighbors = gintervals.neighbors(diff.regions.intervs, k4.shaman.intervs, maxdist1=2e5, maxdist2=2e5)
  print('printing statisitics for diff shaman explained by 3-way')
  print(c(nrow(neighbors), nrow(diff.regions), nrow(neighbors) / nrow(diff.regions)))

}


enhancer.meso.ecto.analysis2 <- function(enh.scores, exp.per.cluster, gene.md, all.tracks.conts, enh.thresh=9, is.enh=T) {
  enh.scores$ecto = apply(enh.scores[,paste0('vtrack', c(1, 3, 5))], 1, max)
  enh.scores$meso = apply(enh.scores[,paste0('vtrack', c(2, 4))], 1, max)

  enh.scores$is_meso = enh.scores$meso > enh.thresh & enh.scores$meso > enh.scores$ecto + 3
  enh.scores$is_ecto = enh.scores$ecto > enh.thresh & enh.scores$ecto > enh.scores$meso + 3
  enh.scores$is_selected = enh.scores$is_meso | enh.scores$is_ecto
  unfiltered.enh.scores = enh.scores

  ecto.up.genes = rownames(exp.per.cluster)[exp.per.cluster[,1] - exp.per.cluster[,2] > 1]
  meso.up.genes = rownames(exp.per.cluster)[exp.per.cluster[,1] - exp.per.cluster[,2] < -1]

  de.genes = c(ecto.up.genes, meso.up.genes)
  gene.md$is_ecto_exp = rownames(gene.md) %in% ecto.up.genes
  gene.md$is_meso_exp = rownames(gene.md) %in% meso.up.genes
  exp.genes = rownames(exp.per.cluster)
  enh.to.gene = get.enh.to.genes2(gene.md, ecto.up.genes, meso.up.genes, enh.scores, exp.genes, is.enh=is.enh)
  all.enh = do.call(rbind, enh.to.gene)

  all.shaman.scores.list = mclapply(1:19, function(i) {
    chr.name = paste0('chr', i)
    cur.chrom.enh = filter(all.enh, chrom == chr.name)
    if (nrow(cur.chrom.enh) == 0) {
      return(NULL)
    }
    shaman.scores = as.data.frame(sapply(1:length(all.tracks.conts), function(j) {
      rel.scores = filter(all.tracks.conts[[j]], chrom1 == chr.name)
      knn.ret = get.knnx(data.frame(rel.scores$start1, rel.scores$start2),
                         data.frame(cur.chrom.enh$mid, gene.md[cur.chrom.enh$gene, 'tss']), k=1)
      return(rel.scores[knn.ret$nn.index, 'score', drop=F])
    }))
    colnames(shaman.scores) = paste0('V', 1:length(all.tracks.conts))
    print(c('in bla', i, nrow(cur.chrom.enh)))
    return(cbind(cur.chrom.enh, shaman.scores))
  }, mc.cores=25)
  all.shaman.scores = do.call(rbind, all.shaman.scores.list)

  if (is.enh) {
    return(list(all.enh=all.enh, all.shaman.scores=filter(all.shaman.scores, orig_selection), unfiltered.enh.scores=unfiltered.enh.scores, 
                gene.md=gene.md, all.shaman.scores.with.unselected=all.shaman.scores))
  } else {
    return(list(all.enh=all.enh, all.shaman.scores=all.shaman.scores, unfiltered.enh.scores=unfiltered.enh.scores, gene.md=gene.md))
  }
}

get.enh.to.genes2 <- function(gene.md, ecto.up.genes, meso.up.genes, enh.scores, exp.genes=NULL, is.enh=T) {
  max.enh.dist = 5e5
  min.enh.dist = 5e4
  if (is.null(exp.genes)) {
    exp.genes = gene.md$geneSymbol
  }
  gene.to.enh = mclapply(1:nrow(enh.scores), function(i) {
    cur.scores = enh.scores[i,]
    genes = filter(gene.md, chrom == cur.scores$chrom & 
          between(tss, cur.scores$mid - max.enh.dist, cur.scores$mid + max.enh.dist) &
          !between(tss, cur.scores$mid - min.enh.dist, cur.scores$mid + min.enh.dist) &
	  geneSymbol %in% exp.genes)
    genes$dist = cur.scores$mid - genes$tss
    ecto.up.in.range = which(genes$geneSymbol %in% ecto.up.genes)
    meso.up.in.range = which(genes$geneSymbol %in% meso.up.genes)

    if (nrow(genes) > 8) {
      return(NULL)
    } else {
      ret = NULL
      up.genes = filter(genes, dist > 0)
      if (nrow(up.genes) > 0) {
        cand.gene = up.genes[which.min(up.genes$dist),]
	if (is.enh) {
          if ((cand.gene$geneSymbol %in% genes[ecto.up.in.range, 'geneSymbol'] & cur.scores$is_ecto) |
              (cand.gene$geneSymbol %in% genes[meso.up.in.range, 'geneSymbol'] & cur.scores$is_meso)) {
	    cand.gene$orig_selection = T
          } else {
	    cand.gene$orig_selection = F
	  }
          ret = rbind(ret, cand.gene)
	} else {
          if ((cand.gene$geneSymbol %in% genes[ecto.up.in.range, 'geneSymbol'] & cur.scores$is_meso) |
              (cand.gene$geneSymbol %in% genes[meso.up.in.range, 'geneSymbol'] & cur.scores$is_ecto)) {
            ret = rbind(ret, cand.gene)
          }
	}
      }

      down.genes = filter(genes, dist < 0)
      if (nrow(down.genes) > 0) {
        cand.gene = down.genes[which.max(down.genes$dist),]
	if (is.enh) {
          if ((cand.gene$geneSymbol %in% genes[ecto.up.in.range, 'geneSymbol'] & cur.scores$is_ecto) |
              (cand.gene$geneSymbol %in% genes[meso.up.in.range, 'geneSymbol'] & cur.scores$is_meso)) {
	    cand.gene$orig_selection = T
	  } else {
	    cand.gene$orig_selection = F
	  }
          ret = rbind(ret, cand.gene)
	} else {
          if ((cand.gene$geneSymbol %in% genes[ecto.up.in.range, 'geneSymbol'] & cur.scores$is_meso) |
              (cand.gene$geneSymbol %in% genes[meso.up.in.range, 'geneSymbol'] & cur.scores$is_ecto)) {
            ret = rbind(ret, cand.gene)
          }
	}
      }

      if (!is.null(ret)) {
        ret$mid = cur.scores$mid
	ret$is_ecto = cur.scores$is_ecto
	ret$is_meso = cur.scores$is_meso
      }
      return(ret)
    }
  }, mc.cores=25)
  return(gene.to.enh)
}

plot.selected.enhancers <- function(shaman.scores, fig.dir, fig.name, cols) {
  all.enh = shaman.scores[[3]]
  lim = c(0, 25)
  is.selected = all.enh$is_selected
  png(file.path(fig.dir, fig.name), height=400, width=400)
  plot(all.enh$ecto[!is.selected], all.enh$meso[!is.selected], col='grey', pch=19, cex=0.2, xlim=lim, ylim=lim)
  points(all.enh$ecto[is.selected], all.enh$meso[is.selected], col='black', pch=19, cex=0.5)
  dev.off()
}

plot.diff.enhancers <- function(shaman.scores, fig.dir, fig.suffix, cols) {
  is.ecto = shaman.scores[[2]][,'is_ecto']
  pos.shaman = filter(shaman.scores[[2]], V1 > 0 & V2 > 0)
  is.ecto.pos = pos.shaman[,'is_ecto']
  print('number of points in fig')
  print(length(shaman.scores[[2]][is.ecto,'V1']))
  png(file.path(fig.dir, paste0('ecto_', fig.suffix)), height=400, width=400)
  plot(shaman.scores[[2]][is.ecto,'V1'], shaman.scores[[2]][is.ecto,'V2'], col=cols[1], pch=19,
       main='', xlab='', ylab='', cex=1)
  abline(a=0, b=1, col='black', lwd=3)
  grid(col='grey', lwd=2)
  dev.off()
  print('number of points in fig')
  print(length(shaman.scores[[2]][!is.ecto,'V1']))
  png(file.path(fig.dir, paste0('meso_', fig.suffix)), height=400, width=400)
  plot(shaman.scores[[2]][!is.ecto,'V1'], shaman.scores[[2]][!is.ecto,'V2'], col=cols[2], pch=19,
       main='', xlab='', ylab='', cex=1)
  abline(a=0, b=1, col='black', lwd=3)
  grid(col='grey', lwd=2)
  dev.off()
  ecto.density = density(pos.shaman[is.ecto.pos, 'V1'] - pos.shaman[is.ecto.pos, 'V2'])
  meso.density = density(pos.shaman[!is.ecto.pos, 'V1'] - pos.shaman[!is.ecto.pos, 'V2'])
  max.y = max(c(ecto.density$y, meso.density$y))
  png(file.path(fig.dir, paste0('enh_shaman_diff_', fig.suffix)), height=400, width=400)
  plot(ecto.density, col=cols[1], ylim=c(0, max.y), xlab='', ylab='', main='ecto shaman - meso shaman density', lwd=8)
  lines(meso.density, col=cols[2], lwd=10)
  dev.off()
}

get.diff.scored.regions <- function(score.track.names, threshold=20, overlap.dist=5e4, max.dist=NULL) {
  
  num.chroms = 19
  chrom.sizes = get.chrom.sizes()
  all.selection.regions = do.call(rbind, mclapply(1:num.chroms, function(chrom.num) {
    selected.regions = NULL
    chrom.name = paste0('chr', chrom.num)
    chrom.length = chrom.sizes[chrom.num]
    # jumps of 2MB, filter of 1MB
    for (coord in seq(1, chrom.length, 2e6)) {
      if (coord %% 2e7 == 1) {print(c(chrom.num, coord))}
      interv = gintervals.force_range(data.frame(chrom1=chrom.name, start1=coord - 2e6, end1=coord + 2e6, 
                                                   chrom2=chrom.name, start2=coord - 2e6, end2=coord + 2e6))

      
      point.score1 = gextract(score.track.names[1], interv, colnames="score")
      point.score2 = gextract(score.track.names[2], interv, colnames="score")
      if (is.null(point.score1) | is.null(point.score2)) {
        next
      }
      # this also takes care of duplicates - we keep only contacts where start2 > start1
      point.score1$dist = point.score1$start2 - point.score1$start1
      score.filtered1 = point.score1[point.score1$dist < 1e6 & point.score1$dist > 1e4,]
      point.score2$dist = point.score2$start2 - point.score2$start1
      score.filtered2 = point.score2[point.score2$dist < 1e6 & point.score2$dist > 1e4,]
      if (nrow(score.filtered1) == 0 | nrow(score.filtered2) == 0) {
        next
      }

      score.diff = compute.shaman.diff(score.filtered1, score.filtered2)
      if (!is.null(max.dist)) {
        score.diff = score.diff[score.diff$nn.dist < max.dist, ]
      }
    
      should.continue = T

      while (should.continue) {
        #highest.score.index = which.max(score.diff$score)
        highest.score.index = which.max(abs(score.diff$score))
        highest.scoring.region = score.diff[highest.score.index,]
        highest.score = highest.scoring.region$score
        if (length(highest.score) != 0) {
          if (abs(highest.score) >= threshold) {
            score.diff = score.diff[abs(score.diff$start1 - highest.scoring.region$start1) > overlap.dist |
                                          abs(score.diff$start2 - highest.scoring.region$start2) > overlap.dist,]
            if (!any(abs(selected.regions$start1 - highest.scoring.region$start1) <= overlap.dist &
                                          abs(selected.regions$start2 - highest.scoring.region$start2) <= overlap.dist)) {
              selected.regions = rbind(selected.regions, highest.scoring.region)
            }
          } else {
            should.continue = F
          } 
        } else {
          should.continue = F
        }
      }
    }
    return(selected.regions)
  }, mc.cores=num.chroms)) 
  return(all.selection.regions)
}

plot.enh.dist.distrib <- function(unfiltered.enh.scores, gene.md, data.dir, fig.dir, cls.col=MODEL.CLUSTER.COLORS[-2]) {
  max.enh.dist = 5e5
  min.enh.dist = 5e4
  gene.md = filter(gene.md, is_ecto_exp | is_meso_exp)
  enh.scores = filter(unfiltered.enh.scores, is_meso | is_ecto)
  all.dists.options = get.or.create(file.path(data.dir, 'all_gene_enh_dists'), function() {
    all.dists.options = list()
    for (is.ecto.genes in c(T, F)) {
      cur.gene.md = filter(gene.md, is_ecto_exp == is.ecto.genes)
      for (is.ecto.enh in c(T, F)) {
        cur.dists = c()
        cur.enh = filter(enh.scores, is_ecto == is.ecto.enh)
        for (i in 1:nrow(cur.enh)) {
          chrom.genes = filter(cur.gene.md, chrom == cur.enh[i, 'chrom'])
          dists = chrom.genes$tss - cur.enh[i, 'mid']
          dists = dists[between(abs(dists), min.enh.dist, max.enh.dist)]
          if (any(dists > 0)) {
            cur.dists = c(cur.dists, min(dists[dists > 0]))
          }
          if (any(dists < 0)) {
            cur.dists = c(cur.dists, max(dists[dists < 0]))
          }
          all.dists.options[[sprintf('ecto_gene_%s_ecto_enh_%s', is.ecto.genes, is.ecto.enh)]] = cur.dists
        }
      }
    }
    return(all.dists.options)
  })

  all.dens = lapply(all.dists.options, function(x) density(log2(abs(x))))
  ylim=range(sapply(all.dens, "[", "y"), finite=T)

  cols = c(cls.col[1], 'black', 'saddlebrown', cls.col[2])
  dens.names = c('Ecto genes, ecto enhacners', 'Ecto genes, meso enhancers', 'Meso genes, ecto enhancers', 'Meso genes, meso enhancers')
  for (j in 1:2) {
    fig.name = c('enh_gene_dist_distrib.png', 'enh_gene_dist_distrib_with_legend.png')[j]
    png(file.path(fig.dir, fig.name), width=1000, height=1000)
    for (i in 1:length(all.dists.options)) {
      cur.name = dens.names[i]
      cur.dens = all.dens[[i]]
      if (i == 1) {
        plot(cur.dens, col=cols[i], type='l', lwd=8, main='', ylim=ylim)
      } else {
        lines(cur.dens, col=cols[i], lwd=8)
      }
    }
    grid(col='grey', lwd=0.5)
    if (j == 2) {
      legend('topleft', legend=dens.names, col=cols, lty=1, cex=2, lwd=8)
    }
    dev.off()
  }
}

compute.shaman.diff <- function(shaman.scores1, shaman.scores2) {
  stopifnot(all(shaman.scores1$chrom1 == shaman.scores1$chrom2) & all(shaman.scores2$chrom1 == shaman.scores2$chrom2) & 
            all(shaman.scores1$chrom1 == shaman.scores2$chrom1))
  knn.ret = get.knnx(data.frame(shaman.scores1$start1, shaman.scores1$start2), 
                     data.frame(shaman.scores2$start1, shaman.scores2$start2), k=1)
  orig.scores1 = shaman.scores1[knn.ret$nn.index, 'score']
  orig.scores2 = shaman.scores2$score
  nn.dist = knn.ret$nn.dist
  new.scores = (orig.scores2 - orig.scores1)

  shaman.scores2$score = new.scores
  shaman.scores2$orig.scores1 = orig.scores1
  shaman.scores2$orig.scores2 = orig.scores2
  shaman.scores2$nn.dist = nn.dist

  return(shaman.scores2)
}

compute.shaman.window.pvals <- function(all.shaman.scores, enhancer.data.dir, fig.dir) {
  if (is.null(all.conts.unscored)) {
    unscored.tnames = c('ecto_for_meso_comparison', 'meso_for_ecto_comparison')
    all.conts.unscored = lapply(unscored.tnames, function(track.name) gextract(track.name, gintervals.2d.all()))
  }
  all.shaman.scores$min_coord = pmin(all.shaman.scores$mid, all.shaman.scores$tss)
  all.shaman.scores$max_coord = pmax(all.shaman.scores$mid, all.shaman.scores$tss)
  window.size = 5e4
  intervs = data.frame(chrom1=all.shaman.scores$chrom, start1=all.shaman.scores$min_coord - window.size / 2, end1=all.shaman.scores$min_coord + window.size / 2,
                       chrom2=all.shaman.scores$chrom, start2=all.shaman.scores$max_coord - window.size / 2, end2=all.shaman.scores$max_coord + window.size / 2,
                       interv_id=sprintf('%s_%s_%s', all.shaman.scores$chrom, all.shaman.scores$min_coord, all.shaman.scores$max_coord))
  
  loci.scores = get.or.create(file.path(enhancer.data.dir, 'prev_windows_more_conts'), function() mclapply(paste0('chr', 1:19), function(cur.chrom) {
    chrom.intervs = filter(intervs, chrom1 == cur.chrom)
    track1.conts = filter(all.conts.unscored[[1]], chrom1 == cur.chrom & abs(start1 - start2) < 6e5 & abs(start1 - start2) > 2e4)
    track2.conts = filter(all.conts.unscored[[2]], chrom1 == cur.chrom & abs(start1 - start2) < 6e5 & abs(start1 - start2) > 2e4)
    cur.chrom.counts = do.call(rbind, lapply(1:nrow(chrom.intervs), function(j) {
      if (j %% 1000 == 1) print(paste(cur.chrom, j))
      count1 = sum(track1.conts$start1 >= chrom.intervs$start1[j] & track1.conts$start1 <= chrom.intervs$end1[j] &
                   track1.conts$start2 >= chrom.intervs$start2[j] & track1.conts$start2 <= chrom.intervs$end2[j])
      count2 = sum(track2.conts$start1 >= chrom.intervs$start1[j] & track2.conts$start1 <= chrom.intervs$end1[j] &
                   track2.conts$start2 >= chrom.intervs$start2[j] & track2.conts$start2 <= chrom.intervs$end2[j])
      return(c(count1, count2))
    }))
    return(cur.chrom.counts)
  }, mc.cores=20))

  all.shaman.scores$interv_id = sprintf('%s_%s_%s', all.shaman.scores$chrom, all.shaman.scores$min_coord, all.shaman.scores$max_coord)
  for (i in 1:19) {
    cur.chrom = paste0('chr', 1:19)[i]
    chrom.intervs = filter(intervs, chrom1 == cur.chrom)
    colnames(loci.scores[[i]]) = c('ecto_counts', 'meso_counts')
    loci.scores[[i]] = cbind(loci.scores[[i]], all.shaman.scores[match(chrom.intervs$interv_id, all.shaman.scores$interv_id),])
  }
  all.loci.scores = do.call(rbind, loci.scores)

  all.pvals = list()
  all.pvals[[1]] = list()
  all.pvals[[2]] = list()
  for (i in 1:2) {
    if (i == 1) {
      cur.tissue.loci.scores = filter(all.loci.scores, orig_selection & is_ecto)
      extended.bg.scores = filter(all.loci.scores, !orig_selection | !is_ecto)[,c('ecto_counts', 'meso_counts')]
    } else {
      cur.tissue.loci.scores = filter(all.loci.scores, orig_selection & is_meso)
      extended.bg.scores = filter(all.loci.scores, !orig_selection | !is_meso)[,c('ecto_counts', 'meso_counts')]
    }
    covs.tbl = table(rowSums(extended.bg.scores))
    best.fdrs = list()
    used.covs = as.numeric(setdiff(names(covs.tbl)[covs.tbl > 100], c('0', '1')))
    all.tissue.pvals = c()
    for (cur.cov in used.covs) {
      print(cur.cov)
      cur.bg.conts = extended.bg.scores[rowSums(extended.bg.scores) == cur.cov,,drop=F]
      cur.diff.bg = cur.bg.conts[,i] - cur.bg.conts[,3 - i]
      #cur.diff.thresh = quantile(cur.diff.dist, 0.95)
      cur.tested.conts = cur.tissue.loci.scores[rowSums(cur.tissue.loci.scores[,c('ecto_counts', 'meso_counts')]) == cur.cov,,drop=F]
      if (nrow(cur.tested.conts) == 0) {
        next
      }
      cur.diffs = cur.tested.conts[,i] - cur.tested.conts[,3 - i]
      cur.pvals = sapply(cur.diffs, function(cur.diff) mean(cur.diff <= cur.diff.bg))
      #cur.pvals = sapply(cur.diffs, function(cur.diff) mean(cur.diff < cur.diff.bg))
      names(cur.pvals) = cur.tested.conts$interv_id
      all.tissue.pvals = c(all.tissue.pvals, cur.pvals)
      #all.pvals[[i]] = all.tissue.pvals
      all.pvals[[i]][[cur.cov]] = cur.pvals
    }
  }

  png(file.path(fig.dir, 'ep_pvals_ecdf_ecto.png'), width=1000, height=1000)
  plot(ecdf(unlist(all.pvals[[1]][80:length(all.pvals[[1]])])), col=cls.col[1], lwd=10)
  abline(b=1, a=0, col=1, lwd=10)
  dev.off()

  png(file.path(fig.dir, 'ep_pvals_ecdf_meso.png'), width=1000, height=1000)
  plot(ecdf(unlist(all.pvals[[2]][80:length(all.pvals[[2]])])), col=cls.col[2], lwd=10)
  abline(b=1, a=0, col=1, lwd=10)
  dev.off()


}


