library(RColorBrewer)

ery.analysis <- function() {

  # boxplots for 2C and 2F were already created with the esc comparison
  fig.dir.emb.vs.ery = file.path(get.fig.dir(), 'emb_vs_ery_analysis')
  dir.create(fig.dir.emb.vs.ery, showWarnings=F)
  repl.data.dir = get.data.file.dir()

  exp_per_bin <<- get.or.create(EXP.PER.BIN.PATH, function() get.exp.per.bin(AB.SCORES.C.CLUSTERS.PATH))

  unnorm.total.cov = bin_cell_cov
  decay.metrics = get.common.decay.metrics()
  erys = get.erys()
  esc.cells = rownames(decay.metrics)[grep('scell', rownames(decay.metrics))]
  esc.cells = intersect(esc.cells, colnames(unnorm.total.cov))
  emb.cells = setdiff(rownames(sch_decay_metrics), erys)
  emb.cells = intersect(emb.cells, colnames(unnorm.total.cov))

  full.cluster.assignment = c(rep(1, length(emb.cells)), rep(2, length(erys)))
  names(full.cluster.assignment) = c(emb.cells, erys)
  unnorm.total.cov = unnorm.total.cov[,names(full.cluster.assignment)]

  num.umis = get.or.create(file.path(repl.data.dir, 'umi_counts'), function() get.num.umis(rownames(sch_decay_metrics)))
  print(c('median number of contacts:', median(num.umis[,1])))
  print(c('median fraction of trans contacts:', median(num.umis[,2] / num.umis[,1])))
  print(c('total number of contacts:', sum(num.umis[,1])))
  print(c('total fraction of trans contacts:', sum(num.umis[,2]) / sum(num.umis[,1])))
  print(c('fraction of cells enriched for 2-12Mb distance:', mean(sch_decay_metrics$f_mitotic_band > 0.2)))

  print(c('number of pErys:', length(erys)))

  early.late.ratio = colSums(unnorm.total.cov[early_late_bins$strict_A,]) / colSums(unnorm.total.cov[early_late_bins$strict_B,])
  decay.metrics$early_late_ratio = early.late.ratio[rownames(decay.metrics)]

  ab.form.path = file.path(repl.data.dir, 'ab_form_all_cells')
  ab.form = get.or.create(ab.form.path, function() total.ab.format(unique(bin_cell_ab$cell), bin_cell_ab))
  ab.form$atab = ab.form$atab[names(full.cluster.assignment),]
  ab.form$btab = ab.form$btab[names(full.cluster.assignment),]

  e10.ab.form.path = file.path(repl.data.dir, 'ab_form_e10')
  e10.ab.form = get.or.create(e10.ab.form.path, function() total.ab.format(unique(e10_bin_cell_ab$cell), e10_bin_cell_ab))

  deeply.covered.bins = rownames(unnorm.total.cov)[rowSums(unnorm.total.cov[, names(full.cluster.assignment)]) > 5e3]
  for (chrom in c('chrM', 'chrX', 'chrY')) {
    deeply.covered.bins = deeply.covered.bins[!grepl(chrom, deeply.covered.bins)]
  }
  chrom.cov = compute.chrom.cov(full.cluster.assignment, unnorm.total.cov[deeply.covered.bins,], decay.metrics)
  chrom.ab = compute.chrom.ab(full.cluster.assignment, ab.form, min.cov.discard=300, min.cov.per.clust=100, min.cov.cells=names(full.cluster.assignment))

  e10.cls = rep(1, nrow(e10_decay_metrics))
  e10.cells = rownames(e10_decay_metrics)
  names(e10.cls) = e10.cells
  e10.chrom.ab = compute.chrom.ab(e10.cls, e10.ab.form, min.cov.discard=50, min.cov.per.clust=50, min.cov.cells=e10.cells)

  common.bins = intersect(rownames(chrom.ab[[1]]), rownames(e10.chrom.ab[[1]]))
  png(file.path(fig.dir.emb.vs.ery, 'e10_a_scores_diff.png'))
  plot(chrom.ab[[1]][common.bins, 2] - chrom.ab[[1]][common.bins, 1], 
       e10.chrom.ab[[1]][common.bins, 1] - chrom.ab[[1]][common.bins, 1], col='darkblue', pch=19)
  grid(col='grey', lwd=1)
  dev.off()


  ab.scores.fil = ab_scores_c_clusters[, c('emb', 'ery')]
  ab.scores.fil = ab.scores.fil[apply(ab.scores.fil, 1, function(x) all(is.finite(x))),]
  print('fraction of > 30% changing bins in erys vs emb:')
  print(sum(abs(ab.scores.fil$emb - ab.scores.fil$ery) > 0.3) / nrow(ab.scores.fil))

  emb.ery.colors = c(1, 2)
  plot.basic.ery.vs.emb.scatters(ab_scores_c_clusters, repli_scores_c_clusters, exp_per_bin, fig.dir=fig.dir.emb.vs.ery)
  plot.basic.ery.vs.emb.stats(decay.metrics, erys, emb.cells, fig.dir.emb.vs.ery, cls.col=emb.ery.colors)
  plot.trans.ratio(num.umis, erys, emb.cells, fig.dir.emb.vs.ery, cls.col=emb.ery.colors)

  set.misha(SHAMAN.MISHA.PATH)
  ery.bin.clust.ret = clust_spat_a_ery(exp_per_bin, ab_scores_c_clusters, fig.dir=fig.dir.emb.vs.ery)
  plot.tal.gata.analysis('tal1_bigwig', clust.ret=ery.bin.clust.ret, fig.dir=fig.dir.emb.vs.ery)
  plot.tal.gata.analysis('gata1_bigwig', clust.ret=ery.bin.clust.ret, fig.dir=fig.dir.emb.vs.ery)

  emb.bin.clust.ret = clust_spat_a_ery(exp_per_bin, ab_scores_c_clusters, fig.dir=fig.dir.emb.vs.ery, tag='emb')

  track.names = paste0(c(EMB.POOL.TRACK.NAME, ERY.POOL.TRACK.NAME), '_equal_size')
  score.track.names = paste0(track.names, '_matshuff_scores')
  emb.ery.tracks.conts = lapply(score.track.names, function(track.name) gextract(track.name, gintervals.2d.all(), colnames='score'))

  e10.ery.tracks.conts = lapply(E10.POOL.TRACK.NAMES, function(track.name) gextract(track.name, gintervals.2d.all(), colnames='score'))

  pooled.hic.fig.dir = file.path(fig.dir.emb.vs.ery, 'cluster_pooled_hic_matrices')
  dir.create(pooled.hic.fig.dir, showWarnings=F)
  plot.pooled.shaman(emb.ery.tracks.conts, ery.bin.clust.ret, fig.dir=pooled.hic.fig.dir)

  # matching
  mcell_mc_plot_marks(mc_id='e9_recolor2', gset_id='e9_bs500f', mat_id='e9')
  exp.ret = get.or.create(file.path(repl.data.dir, 'atlases'), function() get.atlases(gene.md))
  atlas = exp.ret$atlas
  e9.exp = exp.ret$e9.exp
  
  matching.fig.dir = file.path(fig.dir.emb.vs.ery, 'matching')
  dir.create(matching.fig.dir, showWarnings=F)
  ery.metacell.threshold = list(c(0.1, 0.1), c(0.08, 0.13))
  exp.data = list(list(log2(e9.exp@mc_fp), e9.exp@colors, 'e9', e9.exp, 'e9_recolor2'), list(log2(atlas@mc_fp), atlas@colors, 'atlas', atlas, 'emb_gotg_bs500f'))
  for (i in 1:length(exp.data)) {
    lfp = exp.data[[i]][[1]]
    colors = exp.data[[i]][[2]]
    name = exp.data[[i]][[3]]
    mc.obj = exp.data[[i]][[4]]
    mc.name = exp.data[[i]][[5]]

    cov.matching.ret = metacell.matching(lfp, chrom.cov, 1:2)
    ab.matching.ret = metacell.matching(lfp, chrom.ab, 1:2)
    plot.metacell.matching(cov.matching.ret, colors, matching.fig.dir, paste0('cov_matching_', name), atlas=mc.obj, mc.name=mc.name)
    plot.metacell.matching(ab.matching.ret, colors, matching.fig.dir, paste0('ab_matching_', name), atlas=mc.obj, mc.name=mc.name)
  }



  ery.bins = data.frame(chrom=c('chr2', 'chr4', 'chr16', 'chr16', 'chr15', 'chr14', 'chr13', 'chr12'), mid=c(45e6, 6e6, 77e6, 65e6, 43e6, 104e6, 84e6, 37e6))
  ery.genes = c('Hbb-y', 'Cpox')

  plot.dir.name = file.path(fig.dir.emb.vs.ery, 'hic_matrices')
  dir.create(plot.dir.name, showWarnings=F)
  e10.plot.dir.name = file.path(fig.dir.emb.vs.ery, 'e10_hic_matrices')
  dir.create(e10.plot.dir.name, showWarnings=F)


  ab.form.hres.path = file.path(repl.data.dir, 'ab_form_hres_erys')
  ab.form.hres = get.or.create(ab.form.hres.path, function() total.ab.format(names(full.cluster.assignment), bin_cell_ab_hres))
  chrom.ab.hres = get.or.create(file.path(repl.data.dir, 'chrom_ab_hres_erys'),
                  function() compute.chrom.ab(full.cluster.assignment, ab.form.hres, min.cov.discard=40, min.cov.per.clust=40, min.cov.cells=names(full.cluster.assignment)))

  gene.md = get.gene.metadata.jumps()
  all.interesting.bins = genes.to.intervs(ery.genes, gene.md, F)
  all.interesting.bins$geneSymbol = as.character(all.interesting.bins$geneSymbol)
  all.interesting.bins$clust = 3

  for (i in 1:nrow(all.interesting.bins)) {
    interv = all.interesting.bins[i, ]
    forced.interv = gintervals.force_range(data.frame(chrom=interv$chrom, start=interv$start, end=interv$end))
    stopifnot(nrow(interv) == nrow(forced.interv))
    forced.interv$orig_start = interv$orig_start
    forced.interv$tss = interv$tss
    forced.interv$geneSymbol = interv$geneSymbol
    forced.interv$clust = interv$clust
    interv = forced.interv
    interv$chrom = as.character(interv$chrom)
    file.name = paste0(interv$geneSymbol, '_', interv$chrom, '_',
                         trimws(format(interv$start, scientific=F)), '.png')
    plot.file.name = file.path(plot.dir.name, file.name)
    plot.whole.chrom(chrom.ab.hres, chrom.ab.hres, chrom.ab.hres, emb.ery.tracks.conts, as.numeric(substr(interv$chrom, 4, 6)), plot.file.name, coord.range=c(interv$start, interv$end), cls.col=emb.ery.colors, smooth.cov=F, height.vec=c(3, 3, 1.5), width=6, sel.attr.plots=c(2), show.coord.text=F, smooth.ab=T, rotate=T, mark.coord=interv$tss)

    e10.file.name = file.path(e10.plot.dir.name, paste0(interv$geneSymbol, '.png'))
    plot.whole.chrom(chrom.ab.hres, chrom.ab.hres, chrom.ab.hres, e10.ery.tracks.conts, as.numeric(substr(interv$chrom, 4, 6)), e10.file.name, coord.range=c(interv$start, interv$end), cls.col=emb.ery.colors, smooth.cov=F, height.vec=c(3, 3, 1.5), width=6, sel.attr.plots=c(2), show.coord.text=F, smooth.ab=T, rotate=T, mark.coord=interv$tss)
  }

  # and now ery bins
  de.ab.bins = lapply(1:2, function(i) get.ab.up.bins(chrom.ab, i, 3 - i, 0.2))
  interesting.bins.list = de.ab.bins
  interesting.bins.df.list = lapply(interesting.bins.list, group.bins)
  for (i in 1:length(interesting.bins.df.list)) {
    interesting.bins.df.list[[i]]$clust = i
  }
  all.interesting.bins = interesting.bins.df.list[[2]]
  all.interesting.bins$orig_start = as.numeric(all.interesting.bins$start)
  all.interesting.bins$start = pmax(as.numeric(all.interesting.bins$start) - 2e6, 0)
  all.interesting.bins$end = as.numeric(all.interesting.bins$end) + 2e6

  for (i in 1:nrow(all.interesting.bins)) {
    interv = all.interesting.bins[i, ]
    forced.interv = gintervals.force_range(data.frame(chrom=interv$chrom, start=interv$start, end=interv$end))
    stopifnot(nrow(interv) == nrow(forced.interv))
    forced.interv$orig_start = interv$orig_start
    forced.interv$seq_len = interv$seq_len
    forced.interv$clust = interv$clust
    forced.interv$is_orphan = interv$is_orphan
    forced.interv$is_orphan2 = interv$is_orphan2
    interv = forced.interv
    interv$chrom = as.character(interv$chrom)

    if (!(any(ery.bins$chrom == interv$chrom & between(ery.bins$mid, interv$start, interv$end)))) {
      next
    }

    file.name = paste0(interv$chrom, '_', trimws(format(interv$start, scientific=F)), '.png')
    plot.file.name = file.path(plot.dir.name, file.name)

    plot.whole.chrom(chrom.ab.hres, chrom.ab.hres, chrom.ab.hres, emb.ery.tracks.conts, as.numeric(substr(interv$chrom, 4, 6)), plot.file.name, coord.range=c(interv$start, interv$end), cls.col=emb.ery.colors, smooth.cov=F, height.vec=c(3, 3, 1.5), width=6, sel.attr.plots=c(2), show.coord.text=F, smooth.ab=T, rotate=T)

    e10.file.name = file.path(e10.plot.dir.name, paste0(interv$orig_start, '.png'))
    plot.whole.chrom(chrom.ab.hres, chrom.ab.hres, chrom.ab.hres, e10.ery.tracks.conts, as.numeric(substr(interv$chrom, 4, 6)), e10.file.name, coord.range=c(interv$start, interv$end), cls.col=emb.ery.colors, smooth.cov=F, height.vec=c(3, 3, 1.5), width=6, sel.attr.plots=c(2), show.coord.text=F, smooth.ab=T, rotate=T)
  }

}

plot.trans.ratio <- function(num.umis, erys, emb.cells, fig.dir, cls.col=c(1, 2)) {
  erys = intersect(erys, rownames(num.umis))
  emb.cells = intersect(emb.cells, rownames(num.umis))
  trans.frac = num.umis[, 2] / num.umis[, 1]
  xlim = range(trans.frac)
  png(file.path(fig.dir, 'emb_ery_trans.png'), width=400, height=400)
  plot(density(trans.frac[erys]), xlim=xlim, col=cls.col[2], lwd=8, xlab='', ylab='', main='')
  lines(density(trans.frac[emb.cells]), xlim=xlim, lwd=8, col=cls.col[1])
  grid(col='grey', lwd=0.5)
  dev.off()
}

plot.basic.ery.vs.emb.stats <- function(decay.metrics, erys, emb.cells, fig.dir, cls.col=c(1, 2, 3)) {
  decay.metrics = decay.metrics[c(emb.cells, erys),]
  col=c(rep(cls.col[1], length(emb.cells)), rep(cls.col[2], length(erys)))
  png(file.path(fig.dir, 'far_vs_mitotic.png'), width=400, height=400)
  plot(decay.metrics[,'f_mitotic_band'], 1 - decay.metrics[,'near_f'], col=col, pch=19, cex=1.2, xlab='', ylab='')
  dev.off()

  print('ery fraction of contacts above 2Mb:')
  print(quantile(1 - decay.metrics[erys, 'near_f'], (0:20)/20))

  png(file.path(fig.dir, 'tight_far_vs_far.png'), width=400, height=400)
  plot(decay.metrics[,'far_tightness'], 1 - decay.metrics[,'near_f'], col=col, pch=19, cex=1.2, xlab='', ylab='')
  dev.off()

  print('ery far tightness:')
  print(quantile(decay.metrics[erys, 'far_tightness'], (0:20)/20))

  png(file.path(fig.dir, 'ery_short_vs_repli.png'), height=300, width=300)
  is.ery = rownames(decay.metrics) %in% erys
  plot(log2(decay.metrics$early_late_ratio)[!is.ery], decay.metrics$f_near_band[!is.ery],  pch=19, cex=1, col='black', ylab="", xlab="repli-score", main="")
  points(log2(decay.metrics$early_late_ratio)[is.ery], decay.metrics$f_near_band[is.ery], col='red', pch=19, cex=1)
  grid(col='grey', lwd=1)
  dev.off()

  png(file.path(fig.dir, 'repli_score_dist.png'), width=400, height=400)
  #min.early = min(decay.metrics$early_late_ratio)
  min.early = 1
  ery.density = density(log2(decay.metrics[erys, 'early_late_ratio']) / min.early)
  emb.density = density(log2(decay.metrics[emb.cells, 'early_late_ratio']) / min.early)
  plot(emb.density, col=cls.col[1], type='l', lwd=8, main='')
  lines(ery.density, col=cls.col[2], type='l', lwd=8)
  grid(col='grey', lwd=0.5)
  dev.off()
}

plot.basic.ery.vs.emb.scatters <- function(ab.scores, rep.scores, exp.per.bin, fig.dir=get.fig.dir()) {

  load(EMB.G1.CELLS.PATH)
  high.cov.loci = rownames(bin_cell_cov_hres_ds)[rowSums(bin_cell_cov_hres_ds[, ery_G1_cells]) > 15]
  stopifnot(all(rownames(bin_cell_cov_hres_ds) == rownames(rep.scores)))
  stats.to.plot = list(list(cur.data=ab.scores, name='ab', lim=c(0, 1), cex=0.3),
                  list(cur.data=rep.scores[high.cov.loci,], name='repli', lim=c(-1.6, 1.6), cex=0.3),
                  list(cur.data=log2(exp.per.bin), name='exp', lim=c(-23, -6), cex=0.5))
  for (i in 1:length(stats.to.plot)) {
    png(file.path(fig.dir, sprintf('emb_ery_scatter_%s.png', stats.to.plot[[i]]$name)))
    # hack to plot the points under the grid
    cur.data = stats.to.plot[[i]]$cur.data
    cur.name = stats.to.plot[[i]]$name
    lim = stats.to.plot[[i]]$lim
    cex = stats.to.plot[[i]]$cex
    if (is.null(lim)) {
      lim = range(c(cur.data$emb, cur.data$ery), finite=T)
    }
    plot(cur.data$emb, cur.data$ery, pch=19, cex=cex, col='white', ylim=lim, xlim=lim)
    grid(col='black')
    points(cur.data$emb, cur.data$ery, pch=19, cex=cex, col='darkblue')
    abline(b=1, a=0, col='red', lwd=3, lty='dashed')
    if (cur.name == 'ab') {
      abline(b=1, a=0.3, col='red', lwd=2, lty='dashed')
      abline(b=1, a=-0.3, col='red', lwd=2, lty='dashed')
    }
    dev.off()
  }
}

get.bs.scores <- function(tname) {
  vtrack1 = gvtrack.create('vtrack1', tname, 'global.percentile')
  intervs= gextract('vtrack1', gintervals.all())
  intervs$score = -log2(1 - intervs$vtrack1)
  intervs = filter(intervs, score > 8)
  return(intervs)
}

plot.tal.gata.analysis <- function(tname, clust.ret=NULL, fig.dir=get.fig.dir(), data.dir=get.data.file.dir()) {
  ab.scores = ab_scores_c_clusters
  if (is.null(clust.ret)) {
    clust.ret = clust_spat_a_ery(exp_per_bin, ab.scores, fig.dir=fig.dir)
  }
  diff.bins = names(clust.ret$clusters)
  splitted = strsplit(diff.bins, '_')
  intervs = data.frame(chrom=sapply(splitted, function(x) x[1]), start=as.numeric(sapply(splitted, function(x) x[2])) - 4e4)
  intervs$end = intervs$start + 12e4 - 1
  intervs = gintervals.force_range(intervs)
  rownames(intervs) = diff.bins

  all.exp = do.call(rbind, lapply(1:nrow(intervs), function(i) {
    bin.names = get.all.bins(intervs[i, 'chrom'], as.numeric(intervs[i, 'start']), as.numeric(intervs[i, 'end']), 4e4)[[1]]
    bin.names = intersect(bin.names, rownames(exp_per_bin))
    return(c(max(exp_per_bin[bin.names, 'emb'], na.rm=T), max(exp_per_bin[bin.names, 'ery'], na.rm=T)))
  }))

  emb.exp = all.exp[,1]
  ery.exp = all.exp[,2]
  exp.diff = log2(emb.exp) - log2(ery.exp)
  exp.bin.type = ifelse(pmax(log2(emb.exp), log2(ery.exp)) < -18, 'inactive', 
                 ifelse(exp.diff < -2, 'ery_active',
                 ifelse(exp.diff > 2, 'emb_active', 'neutral')))
  exp.bin.type = factor(exp.bin.type, levels=c('inactive', 'neutral', 'ery_active', 'emb_active'))
  intervs$type = exp.bin.type
  intervs.splitted = lapply(split(1:nrow(intervs), intervs$type), function(x) intervs[x,])
  bs.scores = get.or.create(file.path(data.dir, paste0(tname, '_bs_scores')), function() get.bs.scores(tname))
  cols = c('black', 'grey', 'red', 'blue', 'saddlebrown')
  plot.orphan.bs.score.cdf(intervs.splitted, bs.scores, rownames(ab_scores_c_clusters), file.path(fig.dir, paste0(tname, '_exp_score_T.png')), is.max=T, col=cols)


  intervs.splitted2 = lapply(split(1:nrow(intervs), clust.ret$clusters), function(x) intervs[x,])
  cols2 = c(brewer.pal(length(table(clust.ret$clusters)), 'Set3'), 'saddlebrown')
  plot.orphan.bs.score.cdf(intervs.splitted2, bs.scores, rownames(ab_scores_c_clusters), file.path(fig.dir, paste0(tname, '_clusts_score_T.png')), is.max=T, col=cols2)
}


clust_spat_a_ery = function(exp_per_bin, ab_scores_c_clusters, fig.dir=get.fig.dir(), tag="ery", seed=42)
{
	set.seed(42)
	ab40 = ab_scores_c_clusters

	ab40$chrom = sub("_\\d+","",rownames(ab40))
	ab40$start = as.numeric(sub("chr.+_","",rownames(ab40)))
	ab40$end = ab40$start + 40000 - 1

	ab40s = ab40[order(as.numeric(as.factor(ab40$chrom))*1e+9+ab40$start),]
	ab40s$end = ab40s$end+1
	ab40s$id = 1:nrow(ab40s)

	f_na = is.na(ab40$ery-ab40$emb)
	f = !f_na & (ab40s$ery-ab40s$emb) > 0.35; 
	if(tag == "emb") {
	f = !f_na & (ab40s$ery-ab40s$emb) < -0.35; 
	}
	f[is.na(f)]=F

	hits = gintervals.force_range(ab40s[f,c("chrom","start", "end", "emb","ery","id")])
	peaks = gintervals.canonic(hits)
	peaks$peak_id = 1:nrow(peaks)

	hits_i = gintervals.neighbors(hits, peaks)

	dlt = hits_i$ery-hits_i$emb

	if (tag == "emb") {
	  peak_max = tapply(dlt, hits_i$peak_id, min)
	} else {
	  peak_max = tapply(dlt, hits_i$peak_id, max)
	}

	hits_i$peak_max = peak_max[hits_i$peak_id]
	
	f_top = hits_i$peak_max == dlt

	top_ids = hits_i$id[f_top]
	print('number of high ery loci')
	print(length(top_ids))

	a_mat_ery = matrix(ab40s$ery[top_ids-10], ncol=1)
	a_mat_emb = matrix(ab40s$emb[top_ids-10], ncol=1)
	for(i in -9:10) {
		a_mat_ery = cbind(a_mat_ery, ab40s$ery[top_ids+i])
		a_mat_emb = cbind(a_mat_emb, ab40s$emb[top_ids+i])
	}

	rna_emb = log2(1e-5+exp_per_bin$emb)
	rna_ery = log2(1e-5+exp_per_bin$ery)
	names(rna_emb) = rownames(exp_per_bin)
	names(rna_ery) = rownames(exp_per_bin)

	rna_mat_ery = matrix(rna_ery[rownames(ab40s)[top_ids-10]], ncol=1)
	rna_mat_emb = matrix(rna_emb[rownames(ab40s)[top_ids-10]], ncol=1)
	for(i in -9:10) {
		rna_mat_ery = cbind(rna_mat_ery, rna_ery[rownames(ab40s)[top_ids+i]])
		rna_mat_emb = cbind(rna_mat_emb, rna_emb[rownames(ab40s)[top_ids+i]])
	}
	rna_mat_ery[is.na(rna_mat_ery)] = -30
	rna_mat_emb[is.na(rna_mat_emb)] = -30

	polarity = rowMeans(a_mat_ery[,1:10],na.rm=T)- rowMeans(a_mat_ery[,11:20],na.rm=T)
	polarity[is.nan(polarity)]=0
	a_mat_ery_p = a_mat_ery
	a_mat_emb_p = a_mat_emb
	a_mat_ery_p[polarity<0,] = a_mat_ery[polarity<0,21:1]
	a_mat_emb_p[polarity<0,] = a_mat_emb[polarity<0,21:1]

	rna_mat_ery_p = rna_mat_ery
	rna_mat_emb_p = rna_mat_emb
	rna_mat_ery_p[polarity<0,] = rna_mat_ery[polarity<0,21:1]
	rna_mat_emb_p[polarity<0,] = rna_mat_emb[polarity<0,21:1]


	km = tglkmeans::TGL_kmeans(cbind(a_mat_ery_p,a_mat_emb_p), k=8, id_column=F, reorder_func=mean, seed=42)
	bins_hclust = hclust(tgs_dist(cbind(a_mat_ery_p,a_mat_emb_p)), "complete")
	hc_ord = rep(0, length(km$cluster))
	hc_ord[bins_hclust$order] = 1:length(km$cluster)

	ord = order(km$cluster+1e-5*hc_ord)
	bin.names = rownames(ab40s)[top_ids]
	names(km$cluster) = bin.names
	clusters = km$cluster[ord]

        shades = shaman_score_pal()
	bin.cls.col = brewer.pal(length(table(km$cluster)), 'Set3')
	ery.a.spat.abs = cbind(a_mat_ery_p, a_mat_emb_p)[ord,]
	rownames(ery.a.spat.abs) = names(clusters)
	pheatmap::pheatmap(ery.a.spat.abs, 
			   cluster_cols=F, cluster_rows=F, 
			   color=shades, breaks = seq(0,1,l=200),
			   filename = file.path(fig.dir, sprintf("%s_A_spat_abs.png",tag)), width=10,height=14,
			   show_rownames=F, show_colnames=F, annotation_names_row=F, annotation_names_col=F, annotation_legend=F,
			   annotation_row=as.data.frame(clusters), annotation_colors=list(clusters=bin.cls.col))

        drna = (rna_mat_emb_p-rna_mat_ery_p)[ord,]
	all.max.exp = apply(cbind(rna_mat_ery_p[ord,], rna_mat_emb_p[ord,]), 1, max)
	print('fraction of expressed loci')
	print(sum(all.max.exp > -16) / length(all.max.exp)) 
	drna.lim = pmax(pmin(drna, 5), -5)
	pheatmap::pheatmap(drna.lim,
					cluster_cols=F, cluster_rows=F, 
					color=c(colorRampPalette(c("darkred", "white", "darkblue"))(102)), 
					breaks = c(-5,seq(-2,2,l=100),5),
			                show_rownames=F, show_colnames=F, annotation_names_row=F, annotation_names_col=F, annotation_legend=F,
					filename = file.path(fig.dir, sprintf("%s_A_spat_drna.png",tag)), width=5,height=14)
	
	should_reverse = polarity < 0
	names(should_reverse) = names(km$cluster)

        return(list(clusters=km$cluster, should_reverse=should_reverse))
}

plot.orphan.bs.score.cdf <- function(intervs.list, bs, bins, fig.path, is.max=F, bin.size=4e4, col=NULL, draw.rand=T, num.rand=5e3) {
  all.cdf = lapply(intervs.list, function(intervs) bs.score.cdf(bs, intervs, is.max))

  rand.bins = sample(bins, num.rand)
  splitted = strsplit(rand.bins, '_')
  rand.intervs3 = data.frame(chrom=sapply(splitted, function(x) x[1]), start=sapply(splitted, function(x) as.numeric(x[2])))
  rand.intervs3$end = rand.intervs3$start + bin.size
  rand.intervs3= gintervals.force_range(rand.intervs3)

  rand3.cdf = bs.score.cdf(bs, rand.intervs3, is.max)
  max.value = max(sapply(c(all.cdf, rand3.cdf), function(x) quantile(x, 1)))
  min.value = min(sapply(c(all.cdf, rand3.cdf), function(x) quantile(x, 0)))
  if (is.null(col)) {
    col = 1:(length(all.cdf) + 1)
  }
  png(fig.path, width=800, height=800)
  plot(all.cdf[[1]], col=col[1], lwd=8, do.points=F, xlim=c(min.value, max.value), verticals=T)
  for (i in 2:length(all.cdf)) {
    lines(all.cdf[[i]], col=col[i], lwd=8, do.points=F, verticals=T)
  }
  if (draw.rand) {
    lines(rand3.cdf, col=col[length(all.cdf) + 1], lwd=8, do.points=F, verticals=T)
  }
  dev.off()
}

bs.score.cdf <- function(bs, intervs, is.max=F) {
  if (is.max) {
    func = max
  } else {
    func = mean
  }
  intervs$bin.names = paste0(intervs$chrom, '_', trimws(format(intervs$start, scientific=F)))
  all.dists = gintervals.neighbors(intervs, bs, maxneighbors=nrow(bs), maxdist=0)
  mean.scores = tapply(1:nrow(all.dists), all.dists$bin.names, function(x) func(all.dists[x, 'score']))
  return(ecdf(mean.scores))
}

get.atlases <- function(gene.md) {
  atlas.obj = get.or.create(ATLAS.PATH, function() error('Did not find expression atlas!'))
  e9.obj = get.or.create(E9.EXP.PATH, function() error('Did not find e9 expression data!'))
  esc.obj = get.or.create(ESC.EXP.PATH, function() error('Did not find esc expression data!'))

  matching = match.gene.names(rownames(atlas.obj@e_gc), gene.md$geneSymbol)
  rownames(atlas.obj@e_gc) = gene.md$geneSymbol[matching]
  atlas.obj@e_gc = atlas.obj@e_gc[!is.na(matching),]
  rownames(atlas.obj@mc_fp) = gene.md$geneSymbol[matching]
  atlas.obj@mc_fp = atlas.obj@mc_fp[!is.na(matching),]

  matching = match.gene.names(rownames(e9.obj@e_gc), gene.md$geneSymbol)
  rownames(e9.obj@e_gc) = gene.md$geneSymbol[matching]
  e9.obj@e_gc = e9.obj@e_gc[!is.na(matching),]
  rownames(e9.obj@mc_fp) = gene.md$geneSymbol[matching]
  e9.obj@mc_fp = e9.obj@mc_fp[!is.na(matching),]

  matching = match.gene.names(rownames(esc.obj@e_gc), gene.md$geneSymbol)
  rownames(esc.obj@e_gc) = gene.md$geneSymbol[matching]
  esc.obj@e_gc = esc.obj@e_gc[!is.na(matching),]
  rownames(esc.obj@mc_fp) = gene.md$geneSymbol[matching]
  esc.obj@mc_fp = esc.obj@mc_fp[!is.na(matching),]

  #load('~/proj/e9schic/data/e9_lfp.Rda', ver=T)
  return(list(atlas=atlas.obj, e9.exp=e9.obj, esc.exp=esc.obj))
}

match.gene.names <- function(genes1, genes2) {
  genes1.splitted = strsplit(genes1, ';')
  genes2.splitted = strsplit(genes2, ';')
  genes2.indices = unlist(lapply(1:length(genes2), function(i) {ret = rep(i, length(genes2.splitted[[i]]));names(ret) = genes2.splitted[[i]];return(ret)}))

  all.matches = sapply(genes1.splitted, function(gene1) {

    indices = genes2.indices[gene1]
    if (all(is.na(indices))) {
      return(NA)
    } else {
      return(min(indices, na.rm=T))
    }

  })
  return(all.matches)
}

metacell.matching <- function(lfp, chrom.ret, clusters.to.compare, gene.md=NULL, max.num.genes=50, de.genes=NULL) {
  if (is.null(gene.md)) {
    gene.md = get.gene.metadata.jumps()
  }
  gene.md$bin.name = paste0(gene.md$chrom, '_', trimws(gene.md$bin))
  gene.md = gene.md[!duplicated(gene.md$geneSymbol),]
  nclust = length(clusters.to.compare)
  #de.genes = lapply(1:nmc, function(j) rownames(lfp)[lfp[,j] > 0.5])
  all.matchings = list()
  de.genes.per.clust = list()
  mc.indices.list = list()
  if (is.null(de.genes)) {
    nmc = ncol(lfp)
    de.genes = lapply(1:nmc, function(j) intersect(names(sort(lfp[,j], decreasing=T))[1:max.num.genes], 
                                                   rownames(lfp)[lfp[,j] > 0.5]))
  } else {
    nmc = length(de.genes)
  }

  for (i in 1:nclust) {
    if (nclust == 2) {
      chrom.ratio = log(chrom.ret[[1]][,i] / chrom.ret[[1]][,3 - i])
    } else {
      chrom.ratio = log(chrom.ret[[1]][,i] / chrom.ret[[2]][,i])
    }
    mc.match = sapply(1:nmc, function(j) {
      gene.matching = match(de.genes[[j]], gene.md$geneSymbol)
      return(mean(chrom.ratio[gene.md[gene.matching, 'bin.name']], na.rm=T))
    })
    all.matchings[[i]] = mc.match
    top.mc.indices = order(mc.match, decreasing=T)[1:3]
    mc.indices.list[[i]] = top.mc.indices
    de.genes.per.clust[[i]] = unique(unlist(de.genes[top.mc.indices]))
  }
  return(list(de.genes=de.genes.per.clust, mc.de.genes=de.genes, all.matchings=all.matchings, top.mc.indices=mc.indices.list))
}

plot.metacell.matching <- function(matching.ret, colors, fig.dir, fig.name, atlas=NULL, mc.name=NULL) {
  nclust = length(matching.ret$all.matchings)

  shades = colorRampPalette(c("darkblue", "blue", "white","red", "yellow"))(200)
  for (i in 1:nclust) {
    matching = matching.ret$all.matchings[[i]]
    cols = shades[max.min.rescale(matching, 1, 200)]
    
    my.metacell.plot(mc.name, cols, fig.dir, paste0(fig.name, '_', i))

    png(file.path(fig.dir, paste0(fig.name, '_', i, '_legend.png')))
    image.plot(matrix(c(max(matching), min(matching))), col=shades)
    dev.off()
  }

  for (i in 1:nclust) {
    de.genes = matching.ret$de.genes[[i]]
    matching = colSums(atlas@mc_fp[de.genes,] > 1.2)
    cols = shades[max.min.rescale(matching, 1, 200)]
    
    my.metacell.plot(mc.name, cols, fig.dir, paste0(fig.name, '_exp_', i))
  }
  my.metacell.plot(mc.name, atlas@colors, fig.dir, paste0(fig.name, '_orig_cols'))
  my.metacell.plot(mc.name, atlas@colors, fig.dir, paste0(fig.name, '_orig_cols_only_mc'), only.mc=T)
}

my.metacell.plot <- function (mc2d_id, cols, fig_dir, fig_suffix, legend_pos = "topleft", only.mc=F, mc.text=F)
{
    mcp_2d_height = 2400
    mcp_2d_width = 2400
    mcp_2d_plot_key = F
    mcp_2d_cex = 3
    mcp_2d_legend_cex = 1
    mc2d = scdb_mc2d(mc2d_id)
    if (is.null(mc2d)) {
        stop("missing mc2d when trying to plot, id ", mc2d_id)
    }
    mc = scdb_mc(mc2d@mc_id)
    if (is.null(mc)) {
        stop("missing mc in mc2d object, id was, ", mc2d@mc_id)
    }
    #fig_nm = scfigs_fn(mc2d_id, "2d_proj")
    fig_nm = scfigs_fn(mc2d_id, fig_suffix, dir=fig_dir)
    png(fig_nm, width = mcp_2d_width, height = mcp_2d_height)
    #cols = mc@colors
    cols[is.na(cols)] = "gray"
    plot(mc2d@sc_x, mc2d@sc_y, pch = 19, col = cols[mc@mc[names(mc2d@sc_x)]], cex=mcp_2d_cex)
    fr = mc2d@graph$mc1
    to = mc2d@graph$mc2
    segments(mc2d@mc_x[fr], mc2d@mc_y[fr], mc2d@mc_x[to], mc2d@mc_y[to])
    if (!only.mc) {
      points(mc2d@mc_x, mc2d@mc_y, cex = 3 * mcp_2d_cex, col = "black",
          pch = 21, bg = cols)
      if (mc.text) {
        text(mc2d@mc_x, mc2d@mc_y, 1:length(mc2d@mc_x), cex = mcp_2d_cex)
      }
    }
    dev.off()
}

max.min.rescale <- function(value, wanted.min, wanted.max) {
  max.val = max(value)
  min.val = min(value)
  
  newvalue = (wanted.max - wanted.min) / (max.val - min.val) * (value - min.val) + wanted.min
  return(newvalue)
}

get.ab.differential.bins <- function(chrom.ab.ret, cluster, compared.to, threshold=0.1) {
  #ab.lfc = log(chrom.ab.ret[[1]][,cluster] / chrom.ab.ret[[2]][,cluster])
  chrom.ab = chrom.ab.ret[[1]][apply(chrom.ab.ret[[1]], 1, function(x) all(!is.na(x))),]
  ab.diff = chrom.ab[,cluster] - chrom.ab[,compared.to]
  ab.bins = names(ab.diff)[abs(ab.diff) > threshold]
  ab.sign = ifelse(ab.diff[ab.bins] > 0, 1, -1)
  return(list(bins=ab.bins, sign=ab.sign))
}

get.bins.sorted <- function(bin.names, chrom) {
  splitted = strsplit(bin.names, '_')
  is.in.chrom = sapply(1:length(splitted), function(i) splitted[[i]][1] == chrom)
  if (sum(is.in.chrom) == 0) {
    return(c())
  }
  ordered.bins = paste(chrom, sort(sapply(splitted[is.in.chrom], function(lst) as.numeric(lst[[2]]))), sep='_')
  return(ordered.bins)
}
