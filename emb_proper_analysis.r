library(zoo)
library(umap)
library(sm)
library(pheatmap)
library(caret)
library(reshape2)


get.lres.repli.scores <- function(clusters) {
  if (!file.exists(REPLI.SCORES.EMB.CLUSTERS.LOW.RES.PATH)) {
    load(EMB.G1.CELLS.PATH)
    cls1_loc_SG1 = log2(rowMeans(bin_cell_cov_ds[, names(clusters)[clusters == 1]]) / rowMeans(bin_cell_cov_ds[, emb_G1_cells]))
    cls2_loc_SG1 = log2(rowMeans(bin_cell_cov_ds[, names(clusters)[clusters == 2]]) / rowMeans(bin_cell_cov_ds[, emb_G1_cells]))
    cls3_loc_SG1 = log2(rowMeans(bin_cell_cov_ds[, names(clusters)[clusters == 3]]) / rowMeans(bin_cell_cov_ds[, emb_G1_cells]))

    repli_scores_emb_clusters_lres = data.frame(repli1=cls1_loc_SG1, repli2=cls2_loc_SG1, repli3=cls3_loc_SG1)
    write.table(repli_scores_emb_clusters_lres, file=REPLI.SCORES.EMB.CLUSTERS.LOW.RES.PATH)
  }
  repli_scores_emb_clusters_lres = read.table(REPLI.SCORES.EMB.CLUSTERS.LOW.RES.PATH)
  return(repli_scores_emb_clusters_lres)
}

get.splitted.cls <- function() {
  # this file was created based on the 3 cluster solution of the model
  # but cluster 2 cells were further splitted into those with "strong" assocation
  # to the the cluster, and two smaller subclusters
  splitted.cls.raw = read.table(CLUSTERING5.PATH)
  splitted.cls = splitted.cls.raw[,1]
  names(splitted.cls) = rownames(splitted.cls.raw)
  return(splitted.cls)
}

get.splitted.repli.scores <- function() {
  if (!file.exists(REPLI.SCORES.EMB.CLUSTERS.SPLITTED.PATH)) {
    load(EMB.G1.CELLS.PATH)
    splitted.cls = get.splitted.cls()
    cls1_loc_SG1 = log2(rowMeans(bin_cell_cov_hres_ds[, names(splitted.cls)[splitted.cls == 1]]) / rowMeans(bin_cell_cov_hres_ds[, emb_G1_cells]))
    cls2_loc_SG1 = log2(rowMeans(bin_cell_cov_hres_ds[, names(splitted.cls)[splitted.cls == 2]]) / rowMeans(bin_cell_cov_hres_ds[, emb_G1_cells]))
    cls3_loc_SG1 = log2(rowMeans(bin_cell_cov_hres_ds[, names(splitted.cls)[splitted.cls == 3]]) / rowMeans(bin_cell_cov_hres_ds[, emb_G1_cells]))

    repli_scores_emb_splitted = data.frame(C21=cls1_loc_SG1, C22=cls2_loc_SG1, C23=cls3_loc_SG1)
    write.table(repli_scores_emb_splitted, file=REPLI.SCORES.EMB.CLUSTERS.SPLITTED.PATH)
  }
  repli_scores_emb_splitted = read.table(REPLI.SCORES.EMB.CLUSTERS.SPLITTED.PATH)
  return(repli_scores_emb_splitted)
}

get.splitted.ab.scores <- function() {
  if (!file.exists(AB.SCORES.EMB.CLUSTERS.SPLITTED.PATH)) {
    splitted.cls = get.splitted.cls()

    # requires prior loading of bin_cell_ab_hres
    ab.form.hres = total.ab.format(names(splitted.cls), bin_cell_ab_hres)
    chrom.ab.hres = compute.chrom.ab(splitted.cls, ab.form.hres, min.cov.discard=40, min.cov.per.clust=40)[[1]]
    ab_scores_emb_splitted = data.frame(C21=chrom.ab.hres[,1], C22=chrom.ab.hres[,2], C23=chrom.ab.hres[,3])
    write.table(ab_scores_emb_splitted, file=AB.SCORES.EMB.CLUSTERS.SPLITTED.PATH)
  }
  ab_scores_emb_splitted = read.table(AB.SCORES.EMB.CLUSTERS.SPLITTED.PATH)[,1:3]
  return(ab_scores_emb_splitted)
}


emb.proper.analysis <- function() {
  decay.metrics = get.common.decay.metrics()
  good.clusters = c(1, 3)

  embryo.cls.path = file.path(repl.data.dir, 'embryo_cls')
  embryo.cls.ret = get.or.create(embryo.cls.path, cluster.embryo.cells)
  model.params = embryo.cls.ret$model.params
  unnorm.total.cov = embryo.cls.ret$unnorm.total.cov
  clusters = model.to.clusters(model.params)

  fig.dir = get.fig.dir()
  repl.fig.dir = file.path(fig.dir, 'repl_model')
  dir.create(repl.fig.dir, showWarnings=F)
  repl.data.dir = get.data.file.dir()


  # and now classification - extend clustering into non replicating cells
  emb.decay.metrics = sch_decay_metrics
  cells.in.ab.form = setdiff(intersect(colnames(bin_cell_cov), rownames(emb.decay.metrics)), get.erys())
  ab.form.path = file.path(repl.data.dir, 'ab_form_emb_proper')
  ab.form = get.or.create(ab.form.path, function() total.ab.format(cells.in.ab.form, bin_cell_ab))
  de.bins = get.de.replicating.bins(model.params, good.clusters)
  g1.norm.total.cov = unnorm.total.cov / model.params$bin.probs
  g1.norm.total.cov = t(t(g1.norm.total.cov) / colSums(g1.norm.total.cov))

  repli_scores_emb_clusters_lres = get.lres.repli.scores(clusters)
  plot.de.bins(de.bins, ab.form, g1.norm.total.cov, repli_scores_emb_clusters_lres,
               names(clusters)[clusters == good.clusters[1]], names(clusters)[clusters == good.clusters[2]], repl.fig.dir, cols=MODEL.CLUSTER.COLORS[good.clusters])

  cells.to.calc.scores = cells.in.ab.form
  cells.to.cluster = union(cells.to.calc.scores[decay.metrics[cells.to.calc.scores, 'group'] == 3], colnames(model.params$total.cov))

  bin_cell_ab_trans = get.or.create(BIN.CELL.AB.TRANS.PATH, function() stop())
  ab.form.trans.path = file.path(repl.data.dir, 'ab_form_emb_proper_trans')
  ab.form.trans = get.or.create(ab.form.trans.path, function() total.ab.format(cells.in.ab.form, bin_cell_ab_trans))


  # classify more cells
  all.a.scores = get.bin.a.scores(names(clusters), rownames(model.params$total.cov), ab.form.trans, TRUE)
  scores.by.sampling = lapply(de.bins, function(cls.bins) {get.scores.by.sampling(cells.to.calc.scores, all.a.scores, cls.bins, nbreaks=6, ab.form.trans)})
  scores.norm.all.cells = lapply(scores.by.sampling, function(cls.scores.df) (cls.scores.df[1,] - colMeans(cls.scores.df[2:nrow(cls.scores.df),],)))
  scores.norm = lapply(scores.norm.all.cells, function(x) x[cells.to.cluster])
  scores.unnorm.all.cells = lapply(scores.by.sampling, function(cls.scores.df) (cls.scores.df[1,]))
  scores.unnorm = lapply(scores.unnorm.all.cells, function(x) x[cells.to.cluster])
  ab.cluster.ret = get.ab.cluster.assignments(model.params, scores.norm, good.clusters)
  plot.a.scores.cell.cycle(scores.norm.all.cells, scores.unnorm.all.cells, model.params, sch_decay_metrics, ab.cluster.ret[3:4], repl.fig.dir)

  ab.model.params = ab.cluster.ret[[1]]
  all.cluster.assignment = ab.cluster.ret[[2]]
  
  new.cells = setdiff(names(scores.norm[[1]]), names(clusters))
  new.assigned.cells = setdiff(names(all.cluster.assignment), names(clusters))
  new.cluster.assignment = rep(0, length(new.cells))
  names(new.cluster.assignment) = new.cells
  new.cluster.assignment[new.assigned.cells] = all.cluster.assignment[new.assigned.cells]
  ab.model.params$cols = MODEL.CLUSTER.COLORS[good.clusters]
  
  plot.ab.scores(clusters, scores.norm, ab.cluster.ret[3:4], 'ab_scores_orig_model.png', fig.dir=repl.fig.dir, MODEL.CLUSTER.COLORS)
  plot.ab.scores(new.cluster.assignment + 1, scores.norm, ab.cluster.ret[3:4], 'ab_scores.png', fig.dir=repl.fig.dir, c('grey', ab.model.params$cols))
  
  #save(ab.model.params, file=MODEL.PATH)
  #save(all.cluster.assignment, file=CLUSTERING.PATH)

  e14.analysis()

  full.cluster.assignment = all.cluster.assignment
  full.cluster.assignment[repl.erys] = 3
  full.cluster.assignment[repl.esc] = 4

  # now shaman plots
  ab.form.hres.path = file.path(repl.data.dir, 'ab_form_hres_emb_proper')
  ab.form.hres = get.or.create(ab.form.hres.path, function() total.ab.format(names(full.cluster.assignment), bin_cell_ab_hres))
  chrom.ab.hres = get.or.create(file.path(repl.data.dir, 'chrom_ab_hres_emb_proper'),
                  function() compute.chrom.ab(full.cluster.assignment, ab.form.hres, min.cov.discard=40, min.cov.per.clust=40, min.cov.cells=names(full.cluster.assignment)))

  set.misha(SHAMAN.MISHA.PATH)
  track.names = paste0(EMB.CLS.ANALYSIS.TRACK.NAMES, '_equal_size')
  score.track.names = paste0(track.names[1:2], '_matshuff_scores')
  all.tracks.conts = lapply(score.track.names, function(track.name) gextract(track.name, gintervals.2d.all(), colnames='score'))

  embryo.proper.genes = c('Sox2', 'Twist1', 'Igf2', 'Crabp2')
  gene.md = get.gene.metadata.jumps()
  gene.md = gene.md[!duplicated(gene.md$geneSymbol),]
  rownames(gene.md) = gene.md$geneSymbol
  all.interesting.bins = genes.to.intervs(embryo.proper.genes, gene.md, F, delta=c(1e6, 1e6, 2e6, 2e6))
  all.interesting.bins$geneSymbol = as.character(all.interesting.bins$geneSymbol)
  all.interesting.bins$clust = 3

  plot.dir.name = file.path(repl.fig.dir, 'hic_matrices')
  dir.create(plot.dir.name, showWarnings=F)

  for (i in 1:nrow(all.interesting.bins)) {
    interv = all.interesting.bins[i, ]
    forced.interv = gintervals.force_range(data.frame(chrom=interv$chrom, start=interv$start, end=interv$end))
    stopifnot(nrow(interv) == nrow(forced.interv))
    forced.interv$orig_start = interv$tss
    forced.interv$tss = interv$tss
    forced.interv$clust = interv$clust
    forced.interv$geneSymbol = interv$geneSymbol
    interv = forced.interv
    interv$chrom = as.character(interv$chrom)

    plot.file.name = file.path(plot.dir.name, paste0(interv$geneSymbol, '_', interv$chrom, '_', 
                       trimws(format(interv$start, scientific=F)), '.png'))
    if (interv$geneSymbol == 'Twist1') {
      plot.whole.chrom(chrom.ab.hres, chrom.ab.hres, chrom.ab.hres, all.tracks.conts, as.numeric(substr(interv$chrom, 4, 6)), plot.file.name, coord.range=c(interv$start, interv$end), mark.coord=interv$orig_start, cls.col=MODEL.CLUSTER.COLORS[good.clusters], smooth.cov=F, height.vec=c(3, 3, 1.5), width=6, sel.attr.plots=c(2), show.coord.text=F, smooth.ab=T, rotate=T, coord.on.bottom=c(interv$orig_start, 34906190))
    } else {
      plot.whole.chrom(chrom.ab.hres, chrom.ab.hres, chrom.ab.hres, all.tracks.conts, as.numeric(substr(interv$chrom, 4, 6)), plot.file.name, coord.range=c(interv$start, interv$end), mark.coord=interv$orig_start, cls.col=MODEL.CLUSTER.COLORS[good.clusters], smooth.cov=F, height.vec=c(3, 3, 1.5), width=6, sel.attr.plots=c(2), show.coord.text=F, smooth.ab=T, rotate=T, coord.on.bottom=c(interv$orig_start)) 
    }
  }

}

plot.ab.scores <- function(clusters, ab.scores, lines=list(), fig.name, fig.dir=get.fig.dir(), cols=1:length(unique(clusters))) {
	png(file.path(fig.dir, fig.name))
	plot(ab.scores[[1]][names(clusters)], ab.scores[[2]][names(clusters)], col=cols[clusters], pch=19, cex=2)
        for (line in lines) {
          abline(b=line[[1]], a=line[[2]], col='grey', lwd=4)
        }
	dev.off()	
}

get.bin.a.scores <- function(cells, bins, ab.form, norm.scores=T) {
  if (norm.scores) {
    cls.score = (colSums(ab.form$atab[cells,]) / colSums((ab.form$atab[cells,] + ab.form$btab[cells,]))) / 
                (colSums(ab.form$atab) / colSums((ab.form$atab + ab.form$btab)))
  } else {
    cls.score = colSums(ab.form$atab[cells,]) / colSums((ab.form$atab[cells,] + ab.form$btab[cells,]))
  }
  return(cls.score[bins])
}

get.cell.a.scores <- function(cells, bins, ab.form) {
  cls.score = (rowSums(ab.form$atab[,bins]) / rowSums((ab.form$atab[,bins] + ab.form$btab[,bins]))) / 
                (rowSums(ab.form$atab) / rowSums((ab.form$atab + ab.form$btab)))
  return(cls.score[cells])
}

sample.bins.with.scores <- function(all.a.scores, selected.bins, nbreaks) {
  a.scores = all.a.scores[selected.bins]
  breaks = quantile(all.a.scores, (0:nbreaks) / nbreaks)
  # hack to avoid problems with the minimum
  breaks[1] = breaks[1] - 0.001
  bin.sizes = table(cut(a.scores, breaks))
  binned.bins = cut(all.a.scores, breaks)
  names(binned.bins) = names(all.a.scores)
  
  all.sampled.bins = unlist(lapply(1:length(bin.sizes), function(i) {
    bins.to.sample.from = names(binned.bins)[binned.bins == names(bin.sizes)[i]]
    sampled.bins = sample(bins.to.sample.from, bin.sizes[i])
    return(sampled.bins)
  }))
  return(all.sampled.bins)
}

get.scores.by.sampling <- function(cells, all.a.scores, selected.bins, nbreaks, ab.form, nsamples=100) {
  all.sampled.bins = mclapply(1:nsamples, function(i) {sample.bins.with.scores(all.a.scores, selected.bins, nbreaks)}, mc.cores=25)
  all.scores = get.cluster.score.per.cell(cells, c(list(selected.bins), all.sampled.bins), ab.form, ncores=25)
  scores.df = do.call(rbind, all.scores)
  return(scores.df)
}

get.cluster.score.per.cell <- function(cells, de.bins, ab.form, ncores=1) {
  return(mclapply(1:length(de.bins), function(i) {
    cls.bins = de.bins[[i]]
    return(get.cell.a.scores(cells, cls.bins, ab.form))
  }, mc.cores=ncores))
}

get.ab.cluster.assignments <- function(model.params, cluster.scores, cluster.indices) {

  cluster.scores1 = cluster.scores[[1]][colnames(model.params$total.cov)]
  cluster.scores2 = cluster.scores[[2]][colnames(model.params$total.cov)]

  df1 = data.frame(score1=cluster.scores1, score2=cluster.scores2, cls=as.factor(model.to.clusters(model.params) == cluster.indices[1]))
  df2 = data.frame(score1=cluster.scores1, score2=cluster.scores2, cls=as.factor(model.to.clusters(model.params) == cluster.indices[2]))
  df3 = data.frame(score1=cluster.scores[[1]], score2=cluster.scores[[2]])

  is_cls1 = as.factor(model.to.clusters(model.params) == cluster.indices[1])
  is_cls2 = as.factor(model.to.clusters(model.params) == cluster.indices[2])
  fit1 = train(cls ~., data=df1, method='svmLinear')
  fit2 = train(cls ~., data=df2, method='svmLinear')
  preds1 = as.logical(predict(fit1, newdata=df3, type='raw'))
  preds2 = as.logical(predict(fit2, newdata=df3, type='raw'))
  line1 = linear.svm.to.line(fit1)
  line2 = linear.svm.to.line(fit2)

  new.cluster.assignment = rep(0, nrow(df3))
  names(new.cluster.assignment) = names(cluster.scores[[1]])
  new.cluster.assignment[preds1] = new.cluster.assignment[preds1] + 1
  new.cluster.assignment[preds2] = new.cluster.assignment[preds2] + 2

  stopifnot(sum(new.cluster.assignment == 3) <= 5)
  new.cluster.assignment[new.cluster.assignment == 3] = 1

  all.cluster.assignment = new.cluster.assignment[new.cluster.assignment != 0]
  new.cluster.assignment = new.cluster.assignment[colnames(model.params$total.cov)]

  new.model.params = model.params
  new.model.params$total.cov = model.params$total.cov[, new.cluster.assignment %in% c(1, 2)]
  new.model.params$var.mult = model.params$var.mult[, new.cluster.assignment %in% c(1, 2)]
  new.model.params$s.scores = model.params$s.scores[new.cluster.assignment %in% c(1, 2)]
  
  new.bin.clusters = matrix(NA, ncol=2, nrow=nrow(model.params$total.cov))
  new.bin.clusters[, 1] = model.params$bin.clusters[, cluster.indices[1]]
  new.bin.clusters[, 2] = model.params$bin.clusters[, cluster.indices[2]]
  new.model.params$bin.clusters = new.bin.clusters
  
  new.e.z = matrix(0, ncol=2, nrow=ncol(model.params$total.cov))
  new.e.z[new.cluster.assignment == 1, 1] = 1
  new.e.z[new.cluster.assignment == 2, 2] = 1
  new.e.z = new.e.z[new.cluster.assignment %in% c(1, 2), ]
  new.model.params$e.z = new.e.z
  new.model.params$mixture.fractions = colSums(new.e.z) / sum(new.e.z)
  return(list(new.model.params, all.cluster.assignment, line1=line1, line2=line2, new.cluster.assignment=new.cluster.assignment))

}

linear.svm.to.line <- function(fit.model) {
  # A very horrible patch to find the line equation of a fitted model.
  # Did not find how the original coefficients can be extracted from the model.
  search.range = seq(-1, 1, 0.001)
  val1 = NA
  val2 = NA
  for (i in 1:(length(search.range) - 1)) {
    df1 = data.frame(score1=c(0, 0), score2=c(search.range[i], search.range[i+1]))
    df2 = data.frame(score1=c(-0.1, -0.1), score2=c(search.range[i], search.range[i+1]))
    preds1 = as.logical(predict(fit.model, newdata=df1, type='raw'))
    preds2 = as.logical(predict(fit.model, newdata=df2, type='raw'))
    if (preds1[1] != preds1[2]) {
      val1 = search.range[i]
    }
    if (preds2[1] != preds2[2]) {
      val2 = search.range[i]
    }
  }
  stopifnot(!is.na(val1))
  stopifnot(!is.na(val2))
  return(list(slope=10*(val1 - val2), intercept=val1))
}

plot.a.scores.cell.cycle <- function(scores.norm.all.cells, scores.unnorm.all.cells, model.params, decay.metrics, lines, fig.dir) {
  all.cells = names(scores.norm.all.cells[[1]])
  model.cells = names(model.to.clusters(model.params))
  s.scores = model.params$s.scores
  names(s.scores) = model.cells
  cols = rep(NA, length(all.cells))
  names(cols) = all.cells
  set2.cols = brewer.pal(6, 'Set2')
  col.options = CELL.CYCLE.COLS
  col.options[3] = 'lightgrey'
  col.options[6] = 'grey'

  cell.cycle.phase = sch_decay_metrics[all.cells, 'group']
  cols = col.options[cell.cycle.phase]
  names(cols) = all.cells
  cols[names(cols) %in% model.cells] = col.options[6]
  is.g1.cells = cell.cycle.phase == 2
  is.s.cells = cell.cycle.phase == 3
  is.g2.cells = cell.cycle.phase == 4
  png(file.path(fig.dir, 'bin_norm.png'), width=1600, height=1000)
  par(mfrow=c(2, 3))
  for (i in 1:2) {
    cur.scores = list(scores.norm.all.cells, scores.unnorm.all.cells)[[i]]
    xrange = range(cur.scores[[1]])
    yrange = range(cur.scores[[2]])

    for (j in 1:3) {
      cells.to.draw = list(is.g1.cells, is.g2.cells, !is.s.cells)[[j]]
      plot(cur.scores[[1]][is.s.cells], cur.scores[[2]][is.s.cells], pch=19, xlim=xrange, ylim=yrange, col=cols[is.s.cells], labels=j == 1, ylab='')
      cex = ifelse(j==3, 1, 1.5)
      points(cur.scores[[1]][cells.to.draw], cur.scores[[2]][cells.to.draw], pch=19, xlim=xrange, ylim=yrange, col=cols[cells.to.draw])
      if (i == 1) {
        for (line in lines) {
          abline(b=line[[1]], a=line[[2]], col='grey', lwd=2)
        }
      }
    }
  }
  dev.off()

}

e14.analysis <- function() {
  # after loading all
  gsetroot(SCHIC.MISHA.PATH)
  strict.ab = early_late_bins
  bins.to.df <- function(bins) {
    splitted = strsplit(bins, '_')
    chrom = sapply(splitted, function(x) x[[1]])
    start = as.numeric(sapply(splitted, function(x) x[[2]]))
    end = pmin(start + 2e5, get.chrom.sizes()[as.numeric(substr(chrom, 4, 6))])
    return(data.frame(chrom=chrom, start=start, end=end))
  }
  strict.a.df = bins.to.df(strict.ab[[1]])
  strict.b.df = bins.to.df(strict.ab[[2]])
  strict.a.df$ab_tor = 'a'
  strict.b.df$ab_tor = 'b'
  intervs_ab = rbind(strict.a.df, strict.b.df)
  
  
  bonev.tname = "hic.SC.bonev_npc.NPC_Bonev"
  hsc.tname = "hic.SC.chen_hsc.HSC_HSC"

  repl.data.dir = get.data.file.dir()
  e14.data.dir = file.path(repl.data.dir, 'e14')
  dir.create(e14.data.dir, showWarnings=F)

  ab1.path = file.path(e14.data.dir, 'hsc_bin_ab_ret_5e5')
  ab2.path = file.path(e14.data.dir, 'hsc_bin_ab_ret2_5e5')
  ab3.path = file.path(e14.data.dir, 'hsc_bin_ab_ret3_5e5')
  ab4.path = file.path(e14.data.dir, 'hsc_bin_ab_ret4_5e5')
  if (!file.exists(ab1.path)) {
    bin_ab_ret = schic_init_cell_bin_ab_direct(bonev.tname, intervs_ab, cov_binsize=4e4, min_dist=5e5)
    save(bin_ab_ret,  file=ab1.path)
  }
  if (!file.exists(ab2.path)) {
    bin_ab_ret2 = schic_init_cell_bin_ab_direct(hsc.tname, intervs_ab, cov_binsize=4e4, min_dist=5e5)
    save(bin_ab_ret2, file=ab2.path)
  }
  if (!file.exists(ab3.path)) {
    set.misha(SHAMAN.MISHA.PATH)
    bin_ab_ret3 = schic_init_cell_bin_ab_direct("cluster_ds_paper_1_contacts", intervs_ab, cov_binsize=4e4, min_dist=5e5)
    save(bin_ab_ret3, file=ab3.path)
  }
  if (!file.exists(ab4.path)) {
    set.misha(SHAMAN.MISHA.PATH)
    bin_ab_ret4 = schic_init_cell_bin_ab_direct("cluster_ds_paper_2_contacts", intervs_ab, cov_binsize=4e4, min_dist=5e5)
    save(bin_ab_ret4, file=ab4.path)
  }
  gsetroot(SCHIC.MISHA.PATH)

  load(ab1.path, v=T)
  load(ab2.path, v=T)
  load(ab3.path, v=T)
  load(ab4.path, v=T)


  used.bins = Reduce(intersect, lapply(list(bin_ab_ret, bin_ab_ret2, bin_ab_ret3, bin_ab_ret4), function(x) x$f1))
  bin_ab_ret_concat = do.call(rbind, lapply(list(bin_ab_ret, bin_ab_ret2, bin_ab_ret3, bin_ab_ret4), function(x) x[x$f1 %in% used.bins,]))
  ab.form.concat = total.ab.format(sapply(list(bin_ab_ret, bin_ab_ret2, bin_ab_ret3, bin_ab_ret4), function(x) unique(x$cell)), bin_ab_ret_concat)

  cov.per.bin = apply(ab.form.concat$atab + ab.form.concat$btab, 2, min)
  new.atab = c()
  new.btab = c()
  # be wary of parallelization because of scm_downsamp
  for (i in 1:4) {
    print(i)
    cur.tabs = data.frame(cur.atab = unlist(ab.form.concat$atab[i,,drop=T]), 
                          cur.btab = unlist(ab.form.concat$btab[i,,drop=T]))
    # next line is a terrible hack - I do downsampling twice (rep(j, 5)) because otherwise the function fails
    # I do 5 because with 2 there's a chance that the output matrix will be of format "dtrMatrix"
    #cur.tabs.ds = do.call(rbind, lapply(1:nrow(cur.tabs), function(j) {
    #  if (j %% 1e4 == 1) print(j)
    #  scm_downsamp(t(cur.tabs[rep(j, 5),,drop=F]), cov.per.bin[j])[,1]
    #}))
    cur.tabs.ds = do.call(rbind, lapply(1:nrow(cur.tabs), function(j) {
      if (j %% 1e3 == 1) print(j)
      table(factor(sample(c(rep('a', cur.tabs[j, 1]), rep('b', cur.tabs[j, 2])), cov.per.bin[j], replace=F), levels=c('a', 'b')))
    }))
    stopifnot(all(rowSums(cur.tabs.ds) == cov.per.bin))
    new.atab = cbind(new.atab, cur.tabs.ds[,1])
    new.btab = cbind(new.btab, cur.tabs.ds[,2])
  }

  rownames(new.atab) = colnames(ab.form.concat$atab)
  rownames(new.btab) = colnames(ab.form.concat$btab)
  colnames(new.atab) = rownames(ab.form.concat$atab)
  colnames(new.btab) = rownames(ab.form.concat$btab)
  ab.form.concat$prev.atab = ab.form.concat$atab
  ab.form.concat$prev.btab = ab.form.concat$btab
  ab.form.concat$atab = t(new.atab)
  ab.form.concat$btab = t(new.btab)

  save(ab.form.concat, file=file.path(e14.data.dir, 'e14_ab_form_concat'))
  load(file.path(e14.data.dir, 'e14_ab_form_concat'), v=T)


  clusters = 1:4
  cnames = rownames(ab.form.concat$atab)
  names(clusters) = cnames
  chrom.ab = compute.chrom.ab(clusters, ab.form.concat, min.cov.per.clust=80, min.cov.discard=0)
  chrom.ab.fil = chrom.ab[[1]][rowSums(is.na(chrom.ab[[1]])) == 0,]
  ecto.bins = rownames(chrom.ab.fil)[chrom.ab.fil[,3] - chrom.ab.fil[,4] > 0.15]
  meso.bins = rownames(chrom.ab.fil)[chrom.ab.fil[,4] - chrom.ab.fil[,3] > 0.15]
  bin.cols = ifelse(rownames(chrom.ab.fil) %in% ecto.bins, MODEL.CLUSTER.COLORS[1], ifelse(rownames(chrom.ab.fil) %in% meso.bins, MODEL.CLUSTER.COLORS[3], 1))

  fig.dir = file.path(get.fig.dir(), 'e14_analysis')
  dir.create(fig.dir, showWarnings=F)


  sel.bins = bin.cols != 1
  png(file.path(fig.dir, 'diff_boxplot.png'))
  boxplot(split(chrom.ab.fil[sel.bins, 1] - chrom.ab.fil[sel.bins, 2], bin.cols[sel.bins]), col=MODEL.CLUSTER.COLORS[-2])
  dev.off()

  diff1 = chrom.ab.fil[,2] - chrom.ab.fil[,1]
  diff2 = chrom.ab.fil[,4] - chrom.ab.fil[,3]
  abs.diff1 = abs(diff1)
  abs.diff2 = abs(diff2)

  png(file.path(fig.dir, 'abs_diff_hist.png'))
  hist(abs.diff1, freq=F, col=rgb(1, 0, 0, 1/4), ylim=c(0, 25), breaks=seq(0, 0.6, by=0.025), xlim=c(0, 0.4))
  hist(abs.diff2, freq=F, col=rgb(0, 0, 1, 1/4), add=T, breaks=seq(0, 0.6, by=0.025))
  dev.off()

  png(file.path(fig.dir, 'abs_diff_hist_legend.png'))
  plot(1, 1, bg=rgb(1, 0, 0, 1/4), xlim=c(0, 3), ylim=c(0, 3), pch=21, cex=4)
  points(2, 2, bg=rgb(0, 0, 1, 1/4), pch=21, cex=4)
  dev.off()

  png(file.path(fig.dir, 'e14_scores.png'))
  plot(chrom.ab.fil[bin.cols == 1, 1], chrom.ab.fil[bin.cols == 1, 2], bg=1, pch=21, xlim=c(0, 1), ylim=c(0, 1))
  points(chrom.ab.fil[bin.cols != 1, 1], chrom.ab.fil[bin.cols != 1, 2], bg=bin.cols[bin.cols != 1], pch=21)
  dev.off()

  png(file.path(fig.dir, 'e9_scores.png'))
  plot(chrom.ab.fil[bin.cols == 1, 3], chrom.ab.fil[bin.cols == 1, 4], bg=1, pch=21, xlim=c(0, 1), ylim=c(0, 1))
  points(chrom.ab.fil[bin.cols != 1, 3], chrom.ab.fil[bin.cols != 1, 4], bg=bin.cols[bin.cols != 1], pch=21)
  dev.off()

  print('A-score diffs in E14.5:')
  print(mean(abs.diff1 > 0.2))
  print('A-score diffs in E9.5:')
  print(mean(abs.diff2 > 0.2))

}

get.repl.esc <- function(esc.decay.metrics) {
  repl.cells = esc.decay.metrics$cell[esc.decay.metrics[, 'early_f'] > 0.52 & esc.decay.metrics[, 'group'] == 3]
  return(repl.cells)
}

get.repl.erys <- function(decay.metrics) {
  erys = get.erys()
  repl.erys = erys[decay.metrics[erys, 'early_f'] > 0.52]
  return(repl.erys)
}
