library(zoo)
library(umap)
library(sm)
library(pheatmap)
library(caret)
library(reshape2)

norm.cells.and.correlate <- function(total.cov, cells.ord, feat.min, feat.max, ratio.min=2, rollmean.k=20, norm.rollmean=T) {
  total.cov.ord = total.cov[, cells.ord]
  feat.means = rowMeans(total.cov.ord) 
  feat.vars = apply(total.cov.ord, 1, var)

  selected.feats = feat.means > feat.min & feat.means < feat.max & feat.vars / feat.means > ratio.min
  total.cov.ord = total.cov.ord[selected.feats,]
  if (norm.rollmean) {
    total.cov.norm = total.cov.ord - t(apply(total.cov.ord, 1, rollmean, rollmean.k, fill='extend'))
  } else {
    total.cov.norm = total.cov.ord - rowMeans(total.cov.ord) #t(apply(total.cov.ord, 1, mean))
  }
  
  cell.cor = tgs_cor(total.cov.norm)
  hc = hclust(dist(cell.cor), "ward.D2") 
  return(list(total.cov.norm=total.cov.norm, cell.cor=cell.cor, hc=hc, feats=selected.feats))
}

init.clustering <- function(total.cov, cells.ord, nclust, feat.min, feat.max, ratio.min=2, rollmean.k=20) {

  hc = norm.cells.and.correlate(total.cov, cells.ord, feat.min, feat.max, ratio.min=ratio.min, rollmean.k=rollmean.k)$hc
  clust = cutree(hc, k=nclust)
  clust = clust[colnames(total.cov)]
  return(clust)
}

cluster.embryo.cells <- function(nclust=3, repl.duration=11, num.bin.clusters=12, lambda=40, selected.cells=NULL, use.init.func=T, cell.min.cov=8, bin.min.cov=8) {
  ret = get.embryo.cells.model.params(nclust, repl.duration, num.bin.clusters, lambda, selected.cells, use.init.func, cell.min.cov=cell.min.cov, bin.min.cov=bin.min.cov)
  model.params = ret$model.params
  unnorm.total.cov = ret$unnorm.total.cov
  total.cov = ret$total.cov
  init.model.params = model.params
  
  model.ret = run.clustering.iteration.ds(total.cov, model.params, FALSE)
  ret.model.params = model.ret[[1]]
  total.cov = ret.model.params$total.cov
  return(list(model.params=ret.model.params, init.model.params=init.model.params, norm.total.cov=total.cov, unnorm.total.cov=unnorm.total.cov))
}

embryo.init.params = list(feat.min=6e-5, feat.max=1.3e-4, ratio.min=1.75e-5, rollmean.k=20)

embryo.init.clusters.func = function(total.cov, cells.ord, nclust) {
  do.call(init.clustering, c(list(total.cov=total.cov, cells.ord=cells.ord, nclust=nclust), embryo.init.params))
}

get.embryo.cells.model.params <- function(nclust=3, repl.duration=4, num.bin.clusters=6, lambda=0, 
                                          selected.cells=NULL, use.init.func=T, cell.min.cov=8, bin.min.cov=8) {
  total.cov = bin_cell_cov
  decay.metrics = sch_decay_metrics
  if (is.null(selected.cells)) {
    selected.cells = colnames(total.cov)
  }
  
  erys = get.erys()

  g1.cells = rownames(decay.metrics)[decay.metrics$group == 2]
  g1.cells = setdiff(g1.cells, erys)
  g1.cells = intersect(g1.cells, colnames(total.cov))
  
  if (use.init.func) {
    init.clusters.func = embryo.init.clusters.func
  } else {
    init.clusters.func = NULL
  }

  sel.cells = filter.embryo.cells(decay.metrics, erys)
  sel.cells = intersect(sel.cells, colnames(total.cov))
  
  
  sel.cells = intersect(sel.cells, selected.cells)

  ret = get.model.params(total.cov, sel.cells, decay.metrics, g1.cells, lambda=lambda,
                init.clusters.func=init.clusters.func, nclust=nclust, repl.duration=repl.duration, 
		num.bin.clusters=num.bin.clusters, cell.min.cov=cell.min.cov, bin.min.cov=bin.min.cov)
  return(ret)
  
}

get.model.params <- function(total.cov, cells, decay.metrics, g1.cells, init.clusters = NULL, init.clusters.func = NULL,
                          nclust = 4, cell.min.cov = 8, bin.min.cov = 8, 
                          repl.duration = 12, num.bin.clusters = 20, lambda = 0, 
                          bin.min.cov.g1 = bin.min.cov, order.cells.by = 'near_f') {
  set.seed(42)
  stopifnot(all(cells %in% colnames(total.cov)))
  if (!is.null(init.clusters)) {
    stopifnot(all(cells %in% names(init.clusters)))  
  }
  stopifnot(all(cells %in% rownames(decay.metrics)))
  stopifnot(all(g1.cells %in% colnames(total.cov)))
  sel.cells = cells[colMeans(total.cov[,cells]) >= cell.min.cov]
  ncells = length(sel.cells)

  s.scores = max.min.rescale(1:length(sel.cells), 1.2, 1.8)
  cells.ord = sel.cells[order(decay.metrics[sel.cells, order.cells.by])]
  names(s.scores) = cells.ord
  s.scores = s.scores[sel.cells]

  # bin filtering
  for (chrom in c('chrX', 'chrY', 'chrM')) {
    chrom.bins = grepl(chrom, rownames(total.cov))
    total.cov = total.cov[!chrom.bins,]
  }

  covered.locs = rownames(total.cov)[rowMeans(total.cov[, sel.cells]) > bin.min.cov]
  g1.covered.locs = rownames(total.cov)[rowMeans(total.cov[, g1.cells]) > bin.min.cov.g1]
  selected.locs = intersect(intersect(rownames(total.cov), covered.locs), g1.covered.locs)
  total.cov = total.cov[selected.locs, ]

  loc.prob.vector = rowSums(total.cov[, g1.cells])
  loc.prob.vector = loc.prob.vector / sum(loc.prob.vector)

  total.cov = total.cov[selected.locs, sel.cells]
  unnorm.total.cov = total.cov

  

  cell.total.cov = colSums(total.cov)
  total.cov = t(t(total.cov) / cell.total.cov) 

  # cluster initialization
  if (is.null(init.clusters)) {
    if (is.null(init.clusters.func)) {
      init.clusters = sample(1:nclust, ncells, replace=T)
    } else {
      init.clusters = init.clusters.func(total.cov, cells.ord, nclust)
    }
    names(init.clusters) = sel.cells  
    
  } else {
    init.clusters = init.clusters[sel.cells]
    nclust = length(unique(init.clusters))
  }
  names(init.clusters) = sel.cells
  
  e.z = get.init.e.z(sel.cells, nclust, init.clusters)
  mixture.fractions = colMeans(e.z)

  # create model parameters
  model.params = list()
  model.params$mixture.fractions = mixture.fractions
  model.params$max.copy.num = 2
  model.params$repl.duration = repl.duration
  model.params$num.bin.clusters = num.bin.clusters
  model.params$e.z = e.z
  model.params$s.scores = s.scores
  model.params$var.mult = t(1 / (cell.total.cov %*% t(rep(1, nrow(total.cov)))))
  model.params$bin.probs = loc.prob.vector
  model.params$lambda = lambda
  return(list(model.params=model.params, total.cov=total.cov, unnorm.total.cov=unnorm.total.cov))
  
}

get.init.e.z <- function(sel.cells, nclust, init.clusters) {
  ncells = length(sel.cells)
  if (nclust > 2) {
    e.z = matrix(0.5 / (nclust - 1), nrow=ncells, ncol=nclust)  
  } else {
    e.z = matrix(1 / 3, nrow=ncells, ncol=nclust)  
  }
  
  rownames(e.z) = sel.cells
  colnames(e.z) = names(table(init.clusters))

  main.clust.prob = ifelse(nclust > 2, 0.5, 2 / 3)
  for (i in seq_along(sel.cells)) {
    e.z[i, as.character(init.clusters[i])] = main.clust.prob
  }
  return(e.z)
}

filter.embryo.cells <- function(decay.metrics, erys=NULL) {
  if (is.null(erys)) {
    erys = get.erys()
  }
  sel.cells = rownames(decay.metrics[between(decay.metrics$near_f, 0.72, 0.84),])
  #sel.cells = intersect(sel.cells, rownames(decay.metrics[decay.metrics$early_f < 0.6,]))
  sel.cells = intersect(sel.cells, rownames(decay.metrics[decay.metrics$early_f > 0.54,]))
  sel.cells = intersect(sel.cells, rownames(decay.metrics)[decay.metrics$group %in% c(3, 4)])
  sel.cells = sel.cells[!(sel.cells %in% erys)]
  return(sel.cells)
}

predict.cov <- function(model.params) {
  pred.total.cov = model.params$total.cov
  pred.total.cov[,] = NA
  ncells = ncol(pred.total.cov)

  for (i in 1:ncells) {
    bin.exps = model.params$bin.probs * get.bin.expectation(model.params$s.scores[i], model.params$bin.clusters, 
                                   model.params$max.copy.num, model.params$num.bin.clusters, model.params$repl.duration)
    bin.exps = t(t(bin.exps) / colSums(bin.exps))
    pred.total.cov[, i] = rowSums(t(model.params$e.z[i,] * t(bin.exps)))
  }

  cors = sapply(1:ncells, function(i) cor(model.params$total.cov[, i], pred.total.cov[, i]))
  return(list(pred=pred.total.cov, cors=cors))
}


test.model <- function() {
  decay.metrics = get.common.decay.metrics()
  fig.dir = get.fig.dir()
  fig.dir.test.model = file.path(fig.dir, 'test_model')
  dir.create(fig.dir.test.model, showWarnings=F)

  embryo.cls.path = file.path(repl.data.dir, 'embryo_cls')
  embryo.cls.ret = get.or.create(embryo.cls.path, cluster.embryo.cells)
  model.params = embryo.cls.ret$model.params
  unnorm.total.cov = embryo.cls.ret$unnorm.total.cov

  # simulations
  sim.ret.path = file.path(repl.data.dir, 'sim_results')
  sim.ret = get.or.create(sim.ret.path, function(x) run.simulation(ncells=240, nbins=9600, nclusters=5, repl.duration=1, num.bin.clusters=3))
  sim.model.params = sim.ret$clust.ret[[1]]
  pred.s.scores = sim.model.params$s.scores
  true.s.scores = sim.ret$s.scores
  png(file.path(fig.dir.test.model, 'sim_pred_s_scores.png'))
  plot(pred.s.scores, true.s.scores, cex=1.5, pch=19)
  dev.off()

  # parameter tuning
  cross.val.error.path = file.path(repl.data.dir, 'cross_val_error')
  repli.cross.val.errors = get.or.create(cross.val.error.path, function() get.replication.cross.val.errors(max.num.bin.clusters=6))
  attempted.num.bin.clusters = c(3:15, seq(20, 45, 5))
  cross.val.error.specific.path = file.path(repl.data.dir, 'cross_val_error_specific')
  repli.cross.val.errors.specific = get.or.create(cross.val.error.specific.path, function() sapply(attempted.num.bin.clusters, function(i) get.replication.cross.val.errors.specific(i, i - 1)))
  names(repli.cross.val.errors.specific) = attempted.num.bin.clusters
  plot.repli.cross.val.errors(repli.cross.val.errors, repli.cross.val.errors.specific, fig.dir.test.model)
  
  cross.val.error.second.path = file.path(repl.data.dir, 'cross_val_error_lambda')
  clust.options = 2:5
  lambda.options = seq(0, 60, 20)
  lambda.cross.val.errors = get.or.create(cross.val.error.second.path, function() get.clustering.cross.val.errors(model.params$num.bin.clusters, model.params$repl.duration, clust.options, lambda.options))
  rownames(lambda.cross.val.errors) = clust.options
  colnames(lambda.cross.val.errors) = lambda.options
  plot.lambda.cross.val.errors(lambda.cross.val.errors, fig.dir.test.model)
  cross.val.error.third.path = file.path(repl.data.dir, 'cross_val_error_num_cls')
  clust.cross.val.errors2 = get.or.create(cross.val.error.third.path, function() get.clustering.cross.val.errors(model.params$num.bin.clusters, model.params$repl.duration, 1:9, 0))
  plot.cv.num.cls(clust.cross.val.errors2, fig.dir.test.model)
  
  cross.val.ret = get.or.create(file.path(repl.data.dir, 'cross_val_ret'), function() {
    cross.val.embryo(model.params$num.bin.clusters, model.params$repl.duration, 
                     length(model.params$mixture.fractions), model.params$lambda)
  })
  plot.cross.val.cor(cross.val.ret, unnorm.total.cov, fig.dir=fig.dir.test.model)
  
  # stability tests
  sampling.stability.ret = get.or.create(file.path(repl.data.dir, 'stability_ret'), function() test.embryo.stability(model.params, 5))
  plot.stability(model.params, sampling.stability.ret, fig.dir.test.model, 'sampling_s_stability_with_init')
}

repl.model.analysis <- function() {
  set.misha(SCHIC.MISHA.PATH)
  decay.metrics = get.common.decay.metrics()
  fig.dir = get.fig.dir()
  repl.fig.dir = file.path(fig.dir, 'repl_model')
  dir.create(repl.fig.dir, showWarnings=F)
  repl.data.dir=get.data.file.dir()

  test.model()


  embryo.cls.path = file.path(repl.data.dir, 'embryo_cls')
  embryo.cls.ret = get.or.create(embryo.cls.path, cluster.embryo.cells)

  embryo.cls.10.path = file.path(repl.data.dir, 'embryo_cls_10')
  embryo.cls.10.ret = get.or.create(embryo.cls.10.path, function() cluster.embryo.cells(nclust=10))

  model.params = embryo.cls.ret$model.params
  model.params.10 = embryo.cls.10.ret$model.params
  unnorm.total.cov = embryo.cls.ret$unnorm.total.cov
  cov.per.cell = colSums(unnorm.total.cov)

  clusters = model.to.clusters(model.params)
  clusters10 = model.to.clusters(model.params.10)
  stopifnot(length(model.params$mixture.fractions))

  cell.norm.cor.ret = do.call(norm.cells.and.correlate, c(list(model.params$total.cov, names(clusters)[order(model.params$s.scores)]), embryo.init.params))
  cell.norm.cor.ret.unnorm = do.call(norm.cells.and.correlate, c(list(model.params$total.cov, names(clusters)[order(model.params$s.scores)]), embryo.init.params, F))
  s.scores = model.params$s.scores
  names(s.scores) = names(clusters)
  plot.cell.cor.and.2d.given.cors(clusters, cell.norm.cor.ret$cell.cor, 'cov_norm', 
                                  umap.scores=s.scores, fig.dir=repl.fig.dir, cls.col=MODEL.CLUSTER.COLORS)
  plot.cell.cor.and.2d.given.cors(clusters, cell.norm.cor.ret.unnorm$cell.cor, 'cov_unnorm', 
                                  umap.scores=s.scores, fig.dir=repl.fig.dir, cls.col=MODEL.CLUSTER.COLORS)

  plot.model.params.comparison(clusters, clusters10, fig.dir=repl.fig.dir)
  plot.cluster.density(clusters, model.params$s.scores, cov.per.cell, '3_cls', repl.fig.dir, MODEL.CLUSTER.COLORS)
  plot.cluster.density(clusters10, model.params.10$s.scores, cov.per.cell, '10_cls', repl.fig.dir, brewer.pal(10, 'Set3'))

  good.clusters = c(1, 3)
  bin.trends.fig.dir = file.path(repl.fig.dir, 'bin_trends')
  dir.create(bin.trends.fig.dir, showWarnings=F)
  # only the figure of cluster 2 is shown in the paper
  plot.bin.trends(clusters, model.params, decay.metrics, unnorm.total.cov, 
                  early.bin.cells=names(clusters)[clusters %in% good.clusters], fig.dir=bin.trends.fig.dir, cols=MODEL.CLUSTER.COLORS)

  repli_scores_emb_splitted = get.splitted.repli.scores()
  plot.rs.for.clusters(repli_scores_emb_splitted, repl.fig.dir)
  plot.cell.bin.trends(clusters, model.params, repl.fig.dir, MODEL.CLUSTER.COLORS)

  # looking for higher resolution clustering, and at our power to detect such - hierarchical clustering and sampling experiments
  hier.fig.dir = file.path(repl.fig.dir, 'hier_clustering')
  hier.analysis(hier.fig.dir)


  # subsampling experiments
  subsampled.cls.options = as.list(good.clusters)
  subsample.results.dir = file.path(repl.data.dir, 'subsample_models')
  dir.create(subsample.results.dir, showWarnings=F)
  fig.dir.subsample = file.path(repl.fig.dir, 'subsampling')
  dir.create(fig.dir.subsample, showWarnings=F)
  for (i in 1:length(subsampled.cls.options)) {
    cluster.to.sample = subsampled.cls.options[[i]]
    num.cells.options = c(100, 75, 50, 25)
    fig.name = sprintf('cls_subsampling_%s', cluster.to.sample)
    cluster.results.list = lapply(num.cells.options, function(cur.num.cells) {
      use.cur.params = num.cells.options == sum(clusters == good.clusters[i])
      cur.model.file.path = file.path(subsample.results.dir, sprintf('cls_%s_num_cells_%s', i, cur.num.cells))
      seed = cur.num.cells + i
      subsample.and.cluster(model.params, decay.metrics, num.sampled.cells=cur.num.cells, cluster.to.sample=cluster.to.sample, 
                            model.file.path=cur.model.file.path, seed=seed)
    })
    plot.subsample.and.cluster.results(cluster.results.list, model.params, clusters, 
                                       fig.dir.subsample, fig.name, MODEL.CLUSTER.COLORS)
  }


  # compare to other scHi-C methods
  # code is in method_comparison.r

  # annotation
  annot.fig.dir1 = file.path(repl.fig.dir, 'annotation1')
  annot.fig.dir2 = file.path(repl.fig.dir, 'annotation2')
  dir.create(annot.fig.dir1)
  dir.create(annot.fig.dir2)

  atlas.obj = get.or.create(ATLAS.PATH, function() error('Did not find expression atlas!'))
  e9.obj = get.or.create(E9.EXP.PATH, function() error('Did not find e9 expression data!'))

  ab_scores_emb_splitted = get.splitted.ab.scores()

  set.misha(SCHIC.MISHA.PATH)
  env = schic_init_env()
  env = schic_init_perlim_cell_groups(env)
  env = schic_init_tss_intervs(env)
  genes.to.bins = get.genes.to.bins(env, bin.size=4e4)
  rs = repli_scores_emb_splitted
  ab = ab_scores_emb_splitted
  rs = rs[apply(rs, 1, function(x) all(is.finite(x))),]
  ab = ab[apply(ab, 1, function(x) all(is.finite(x))),]
  common.bins = intersect(rownames(rs), rownames(ab))
  rs = rs[common.bins,]
  ab = ab[common.bins,]
  rs = schic.norm.rs(rs)
  ab = schic.norm.rs(ab)

  emb.gene.module.analysis(e9.obj, rs, ab, genes.to.bins, annot.fig.dir1, thresh=-16, num.gene.mods=15)
  emb.gene.module.analysis(atlas.obj, rs, ab, genes.to.bins, annot.fig.dir2, thresh=-16.7, num.gene.mods=20)

}

get.genes.to.bins <- function(env, bin.size=2e5) {
  tss = env$tss_intervs
  genes.to.bins = unlist(lapply(1:nrow(tss), function(i) {
    all.genes = strsplit(rownames(tss)[i], ';')[[1]]
    coord = as.numeric(tss$start[i])
    bin.coord = floor((coord - 1) / bin.size) * bin.size + 1
    bin.name = paste0(tss$chrom[i], '_', trimws(format(bin.coord, scientific=F)))
    ret = rep(bin.name, length(all.genes))
    names(ret) = all.genes
    return(ret)
  }))
  return(genes.to.bins)
}

schic.norm.rs = function(rs40, col.to.use=1:3)
{
	rs40 = rs40[,col.to.use]
	rs_ord = order(rowMeans(rs40,na.rm=T))
	rs40o = rs40[rs_ord,]

	all_trends = lapply(1:ncol(rs40), function(i) rollmean(rs40[rs_ord, i],200, fill='extend'))
	rsnorm = do.call(cbind, lapply(1:length(all_trends), function(i) rs40o[,i]+(all_trends[[1]]-all_trends[[i]])))

	rownames(rsnorm) = rownames(rs40)[rs_ord]
	f = !is.na(rowSums(rsnorm)) & !is.infinite(rowSums(rsnorm))
	return(rsnorm[f,])
}

emb.gene.module.analysis <- function(mc, rs, ab, genes.to.bins, fig.dir, thresh=-16, num.gene.mods=15) {
  rs = rs[apply(rs, 1, function(x) all(is.finite(x))),]
  ab = ab[apply(ab, 1, function(x) all(is.finite(x))),]
  common.bins = intersect(rownames(rs), rownames(ab))
  rs = rs[common.bins,]
  ab = ab[common.bins,]

  lfp = log2(mc@mc_fp)
  egc = mc@e_gc
  sel.genes = rownames(egc)[apply(egc, 1, max) > 2e-5 & log2(apply(egc, 1, var) / rowMeans(egc)) > thresh]
  lfp.cor = cor(t(lfp[sel.genes,]))
  hc = hclust(as.dist(1 - lfp.cor), 'ward.D2')
  gene.mods = cutree(hc, num.gene.mods)
  gene.splitted = split(names(gene.mods), gene.mods)

  all.cur.bins.length = c()
  rs.scores = do.call(rbind, lapply(gene.splitted, function(genes) {
    cur.genes = sapply(strsplit(genes, ';'), function(x) x[1])
    bins = genes.to.bins[cur.genes]
    cur.bins = intersect(bins, rownames(rs))
    cur.len = length(cur.bins)
    all.cur.bins.length[[length(all.cur.bins.length) + 1]] = cur.len
    all.cur.bins.length <<- all.cur.bins.length
    colMeans(rs[cur.bins,])
  }))

  ab.scores = do.call(rbind, lapply(gene.splitted, function(genes) {
    cur.genes = sapply(strsplit(genes, ';'), function(x) x[1])
    bins = genes.to.bins[cur.genes]
    cur.bins = intersect(bins, rownames(ab))
    colMeans(ab[cur.bins,])
  }))

  clusters = gene.mods[hc$ord]
  perm = rle(clusters)$values
  perm.orig = perm
  names(perm) = as.character(perm) 
  perm[] = 1:num.gene.mods

  rs.scores.norm = rs.scores - rowMeans(rs.scores)
  ab.scores.norm = ab.scores - rowMeans(ab.scores)
  rs.ord = rev(perm.orig)

  shades = colorRampPalette(c("darkblue", "blue", "white","red", "yellow"))(200)

  fig.path = file.path(fig.dir, 'all_gene_module_scores.png')
  x1 = pheatmap(rs.scores[rs.ord,], col=shades, cluster_rows=F, cluster_cols=F)
  x2 = pheatmap(ab.scores[rs.ord,], col=shades, cluster_rows=F, cluster_cols=F)
  x3 = pheatmap(rs.scores.norm[rs.ord,], col=shades, breaks=seq(-0.05, 0.05, length.out=201), cluster_rows=F, cluster_cols=F)
  x4 = pheatmap(ab.scores.norm[rs.ord,], col=shades, breaks=seq(-0.02, 0.02, length.out=201), cluster_rows=F, cluster_cols=F)

  arranged = grid.arrange(x1[[4]], x2[[4]], x3[[4]], x4[[4]], nrow=1, widths=c(1, 1, 1, 1), heights=4)
  ggsave(fig.path, arranged)

  png(file.path(fig.dir, 'gene_correlations.png'))
  image.plot(lfp.cor[hc$ord, hc$ord], col=shades, breaks=seq(-1, 1, length.out=201))
  dev.off()

}

plot.subsample.and.cluster.results <- function(cluster.results.list, model.params, clusters, fig.dir, fig.name, cols=1:length(unique(clusters))) {
  num.subsamples = length(cluster.results.list)
  names(model.params$s.scores) = colnames(model.params$total.cov)
  ylim = range(model.params$s.scores)

  #shades = colorRampPalette(c("darkblue", rep("blue", 4), "white", rep("red", 4), "yellow"))(200)
  shades = colorRampPalette(c("darkblue", "blue", "white", "red", "yellow"))(200)
  for (i in 1:num.subsamples) {
    new.clusters = model.to.clusters(cluster.results.list[[i]]$new.model.params)
    cur.cls.res = cluster.results.list[[i]]$cor.ret
    cor.mat = cur.cls.res$cell.cor[cur.cls.res$hc$order, cur.cls.res$hc$order]
    new.ordering = c()
    for (j in 1:max(new.clusters)) {
      new.ordering = c(new.ordering, rownames(cor.mat)[rownames(cor.mat) %in% names(new.clusters)[new.clusters == j]])
    }
    cor.mat = cor.mat[new.ordering, new.ordering]
    diag(cor.mat) = 0
    #image(cor.mat, col=shades, labels=F)
    cor.mat.lim = pmin(cor.mat, 0.3)
    cor.mat.lim = pmax(cor.mat.lim, -0.3)
    pheatmap(cor.mat.lim, col=shades, cluster_rows=F, cluster_cols=F, breaks = seq(-0.3,0.3,l=200),
             show_rownames=F, show_colnames=F, annotation_names_row=F, annotation_names_col=F, annotation_legend=F,
             filename=file.path(fig.dir, paste0(fig.name, '_heatmap_', i, '.png')),
	     annotation_col=as.data.frame(new.clusters), annotation_colors=list(new.clusters=brewer.pal(7, 'Set1')[c(1, 2, 7)]))
    png(file.path(fig.dir, paste0(fig.name, '_s_scores_', i, '.png')), width=600, height=300)
    plot(model.params$s.scores[colnames(cor.mat)], col=cols[clusters[colnames(cor.mat)]], pch=19, ylab='', ylim=ylim, cex=1.5, xaxs='i')
    dev.off()
  }
}

subsample.and.cluster <- function(model.params, decay.metrics, num.sampled.cells, cluster.to.sample, model.file.path, seed=42) {
  set.seed(seed)
  names(model.params$s.scores) = colnames(model.params$total.cov)
  clusters = model.to.clusters(model.params)
  cls.cells = names(clusters)[clusters == cluster.to.sample]
  non.cls.cells = names(clusters)[clusters != cluster.to.sample]
  sel.cells = c(non.cls.cells, sample(cls.cells, num.sampled.cells))

  new.model.params = get.or.create(model.file.path, function() cluster.embryo.cells(selected.cells=sel.cells)[[1]])
  new.clusters = model.to.clusters(new.model.params)
  ord.cells = names(new.clusters)[order(new.model.params$s.scores)]

  ret = do.call(norm.cells.and.correlate, c(list(new.model.params$total.cov[, sel.cells], ord.cells), embryo.init.params, T))
  return(list(cor.ret=ret, new.model.params=new.model.params))
}

hier.analysis <- function(hier.fig.dir) {
  dir.create(hier.fig.dir, showWarnings=F)

  gene.md = get.gene.metadata.jumps()
  gene.md$bin_formatted = sprintf('%s_%s', gene.md$chrom, gene.md$bin)
  gene.md = gene.md[!duplicated(gene.md$geneSymbol),]
  rownames(gene.md) = gene.md$geneSymbol

  decay.metrics = get.common.decay.metrics()
  repl.data.dir = get.data.file.dir()
  embryo.cls.path = file.path(repl.data.dir, 'embryo_cls')
  embryo.cls.ret = get.or.create(embryo.cls.path, cluster.embryo.cells)
  model.params = embryo.cls.ret$model.params
  clusters = model.to.clusters(model.params)
  cls1 = names(clusters)[clusters == 1]
  cls3 = names(clusters)[clusters == 3]
  cls1.ret.path = file.path(repl.data.dir, 'ecto_hier_cls')
  cls3.ret.path = file.path(repl.data.dir, 'meso_hier_cls')
  cls1.ret = get.or.create(cls1.ret.path, function() cluster.embryo.cells(selected.cells=cls1, cell.min.cov=8, bin.min.cov=8))
  cls3.ret = get.or.create(cls3.ret.path, function() cluster.embryo.cells(selected.cells=cls3, cell.min.cov=8, bin.min.cov=8))

  exp.ret = get.or.create(file.path(repl.data.dir, 'atlases'), function() error('no atlases!'))
  e9.exp = exp.ret$e9.exp
  unnorm.total.cov = bin_cell_cov

  for (i in 1:2) {
    if (i == 1) {
      tmp.clusters = rep(1, length(cls1))
      names(tmp.clusters) = cls1
      embryo.cls.ret = cls1.ret
      sel.colors = c(brain="#647a4f", neuromeso="#8EC792", crest="#C3C388", spinal="#CDE088", surface="#f7f79e")
    } else {
      tmp.clusters = rep(1, length(cls3))
      names(tmp.clusters) = cls3
      embryo.cls.ret = cls3.ret
      sel.colors = c(endothel='#ff891c', mesen='#cc7818', phary='#C9EBFB', cardio='#B51D8D', paraxial='#8DB5CE', exe='#8870ad', allan='#532C8A', intermed='#139992', somitic='#005579')
    }
    cell.norm.cor.ret = do.call(norm.cells.and.correlate, c(list(embryo.cls.ret$model.params$total.cov, names(tmp.clusters)[order(embryo.cls.ret$model.params$s.scores)]), embryo.init.params))
    cell.norm.cor.ret.unnorm = do.call(norm.cells.and.correlate, c(list(embryo.cls.ret$model.params$total.cov, names(tmp.clusters)[order(embryo.cls.ret$model.params$s.scores)]), embryo.init.params, list(F)))
    shades = colorRampPalette(c("darkblue", "blue", "white","red", "yellow"))(100)

    # first use heatmap just for clustering, without saving the figures
    tmp.ret = pheatmap(pmax(pmin(cell.norm.cor.ret$cell.cor, 0.3), -0.3), breaks=seq(-0.3, 0.3, length.out=100), color=shades, 
                                 show_rownames=F, show_colnames=F, treeheight_row = 0, treeheight_col = 0)
    tmp.ret.unnorm = pheatmap(pmax(pmin(cell.norm.cor.ret.unnorm$cell.cor, 0.3), -0.3), breaks=seq(-0.3, 0.3, length.out=100), color=shades, 
                                 show_rownames=F, show_colnames=F, treeheight_row = 0, treeheight_col = 0)

    # also plot the s scores, with and without normalization
    cur.num.clusters = c(3, 9)[i]
    s.scores = embryo.cls.ret$model.params$s.scores
    names(s.scores) = names(model.to.clusters(embryo.cls.ret$model.params))
    png(file.path(hier.fig.dir, c('ecto_s_scores.png', 'meso_s_scores.png')[i]), width=1000, height=400)
    plot(s.scores[colnames(cell.norm.cor.ret$cell.cor)[tmp.ret$tree_row$order]], xaxs='i', ylim=c(1, 2), pch=19)
    dev.off()
    png(file.path(hier.fig.dir, c('ecto_s_scores_unnorm.png', 'meso_s_scores_unnorm.png')[i]), width=1000, height=400)
    plot(s.scores[colnames(cell.norm.cor.ret.unnorm$cell.cor)[tmp.ret.unnorm$tree_row$order]], xaxs='i', ylim=c(1, 2), pch=19)
    dev.off()
    cur.cls = cutree(tmp.ret$tree_row, cur.num.clusters)
    cur.cls.df = data.frame(cur.cls=cur.cls)
    set3.cols = brewer.pal(cur.num.clusters, 'Set3')

    # now save the figures, which include the cluster annotations
    tmp.ret = pheatmap(pmax(pmin(cell.norm.cor.ret$cell.cor, 0.3), -0.3), breaks=seq(-0.3, 0.3, length.out=100), color=shades, 
                                 show_rownames=F, show_colnames=F, treeheight_row = 0, treeheight_col = 0, 
				 annotation_row=cur.cls.df, annotation_colors=list(cur.cls=set3.cols), 
				 filename=file.path(hier.fig.dir, c('ecto_cor_heatmap.png', 'meso_cor_heatmap.png')[i]))
    tmp.ret.unnorm = pheatmap(pmax(pmin(cell.norm.cor.ret.unnorm$cell.cor, 0.3), -0.3), breaks=seq(-0.3, 0.3, length.out=100), color=shades, 
                                 show_rownames=F, show_colnames=F, treeheight_row = 0, treeheight_col = 0,
				 filename=file.path(hier.fig.dir, c('ecto_cor_heatmap_unnorm.png', 'meso_cor_heatmap_unnorm.png')[i]))

    # the next block is for correlations with the a score
    a.pooled = do.call(cbind, tapply(names(cur.cls), cur.cls, function(cells) colSums(ab.form$atab[cells,])))
    b.pooled = do.call(cbind, tapply(names(cur.cls), cur.cls, function(cells) colSums(ab.form$btab[cells,])))
    sel.rows = rownames(a.pooled)[rowSums(a.pooled + b.pooled) > 300]
    rel.a.scores = do.call(cbind, lapply(1:ncol(a.pooled), function(j) log2((a.pooled[,j] / (a.pooled[,j] + b.pooled[,j])) / (rowSums(a.pooled[,-j]) / rowSums(a.pooled[,-j] + b.pooled[,-j])))))

    sel.mcs = names(e9.exp@colors)[e9.exp@colors %in% sel.colors]
    legc.cur.mcs = log2(1e-5 + e9.exp@e_gc[,sel.mcs])
    cur.legc = do.call(cbind, tapply(colnames(legc.cur.mcs), e9.exp@colors[colnames(legc.cur.mcs)], function(mcs) rowMeans(legc.cur.mcs[,mcs,drop=F])))
    cur.legc.norm = cur.legc - rowMeans(cur.legc)

    deeply.covered.bins = rownames(unnorm.total.cov)[rowSums(unnorm.total.cov[, names(cur.cls)]) > 2000]
    for (chrom in c('chrM', 'chrX', 'chrY')) {
      deeply.covered.bins = deeply.covered.bins[!grepl(chrom, deeply.covered.bins)]
    }
    cur.chrom.cov = compute.chrom.cov(cur.cls, unnorm.total.cov[deeply.covered.bins,], decay.metrics)
    rel.cov = do.call(cbind, lapply(1:ncol(cur.chrom.cov[[1]]), function(j) log2(cur.chrom.cov[[1]][,j] / cur.chrom.cov[[2]][,j])))

    diff.thresh = c(1, 2)[i]
    genes.for.a = intersect(rownames(cur.legc)[gene.md[rownames(cur.legc), 'bin_formatted'] %in% sel.rows], rownames(cur.legc.norm)[apply(cur.legc.norm, 1, max) > diff.thresh])
    bins.for.a = gene.md[genes.for.a, 'bin_formatted']
    a.cors = matrix(NA, nrow=cur.num.clusters, ncol=length(sel.colors))
    for (hic.cls.index in 1:cur.num.clusters) {
      for (rna.cls.index in 1:length(sel.colors)) {
	a.cors[hic.cls.index, rna.cls.index] = cor(rel.a.scores[bins.for.a, hic.cls.index], cur.legc.norm[genes.for.a, rna.cls.index], method='spearman')
      }
    }
    colnames(a.cors) = sapply(colnames(cur.legc.norm), function(x) names(sel.colors)[sel.colors == x])
    rownames(a.cors) = 1:nrow(a.cors)
    cur.color.df = data.frame(cls_col=1:cur.num.clusters)
    rownames(cur.color.df) = 1:cur.num.clusters

    get.col.order <- function(mat) {
      return(order(apply(mat, 2, which.max)))
      selected.cols = apply(mat, 1, which.max)
      return(order(c(selected.cols, setdiff(1:ncol(mat), selected.cols))))


      selected.cols = rep(NA, nrow(mat))
      for (row.index in 1:nrow(mat)) {
	col.options = setdiff(1:ncol(mat), selected.cols[!is.na(selected.cols)])
	row.options = which(is.na(selected.cols))
        which.ret = which(mat[row.options, col.options, drop=F] == max(mat[row.options, col.options, drop=F]), arr.ind=T)
	selected.cols[row.options[which.ret[1,1]]] = col.options[which.ret[1, 2]]
      }
      selected.cols = c(selected.cols, setdiff(1:ncol(mat), selected.cols))
      return(selected.cols)
    }

    a.cors.col.ord = get.col.order(a.cors)
    png(file.path(hier.fig.dir, c('ecto_a_cors.png', 'meso_a_cors.png')[i]))
    pheatmap(pmax(pmin(a.cors, 0.3), -0.3)[,a.cors.col.ord], breaks=seq(-0.3, 0.3, length.out=100), color=shades, 
                                 show_rownames=F, show_colnames=T, treeheight_row = 0, treeheight_col = 0, 
				 cluster_rows=F, cluster_cols=F, annotation_row=cur.color.df, annotation_colors=list(cls_col=set3.cols), annotation_names_row=F)
    dev.off()


    # correlation between coverage and expression
    genes.for.cov = intersect(rownames(cur.legc)[gene.md[rownames(cur.legc), 'bin_formatted'] %in% rownames(cur.chrom.cov[[1]])], rownames(cur.legc.norm)[apply(cur.legc.norm, 1, max) > diff.thresh])
    bins.for.cov = gene.md[genes.for.cov, 'bin_formatted']
    cov.cors = matrix(NA, nrow=cur.num.clusters, ncol=length(sel.colors))
    for (hic.cls.index in 1:cur.num.clusters) {
      for (rna.cls.index in 1:length(sel.colors)) {
	cov.cors[hic.cls.index, rna.cls.index] = cor(rel.cov[bins.for.cov, hic.cls.index], cur.legc.norm[genes.for.cov, rna.cls.index], method='spearman')
      }
    }
    colnames(cov.cors) = sapply(colnames(cur.legc.norm), function(x) names(sel.colors)[sel.colors == x])
    rownames(cov.cors) = 1:nrow(cov.cors)
    cur.color.df = data.frame(cls_col=1:cur.num.clusters)
    rownames(cur.color.df) = 1:cur.num.clusters
    cov.cors.col.ord = get.col.order(cov.cors)
    png(file.path(hier.fig.dir, c('ecto_cov_cors.png', 'meso_cov_cors.png')[i]))
    pheatmap(pmax(pmin(cov.cors, 0.3)[,cov.cors.col.ord], -0.3), breaks=seq(-0.3, 0.3, length.out=100), color=shades, 
                                 show_rownames=F, show_colnames=T, treeheight_row = 0, treeheight_col = 0, 
				 cluster_rows=F, cluster_cols=F, annotation_row=cur.color.df, annotation_colors=list(cls_col=set3.cols), annotation_names_row=F)
    dev.off()


    # plot selected matchings between RNA cell types and Hi-C subclusters
    if (i == 1) {
      matchings.to.plot= list(c(2, 5))
    } else {
      matchings.to.plot = list(c(5, 9), c(6 , 4))
    }

    for (matching.index in seq_along(matchings.to.plot)) {
      cur.matching = matchings.to.plot[[matching.index]]
      png(file.path(hier.fig.dir, sprintf('%s_%s_%s.png', i, cur.matching[1], cur.matching[2])), width=900, height=400)
      par(mfrow=c(1, 2), mar=c(2, 2, 2, 2))
      cov.cor.val = cor(rel.cov[bins.for.cov, cur.matching[1]], cur.legc.norm[genes.for.cov, cur.matching[2]], method='spearman')
      cov.cor.formatted = format(round(cov.cor.val, 2), nsmall = 2) 
      plot(rel.cov[bins.for.cov, cur.matching[1]], cur.legc.norm[genes.for.cov, cur.matching[2]], pch=19, xlim=c(-0.55, 0.55), main=cov.cor.formatted)
      abline(v=0, h=0, col=2)
      a.cor.val = cor(rel.a.scores[bins.for.a, cur.matching[1]], cur.legc.norm[genes.for.a, cur.matching[2]], method='spearman')
      a.cor.formatted = format(round(a.cor.val, 2), nsmall = 2) 
      plot(rel.a.scores[bins.for.a, cur.matching[1]], cur.legc.norm[genes.for.a, cur.matching[2]], pch=19, xlim=c(-0.55, 0.55), main=a.cor.formatted)
      abline(v=0, h=0, col=2)
      dev.off()
    }

  }
}

plot.de.bins <- function(de.bins, ab.form, total.cov, repli.scores, cells1, cells2, fig.dir=get.fig.dir(), cols=1:2) {
  bins = intersect(intersect(colnames(ab.form$atab), rownames(total.cov)), rownames(repli.scores))
  de.bins1 = intersect(de.bins[[1]], bins)
  de.bins2 = intersect(de.bins[[2]], bins)
  a.scores1 = colSums(ab.form$atab[cells1, bins]) / (colSums(ab.form$atab[cells1, bins]) + colSums(ab.form$btab[cells1, bins]))
  a.scores2 = colSums(ab.form$atab[cells2, bins]) / (colSums(ab.form$atab[cells2, bins]) + colSums(ab.form$btab[cells2, bins]))
  cov.scores1 = rowSums(total.cov[bins, cells1]) / sum(total.cov[bins, cells1])
  cov.scores2 = rowSums(total.cov[bins, cells2]) / sum(total.cov[bins, cells2])
  bins.for.repli = rownames(repli.scores)[apply(repli.scores[,c(1,3)], 1, function(x) all(is.finite(x)))]
  repli.scores1 = repli.scores[bins.for.repli, 1]
  repli.scores2 = repli.scores[bins.for.repli, 3]
  names(repli.scores1) = bins.for.repli
  names(repli.scores2) = bins.for.repli

  for (i in 1:3) {
    cur.scores1 = list(cov.scores1, a.scores1, repli.scores1)[[i]]
    cur.scores2 = list(cov.scores2, a.scores2, repli.scores2)[[i]]
    attr.name = c('cov', 'a', 'repli')[i]
    png(file.path(fig.dir, sprintf('de_bins_%s_density.png', attr.name)))
    dens.all = density(cur.scores1 - cur.scores2)
    dens.de1 = density((cur.scores1 - cur.scores2)[de.bins1])
    dens.de2 = density((cur.scores1 - cur.scores2)[de.bins2])
    dens.ylim = range(sapply(list(dens.all, dens.de1, dens.de2), '[[', 'y'))
    if (i == 3) {
      plot(dens.all, lwd=5, ylim=dens.ylim, col='grey', xlim=c(-1.1, 1.1))
    } else {
      plot(dens.all, lwd=5, ylim=dens.ylim, col='grey')
    }
    lines(dens.de1, lwd=5, ylim=dens.ylim, col=cols[1])
    lines(dens.de2, lwd=5, ylim=dens.ylim, col=cols[2])
    dev.off()
  }
}

get.de.replicating.bins <- function(model.params, clusters.to.compare, timing.threshold=2) {
  stopifnot(length(clusters.to.compare) == 2)
  first.cluster = clusters.to.compare[1]
  second.cluster = clusters.to.compare[2]
  return(list(rownames(model.params$total.cov)[model.params$bin.clusters[, first.cluster] < model.params$bin.clusters[, second.cluster] - timing.threshold],
    rownames(model.params$total.cov)[model.params$bin.clusters[, second.cluster] < model.params$bin.clusters[, first.cluster] - timing.threshold]))
}

plot.cell.cor.and.2d.given.cors <- function(clusters, cell.cor, name, umap.scores, cls.col=1:length(table(clusters)), fig.dir=get.fig.dir()) {
  umap.config = umap.defaults
  umap.config$input = 'dist'
  umap.config$random_state = 42
  umap.ret = umap(1 - cell.cor, umap.config)
  umap.layout = umap.ret$layout

  png(file.path(fig.dir, paste0(name, '_umap.png')), height=1000, width=1200)
  num.cols = 200
  umap.scores = umap.scores[rownames(umap.layout)]
  color.grad = colorRampPalette(c("blue", "red"))(num.cols)
  score.cols = color.grad[round(max.min.rescale(c(umap.scores, 1, 2), 1, num.cols))]
  score.cols = head(score.cols, length(score.cols) - 2)
  layout(matrix(1:2,ncol=2), width = c(6,1),height = c(5,5))
  #par(mfrow=c(1, 2))
  plot(umap.layout[,1], umap.layout[,2], pch=19, cex=4, col=score.cols)
  legend_image <- as.raster(matrix(color.grad, ncol=1))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
  text(x=1.5, y = seq(0,1,l=5), labels = seq(1,2,l=5))
  rasterImage(legend_image, 0, 0, 1,1)
  dev.off()

  png(file.path(fig.dir, paste0(name, '_umap_colored.png')), height=1e3, width=1e3)
  plot(umap.layout[,1], umap.layout[,2], pch=19, cex=4, col=cls.col[clusters[rownames(umap.layout)]])
  dev.off()
}

plot.cell.bin.trends <- function(clusters, model.params, fig.dir, cols=1:length(unique(clusters))) {
  all.predicted.cov = predict.cov(model.params)$pred
  s.scores = model.params$s.scores
  names(s.scores) = names(clusters)
  s.range = range(s.scores)
  cov.range = quantile(rowMeans(model.params$total.cov), c(0.001, 0.999))
  nclusters = length(model.params$mixture.fractions)
  nbin.clusters = model.params$num.bin.clusters
  png(file.path(fig.dir, 'cell_bin_trends.png'), width=200 * nbin.clusters, height=200 * nclusters)
  par(mfrow=c(nclusters, nbin.clusters), mar=rep(0.5, 4))
  for (i in 1:nclusters) {
    cur.cells = clusters == i
    s.scores.ord = sort(s.scores[cur.cells])
    cells.ord = names(s.scores.ord)
    for (j in 1:nbin.clusters) {
      cur.bins = model.params$bin.clusters[,i] == j
      observed.cov = colMeans(model.params$total.cov[cur.bins, cur.cells])
      names(observed.cov) = names(clusters)[cur.cells]
      predicted.cov = colMeans(all.predicted.cov[cur.bins, cur.cells])
      names(predicted.cov) = names(observed.cov)
      #pred.cov = get.bin.expectation(s.score, j, num.bin.clusters=nbin.clusters, repl.duration=model.params$repl.duration)
      plot(s.scores.ord, observed.cov[cells.ord], pch=19, cex=2, col=cols[i], main='', ylab='', xlab='', yaxs='i', xaxs='i', xlim=s.range, ylim=cov.range, labels=F)
      lines(s.scores.ord, predicted.cov[cells.ord], col='black', lwd=3)
    }
  }
  dev.off()
}

plot.rs.for.clusters <- function(rs, fig.dir) {
  cex = 0.5
  rs = rs[apply(rs, 1, function(x) all(is.finite(x))),]
  #lim = range(rs)
  lim = range(c(-2, 2))

  png(file.path(fig.dir, 'emb_cls_rs_scatter.png'))
  plot(rs[,1], rs[,3], pch=19, cex=cex, col='white', ylim=lim, xlim=lim)
  grid(col='black')
  points(rs[,1], rs[,3], pch=19, cex=cex, col='darkblue')
  abline(b=1, a=0, col='red', lwd=3, lty='dashed')
  dev.off()

  png(file.path(fig.dir, 'emb_c22_vs_others_rs_scatter.png'))
  plot((rs[,1] + rs[,3]) / 2, rs[,2], pch=19, cex=cex, col='white', ylim=lim, xlim=lim)
  grid(col='black')
  points((rs[,1] + rs[,3]) / 2, rs[,2], pch=19, cex=cex, col='darkblue')
  abline(b=1, a=0, col='red', lwd=3, lty='dashed')
  dev.off()
}

plot.bin.trends <- function(clusters, model.params, decay.metrics, unnorm.total.cov, early.bin.cells=NULL, fig.dir=get.fig.dir(), cols=1:length(unique(clusters))) {
  g1.norm.total.cov = unnorm.total.cov / model.params$bin.probs
  total.cov = t(t(g1.norm.total.cov) / colSums(g1.norm.total.cov))
  num.bin.clusters = model.params$num.bin.clusters
  s.scores = model.params$s.scores
  names(s.scores) = names(clusters)
  
  nclust = length(unique(clusters))

  for (l in 1:num.bin.clusters) {
    png(file.path(fig.dir, sprintf('early_bin_trends_filtered_%d.png', l)), w=1600, h=500)
    par(mfrow=c(1, nclust))
    for (k in 1:nclust) {      
      early.bins = apply(model.params$bin.clusters, 1, function(x) {x[k] == l & all(x[k] < x[-k])})
      plot(s.scores[early.bin.cells], colSums(total.cov[early.bins, early.bin.cells]) / colSums(total.cov[, early.bin.cells]), col=cols[clusters[early.bin.cells]], pch=19, ylab='', cex=3)
    }
    dev.off()
  }
}

plot.cluster.density <- function(clusters, s.scores, cov.per.cell, name, fig.dir=get.fig.dir(), cols=1:length(unique(clusters))) {
  png(file.path(fig.dir, sprintf('cluster_density_s_score_%s.png', name)), w=2000, h=1000)
  sm.density.compare(s.scores, clusters, col=cols, lwd=15)
  dev.off()

  png(file.path(fig.dir, sprintf('cluster_density_num_umis_%s.png', name)), w=2000, h=1000)
  sm.density.compare(cov.per.cell[names(clusters)], clusters, col=cols, lwd=15)
  dev.off()
}

plot.model.params.comparison <- function(cls1, cls2, fig.dir=get.fig.dir()) {
  cls.table = table(cls1, cls2)
  png(file.path(fig.dir, 'cluster_comparison.png'), width=(ncol(cls.table) + 0.2) * 100, height=(nrow(cls.table) + 0.2) * 100)
  pheatmap(log2(1+cls.table), cluster_cols=T, cluster_rows=F, show_rownames=F, show_colnames=F, legend=T, clustering_distance_cols='correlation', cellwidth=100, cellheight=100, treeheight_col=0)
  dev.off()
}

