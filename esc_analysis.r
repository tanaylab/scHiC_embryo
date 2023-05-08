library(parallel)
library(dplyr)
library(shaman)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(metacell)
library(tidyr)
library(fields)
library(R.utils)
library(FNN)

get.common.decay.metrics <- function() {
  total.cov = bin_cell_cov
  emb.decay.metrics = sch_decay_metrics
  esc.decay.metrics = esc_sch_decay_metrics
  stopifnot('cell' %in% colnames(esc.decay.metrics))
  rownames(esc.decay.metrics) = esc.decay.metrics$cell
  esc.decay.metrics = esc.decay.metrics[, -1]
  common.cols = intersect(colnames(esc.decay.metrics), colnames(emb.decay.metrics))

  decay.metrics = rbind(emb.decay.metrics[, common.cols], esc.decay.metrics[, common.cols])
  return(decay.metrics)
}

get.erys <- function() {
  ery.path = ERYS.PATH
  erys = read.table(ery.path, stringsAsFactors=F)[,1]
  return(erys)
}

esc.analysis <- function() {
  set.misha(SCHIC.MISHA.PATH)
  exp_per_bin <<- get.or.create(EXP.PER.BIN.PATH, function() get.exp.per.bin(AB.SCORES.C.CLUSTERS.PATH))
  gene.md = get.gene.metadata.jumps()
  unnorm.total.cov = bin_cell_cov
  decay.metrics = get.common.decay.metrics()
  erys = get.erys()
  esc.cells = rownames(decay.metrics)[grep('scell', rownames(decay.metrics))]
  esc.cells = intersect(esc.cells, colnames(unnorm.total.cov))
  emb.cells = setdiff(rownames(sch_decay_metrics), erys)
  emb.cells = intersect(emb.cells, colnames(unnorm.total.cov))
  full.cluster.assignment = c(rep(1, length(emb.cells)), rep(2, length(esc.cells)))
  names(full.cluster.assignment) = c(emb.cells, esc.cells)

  fig.dir.emb.vs.esc = file.path(get.fig.dir(), 'emb_vs_esc_analysis')
  dir.create(fig.dir.emb.vs.esc, showWarnings=F)

  repl.data.dir = get.data.file.dir()
  dir.create(repl.data.dir, showWarnings=F)

  ab.form.path = file.path(repl.data.dir, 'ab_form_all_cells')
  ab.form = get.or.create(ab.form.path, function() total.ab.format(unique(bin_cell_ab$cell), bin_cell_ab))
  ab.form$atab = ab.form$atab[names(full.cluster.assignment),]
  ab.form$btab = ab.form$btab[names(full.cluster.assignment),]

  deeply.covered.bins1 = rownames(unnorm.total.cov)[rowMeans(unnorm.total.cov[, emb.cells]) > 8]
  deeply.covered.bins2 = rownames(unnorm.total.cov)[rowMeans(unnorm.total.cov[, esc.cells]) > 8]
  deeply.covered.bins = intersect(deeply.covered.bins1, deeply.covered.bins2)
  for (chrom in c('chrM', 'chrX', 'chrY')) {
    deeply.covered.bins = deeply.covered.bins[!grepl(chrom, deeply.covered.bins)]
  }

  sel.cells.for.cov = unlist(lapply(1:length(unique(decay.metrics$group)), function(i) {
    group.esc.cells = esc.cells[decay.metrics[esc.cells, 'group'] == i]
    group.emb.cells = emb.cells[decay.metrics[emb.cells, 'group'] == i]
    num.sel = min(length(group.esc.cells), length(group.emb.cells))
    sel.esc.cells = sample(group.esc.cells, num.sel)
    sel.emb.cells = sample(group.emb.cells, num.sel)
    return(c(sel.esc.cells, sel.emb.cells))
  }))

  full.cluster.assignment.repl = full.cluster.assignment[intersect(names(full.cluster.assignment), sel.cells.for.cov)]
  chrom.cov.esc.norm = compute.chrom.cov(full.cluster.assignment.repl, unnorm.total.cov[deeply.covered.bins,], decay.metrics[esc.cells,])
  chrom.cov.emb.norm = compute.chrom.cov(full.cluster.assignment.repl, unnorm.total.cov[deeply.covered.bins,], decay.metrics[emb.cells,])
  chrom.cov = chrom.cov.emb.norm
  chrom.cov[[1]][,2] = chrom.cov.esc.norm[[1]][,2]
  chrom.cov[[2]] = chrom.cov[[1]][,c(2, 1)]
  chrom.ab = compute.chrom.ab(full.cluster.assignment, ab.form, min.cov.discard=300, min.cov.per.clust=100, min.cov.cells=names(full.cluster.assignment))

  not.na.bins = intersect(names(emb_loc_SG1)[!is.na(emb_loc_SG1) & !is.na(esc_loc_SG1)], rownames(chrom.cov[[1]]))
  chrom.cov.norm = list(data.frame(emb_loc_SG1[not.na.bins], esc_loc_SG1[not.na.bins]), data.frame(esc_loc_SG1[not.na.bins], emb_loc_SG1[not.na.bins]))
  plot.loci.scatters.esc(chrom.cov.norm, chrom.ab, fig.dir.emb.vs.esc)

  track.names = paste0(c(EMB.POOL.TRACK.NAME, ESC.POOL.TRACK.NAME), '_equal_size')
  score.track.names = paste0(track.names, '_matshuff_scores')
  set.misha(SHAMAN.MISHA.PATH)
  emb.esc.tracks.conts = lapply(score.track.names, function(track.name) gextract(track.name, gintervals.2d.all(), colnames='score'))

  bins.for.4c = rownames(ab_scores_c_clusters)
  # Takes about 20 minutes
  all.4c = get.or.create(file.path(repl.data.dir, 'all_4c'), function() get.4c(emb.esc.tracks.conts, bins.for.4c))

  diff.per.bin = get.or.create(file.path(repl.data.dir, 'all_4c_shaman_diff'), function() all.4c.to.diff(all.4c))
  print('number of emb high shaman bins')
  print(sum(diff.per.bin < -50))
  print('number of esc high shaman bins')
  print(sum(diff.per.bin > 50))

  plot.esc.vs.emb.shaman(ab_scores_c_clusters, repli_scores_c_clusters, exp_per_bin, diff.per.bin, fig.dir=fig.dir.emb.vs.esc)

  insulation.dir = file.path(repl.data.dir, 'insulation')
  dir.create(insulation.dir, showWarnings=F)
  esc.emb.insulation = get.insulation(insulation.dir, track.names, scale=2e5, res=4e4)
  plot.insulation(esc.emb.insulation, fig.dir.emb.vs.esc)

  # plot Hi-C matrices
  esc.genes = c('Zfp42', 'Dppa2', 'Tet2', 'Arid5b', 'Klf9')
  all.interesting.bins = genes.to.intervs(esc.genes, gene.md, F)
  all.interesting.bins$geneSymbol = as.character(all.interesting.bins$geneSymbol)

  plot.dir.name = file.path(fig.dir.emb.vs.esc, 'hic_matrices')
  dir.create(plot.dir.name, showWarnings=F)
  chrom.ab.hres = list(ab_scores_c_clusters[, c('emb', 'esc')])
  for (i in 1:nrow(all.interesting.bins)) {
    interv = all.interesting.bins[i, ]
    forced.interv = gintervals.force_range(data.frame(chrom=interv$chrom, start=interv$start, end=interv$end))
    stopifnot(nrow(interv) == nrow(forced.interv))
    forced.interv$orig_start = interv$tss
    forced.interv$tss = interv$tss
    forced.interv$geneSymbol = interv$geneSymbol
    forced.interv$clust = 3 
    interv = forced.interv
    interv$chrom = as.character(interv$chrom)
    file.name = paste0(interv$geneSymbol, '_', interv$chrom, '_',
                         trimws(format(interv$start, scientific=F)), '.png')
    plot.file.name = file.path(plot.dir.name, file.name)
    plot.whole.chrom(chrom.ab.hres, chrom.ab.hres, chrom.ab.hres, emb.esc.tracks.conts, as.numeric(substr(interv$chrom, 4, 6)), plot.file.name, 
                     coord.range=c(interv$start, interv$end), cls.col=c('blue', 'darkgreen'), smooth.cov=F, height.vec=c(3, 3, 1.5), width=6, 
		     sel.attr.plots=c(2), show.coord.text=F, smooth.ab=T, rotate=T, mark.coord=interv$orig_start) 
  }

  # boxplots
  my_schic_ab_on_tsss(AB.SCORES.C.CLUSTERS.PATH, 'ab', fig.dir=get.fig.dir())
  my_schic_ab_on_tsss(REPLI.SCORES.C.CLUSTERS.PATH, 'repli', ylim=c(-1.5, 1.5), fig.dir=get.fig.dir())



  print('writing table S1')
  gene.md = get.gene.metadata.jumps()
  gene.md = gene.md[!duplicated(gene.md$geneSymbol),]
  rownames(gene.md) = gene.md$geneSymbol
  exp.per.cluster = get.exp.per.cluster(gene.md)

  gene.md$tss = gene.md$start
  all.gene.shaman = get.or.create(file.path(repl.data.dir, 'shaman_score_all_genes_signed'), function() get.gene.shaman.score.distrib(gene.md, track.names, 
           gene.names=rownames(exp.per.cluster), all.tracks.conts=emb.esc.tracks.conts, num.sampled.genes=0, is.abs=F))
  all.gene.shaman = as.data.frame(all.gene.shaman)

  all.gene.shaman$gene = rownames(exp.per.cluster)
  all.gene.shaman$high_conformation_change = apply(all.gene.shaman, 1, function(x) max(abs(as.numeric(x[!sapply(x, is.null)])), na.rm=T)) > 60
  all.gene.shaman$max_esc_exp = exp.per.cluster[,2]
  all.gene.shaman$max_embryo_exp = exp.per.cluster[,1]
  all.gene.shaman$embryo_esc_exp_diff = exp.per.cluster[,1] - exp.per.cluster[,2]
  all.gene.shaman$is_emb_high = all.gene.shaman$embryo_esc_exp_diff > 1
  all.gene.shaman$is_esc_high = all.gene.shaman$embryo_esc_exp_diff < -1

  all.gene.shaman.fil = all.gene.shaman[apply(all.gene.shaman, 1, function(x) {
    all(sapply(x, function(x) {if(length(x) == 0) return(F);length(x[[1]]) > 0}))
  }),]
  for (i in 1:15) {
    all.gene.shaman.fil[,i] = unlist(all.gene.shaman.fil[,i])
  }
  all.gene.shaman.fil = all.gene.shaman.fil[,c(16:22, 1:15)]
  write.csv(all.gene.shaman.fil[all.gene.shaman.fil$is_emb_high | all.gene.shaman.fil$is_esc_high,], file=file.path(repl.data.dir, 'shaman_score_all_genes_signed_table'))
}

get.exp.per.cluster <- function(gene.md) {
  exp.ret = get.or.create(file.path(repl.data.dir, 'atlases'), function() get.atlases(gene.md))
  e9.exp = exp.ret$e9.exp
  esc.exp = exp.ret$esc.exp

  common.exp.gene.names = union(rownames(e9.exp@e_gc), rownames(esc.exp@e_gc))
  e9.exp.egc = as.data.frame(e9.exp@e_gc)
  e9.exp.egc[setdiff(common.exp.gene.names, rownames(e9.exp.egc)),] = 0
  e9.exp.egc = as.matrix(e9.exp.egc[common.exp.gene.names,])
  esc.exp.egc = as.data.frame(esc.exp@e_gc)
  esc.exp.egc[setdiff(common.exp.gene.names, rownames(esc.exp.egc)),] = 0
  esc.exp.egc = as.matrix(esc.exp.egc[common.exp.gene.names,])
  exp.per.cluster = do.call(cbind, lapply(1:2, function(i) {
    ery.mcs = which(log2(e9.exp@mc_fp['Gata1',]) > 0.5)
    cluster.exp = apply(log2(list(e9.exp.egc[,-ery.mcs], esc.exp.egc)[[i]] + 1e-5), 1, max)
    return(cluster.exp)
  }))
  return(exp.per.cluster)
}

plot.loci.scatters.esc <- function(chrom.cov, chrom.ab, fig.dir, label='all', cex=0.5) {
  chrom.cov = chrom.cov[[1]]
  lim = range(chrom.cov)
  png(file.path(fig.dir, sprintf('esc_vs_emb_cov_%s.png', label)), width=400, height=400)
  plot(chrom.cov[,1], chrom.cov[,2], pch=19, cex=cex, col='darkblue', xlim=lim, ylim=lim)
  grid(col='black', lwd=1)
  abline(b=1, a=0, col='red', lwd=3, lty='dashed')
  dev.off()
  
  chrom.ab = chrom.ab[[1]]
  png(file.path(fig.dir, sprintf('esc_vs_emb_ab_%s.png', label)), width=400, height=400)
  plot(chrom.ab[,1], chrom.ab[,2], pch=19, cex=cex, col='darkblue', xlim=c(0,1), ylim=c(0,1))
  grid(col='black', lwd=1)
  abline(b=1, a=0, col='red', lwd=3, lty='dashed')
  dev.off()
}

plot.esc.vs.emb.shaman <- function(ab.scores, rep.scores, exp.per.bin, shaman.diff.per.bin, fig.dir=get.fig.dir(), data.dir=get.data.file.dir()) {

  exp.diff = log2(exp.per.bin$emb) - log2(exp.per.bin$esc)
  names(exp.diff) = rownames(exp.per.bin)
  exp.bin.type = ifelse(pmax(log2(exp.per.bin$esc), log2(exp.per.bin$emb)) < -18, 'inactive', 
                 ifelse(exp.diff < -2, 'esc_active',
                 ifelse(exp.diff > 2, 'emb_active', 'neutral')))
  names(exp.bin.type) = rownames(exp.per.bin)
  print('number of embryo and esc expressed bins')
  print(table(exp.bin.type))

  is.active = exp.bin.type %in% c('esc_active', 'emb_active')
  png(file.path(fig.dir, 'exp_per_bin_emb_vs_esc.png'))
  # hack to plot the points under the grid
  plot(log2(exp.per.bin$emb), log2(exp.per.bin$esc), pch=19, cex=0.5, col='white')
  grid(col='black')
  points(log2(exp.per.bin$emb)[!is.active], log2(exp.per.bin$esc)[!is.active], pch=19, cex=0.5, col='darkblue')
  points(log2(exp.per.bin$emb)[is.active], log2(exp.per.bin$esc)[is.active], pch=19, cex=0.5, col='darkblue')
  segments(-18, -20, -18, -30, col='red', lwd=2, lty='dashed')
  segments(-18, -20, 0, -2, col='red', lwd=2, lty='dashed')
  segments(-20, -18, -30, -18, col='red', lwd=2, lty='dashed')
  segments(-20, -18, -2, 0, col='red', lwd=2, lty='dashed')
  abline(b=1, a=0, col='red', lwd=3, lty='dashed')
  dev.off()

  ab.diff = ab.scores$emb - ab.scores$esc
  ab.bin.type = ifelse(ab.diff < -0.2, 'esc_active', 
                 ifelse(ab.diff > 0.2, 'emb_active', 'neutral'))
  names(ab.bin.type) = rownames(ab.scores)

  ab.scores.fil = ab.scores[, c('emb', 'esc')]
  ab.scores.fil = ab.scores.fil[apply(ab.scores.fil, 1, function(x) all(is.finite(x))),]
  print('fraction of < 10% non-changing bins')
  print(sum(abs(ab.scores.fil$emb - ab.scores.fil$esc) < 0.1) / nrow(ab.scores.fil))
  print('fraction of > 30% changing bins')
  print(sum(abs(ab.scores.fil$emb - ab.scores.fil$esc) > 0.3) / nrow(ab.scores.fil))

  png(file.path(fig.dir, 'ab_by_exp_type_boxplot.png'))
  ab.score.diff = ab.scores$esc - ab.scores$emb
  names(ab.score.diff) = rownames(ab.scores)
  cur.bins = intersect(rownames(ab.scores), rownames(exp.per.bin)[exp.bin.type %in% c('esc_active', 'emb_active')])
  boxplot(split(ab.score.diff[cur.bins], exp.bin.type[cur.bins]), las=2, col=c('blue', 'darkgreen'), ylim=c(-0.2, 0.2))
  dev.off()
  ab.diff.splitted = split(ab.score.diff[cur.bins], exp.bin.type[cur.bins])
  print('ks for ab difference with expression:')
  print(ks.test(ab.diff.splitted[[1]], ab.diff.splitted[[2]]))

  png(file.path(fig.dir, 'shaman_by_ab_type_density.png'), height=400, width=800)
  cur.bins = intersect(names(ab.bin.type), names(shaman.diff.per.bin))
  shaman.per.bin.type = split(shaman.diff.per.bin[cur.bins], ab.bin.type[cur.bins])
  print(names(shaman.per.bin.type))
  all.dens = lapply(shaman.per.bin.type, density)
  ylim=range(sapply(all.dens, "[", "y"), finite=T)
  plot(all.dens[[1]], col='blue', lwd=8, ylim=ylim)
  lines(all.dens[[2]], col='darkgreen', lwd=8, ylim=ylim)
  lines(all.dens[[3]], col='grey', lwd=8, ylim=ylim)
  dev.off()
  print('ks for shaman by expression')
  print(ks.test(shaman.per.bin.type[[1]], shaman.per.bin.type[[2]]))

  png(file.path(fig.dir, 'repli_by_exp_type_boxplot.png'))
  rep.score.diff = rep.scores$esc - rep.scores$emb
  names(rep.score.diff) = rownames(rep.scores)
  cur.bins = intersect(rownames(rep.scores), rownames(exp.per.bin)[exp.bin.type %in% c('esc_active', 'emb_active')])
  boxplot(split(rep.score.diff[cur.bins], exp.bin.type[cur.bins]), las=2, col=c('blue', 'darkgreen'), ylim=c(-1, 1))
  dev.off()
  rep.diff.splitted = split(rep.score.diff[cur.bins], exp.bin.type[cur.bins])
  print('ks for rep difference with expression:')
  print(ks.test(rep.diff.splitted[[1]], rep.diff.splitted[[2]]))

}

plot.insulation <- function(emb.esc.insulation, fig.dir) {
  png(file.path(fig.dir, 'insulation_unnorm.png'), height=600, width=600)
  #plot(emb.esc.insulation[,1], emb.esc.insulation[,2], col='darkblue', pch=19, ylim=c(-4, -1), xlim=c(-4, -1))
  plot(emb.esc.insulation[,1], emb.esc.insulation[,2], col='darkblue', pch=19, ylim=c(-7, 0), xlim=c(-7, 0))
  grid(col='grey', lwd=1)
  abline(b=1, a=0, col='red', lwd=3, lty='dashed')
  dev.off()
}

get.insulation <- function(insulation.dir, track.names, scale=6e5, res=2e4) {
  insulation.path = file.path(insulation.dir, sprintf('%s_%s', scale, res))
  print(insulation.path)
  all.ins.calc = get.or.create(insulation.path, function() {
    all.ins.calc = do.call(rbind, mclapply(1:19, function(i) {
      chrom.name = paste0('chr', i)
      ins = sch_calc_insu_ds_trian(chrom.name, scale=scale, res=res, track.names=track.names, tmp.insu.dir=insulation.dir)
      ins.calc = t(log(ins$a / (ins$a + ins$b1 + ins$b2)))
      #ins.calc = t(log(ins$a / (ins$b1 + ins$b2)))
      rownames(ins.calc) = paste0(chrom.name, '_', trimws(format(as.numeric(rownames(ins.calc)), scientific=F)))
      return(ins.calc)
    }, mc.cores=19)) 
    return(all.ins.calc)
  })
  bin.names = rownames(all.ins.calc)
  splitted = strsplit(bin.names, '_')
  bin.coords = as.numeric(sapply(splitted, function(x) x[[2]])) + 1
  bin.chroms = sapply(splitted, function(x) x[[1]])
  new.bin.names = paste0(bin.chroms, '_', trimws(format(bin.coords, scientific=F)))
  rownames(all.ins.calc) = new.bin.names
  return(all.ins.calc)
}

##################################
#
#
# insulation schematically: #contacts in A /  (A+B1+B2)
# store results in tables
#
#           --------
#           |   |B2/
#           | A | /
#           |___|/
#           |B1 /
#           |  /
#           | /
#           |/
#      
##################################

###########################################################################
# downsample K contacts independently in each triangle of a locus (A+B1+B2)
sch_calc_insu_ds_trian = function(chrom, scale, res, track.names, tmp.insu.dir, n_ds=30, coords=NULL, coords_nm="", discard_below=1000, return_data=T, rebuild=T) {
  
  sch_table_dir = tmp.insu.dir
  system(sprintf("mkdir -p %s/insu", sch_table_dir))

  tot = list(a = NULL, b1 = NULL, b2 = NULL, ds_a = NULL, ds_b1 = NULL, ds_b2 = NULL)
  fns = list()
  for (tn in names(tot)) {
    if (length(grep("ds", tn)) > 0) {
      fns[tn] = sprintf("%s/%s_insu_%d_every_%d_ds%d_%s%s.txt", sch_table_dir, chrom, scale, res, n_ds, tn, coords_nm)
    }
    else {  
      fns[tn] = sprintf("%s/%s_insu_%d_every_%d_%s%s.txt", sch_table_dir, chrom, scale, res, tn, coords_nm) 
    }
  }

  if (sum(!sapply(fns, file.exists)) == 0 && !rebuild) {
    for (tn in names(tot)) {
      #message(sprintf("Reading %s", fns[[tn]]))
      tot[[tn]] = read.table(fns[[tn]], header=T)
    }
  }
  else {
    chr_coords = gintervals(chrom, 0, -1)
    scope = gintervals.2d(chrom, 0, -1, chrom, 0, -1)

    if (is.null(coords)) {
      coords = seq(chr_coords$start + scale, chr_coords$end - scale - 1, by=res)
    }
    coords = pmin(coords, chr_coords$end - 1)
	iter_2d = gintervals.2d(chroms1 = chrom, 
			starts1=coords, 
			ends1=coords + 1,
			chroms2 = chrom, 
			starts2=coords,
			ends2=coords+1)
    
    for (nm in track.names) {
      message("Processing ", nm, "...")
      
      if(length(gvtrack.ls("obs_a")) == 1) {
		gvtrack.rm("obs_a")
      }
      if(length(gvtrack.ls("obs_b1")) == 1) {
		gvtrack.rm("obs_b1")
      }
      if(length(gvtrack.ls("obs_b2")) == 1) {
		gvtrack.rm("obs_b2")
      }
      gvtrack.create("obs_a", nm, "area")
      gvtrack.create("obs_b1", nm, "area")
      gvtrack.create("obs_b2", nm, "area")

      gvtrack.iterator.2d("obs_a", 
                          sshift1=-scale, eshift1=0, 
                          sshift2=0, eshift2=scale)

      gvtrack.iterator.2d("obs_b1", 
                          eshift1=-res, sshift1=-scale, 
                          eshift2=-res, sshift2=-scale)

      gvtrack.iterator.2d("obs_b2", 
                          eshift1=scale, sshift1=res, 
                          eshift2=scale, sshift2=res)

      ins = gextract("obs_a", "obs_b1", "obs_b2", scope, iterator=iter_2d, band=c(-scale*2,1))
      diag = gextract("obs_a", "obs_b1", "obs_b2", scope, iterator=iter_2d, band=c(-discard_below,1))
      st = data.frame(a = ins$obs_a - diag$obs_a, 
        b1 = ins$obs_b1 - diag$obs_b1,
        b2 = ins$obs_b2 - diag$obs_b2)

      ds_st = apply(st, 1, function(x) { if(!anyNA(x) && sum(x, na.rm=T) >= n_ds) { table(c(names(st), sample(rep(names(st), x), n_ds))) - 1 } else { rep(NA, 3) } })
      ds_st = as.data.frame(t(ds_st))
      colnames(ds_st) = names(st)
      
      tot$a = rbind(tot$a, st$a)
      tot$b1 = rbind(tot$b1, st$b1)
      tot$b2 = rbind(tot$b2, st$b2)
      
      tot$ds_a = rbind(tot$ds_a, ds_st$a)
      tot$ds_b1 = rbind(tot$ds_b1, ds_st$b1)
      tot$ds_b2 = rbind(tot$ds_b2, ds_st$b2)

    }
    for (tn in names(tot)) {
      rownames(tot[[tn]]) = track.names
      colnames(tot[[tn]]) = coords    
      write.table(tot[[tn]], fns[[tn]], quote=F, sep="\t")
    }
  }
  
  if (return_data) {
    return(tot)
  }
}

get.4c <- function(emb.esc.tracks.conts, bins.for.4c) {
  all.4c = do.call(rbind, mclapply(1:19, function(chrom.num) {
    chr.name = paste0('chr', chrom.num)
    cur.chrom.conts = lapply(emb.esc.tracks.conts, function(conts) filter(conts, conts$chrom1 == chr.name & conts$chrom2 == chr.name))
    cur.chr.bins = bins.for.4c[startsWith(bins.for.4c, paste0(chr.name, '_'))]
    print(paste('finished chrom conts', chrom.num))
    cur.chrom.4c = do.call(rbind, lapply(cur.chr.bins, function(bin) {
      splitted = strsplit(bin, '_')
      interv = list(chrom=chr.name, tss=as.numeric(splitted[[1]][2]))
      cur.interv.4c = get.4c.trace.mult.dist(cur.chrom.conts, interv, filter.chrom=F)
      if (nrow(cur.interv.4c) == 0) {
        return(NULL)
      }
      cur.interv.4c$bin = bin
      rownames(cur.interv.4c) = paste0(rownames(cur.interv.4c), '_', cur.interv.4c$bin)
      return(cur.interv.4c)
    }))
    print(paste('finished chrom', chrom.num))
    return(cur.chrom.4c)
  }, mc.cores=25))
  return(all.4c)
}


all.4c.to.diff <- function(all.4c) {
  dist.bins = c(-5e5, -3e5, -2e5, -1e5, -5e4, -2.5e4, 2.5e4, 5e4, 1e5, 2e5, 3e5, 5e5)
  all.bins = unique(all.4c$bin)
  shaman.scores.per.range = do.call(rbind, mclapply(all.bins, function(cur.bin) {
    tryCatch({
      gene.4c = filter(all.4c, bin == cur.bin)
      cont.bin = cut(gene.4c$dist, dist.bins)
      mean.diff.per.range = tapply(1:nrow(gene.4c), cont.bin, function(x) mean(gene.4c[x,2] - gene.4c[x,1], na.rm=T))
      stopifnot(length(mean.diff.per.range) == length(dist.bins) - 1)
      names(mean.diff.per.range) = NULL
      return(mean.diff.per.range)
    }, error=function(e) {
      print(c('error in', cur.bin))
      return(rep(NA, length(dist.bins) - 1))
    })
  }, mc.cores=25))
  max.diffs = apply(shaman.scores.per.range, 1, function(x) x[which.max(abs(x))])
  names(max.diffs) = all.bins
  return(unlist(max.diffs))
}

get.gene.shaman.score.distrib <- function(gene.md, track.names, gene.names=NULL, all.tracks.conts=NULL, num.sampled.genes=1e3, is.abs=T) {
  if (is.null(all.tracks.conts)) {
    all.tracks.conts = lapply(track.names, function(track.name) gextract(track.name, gintervals.2d.all(), colnames='score'))
  }
  set.seed(42)
  gene.options = unique(gene.md$geneSymbol)
  if (!is.null(gene.names)) {
    selected.genes = gene.names
    gene.options = setdiff(gene.options, gene.names)
  } else {
    selected.genes = c()
  }
  selected.genes = c(selected.genes, sample(gene.options, num.sampled.genes))
  gene.md = gene.md[selected.genes, ]

  dist.bins = c(-1e6 - 1, -7.5e5, -5e5, -3e5, -2e5, -1e5, -5e4, -2.5e4, 2.5e4, 5e4, 1e5, 2e5, 3e5, 5e5, 7.5e5, 1e6)
  shaman.background = do.call(rbind, mclapply(1:nrow(gene.md), function(i) {
    tryCatch({
      if (i %% 50 == 1) print(i)
      interv = gene.md[i, ]
      gene.4c = get.4c.trace.mult.dist(all.tracks.conts, interv, use.max.dist=T)
      coords = as.numeric(sapply(strsplit(rownames(gene.4c), '_'), function(x) x[2]))
      gene.4c$dist = coords - interv$tss
      gene.4c = gene.4c[between(abs(gene.4c$dist), -1e6, 1e6),]
      cont.bin = cut(gene.4c$dist, dist.bins)
      max.diff.per.bin = tapply(1:nrow(gene.4c), cont.bin, function(x) {
	gene.4c.diff = gene.4c[x,1] - gene.4c[x,2]
	if (is.abs) {
          return(max(abs(gene.4c.diff)))
	} else {
          return(gene.4c.diff[which.max(abs(gene.4c.diff))])
	}
      })
      stopifnot(length(max.diff.per.bin) == length(dist.bins) - 1)
      print(c('exitting ', i))
      return(max.diff.per.bin)
    }, error=function(e) {
      print(c('error in', i))
      return(rep(NA, length(dist.bins) - 1))
    })
  }, mc.cores=25)) 
  stopifnot(nrow(shaman.background) == length(selected.genes))
  rownames(shaman.background) = selected.genes
  return(shaman.background)
}

get.exp.per.bin <- function(tab_fn, bin.size=4e4) {
  set.misha(SCHIC.MISHA.PATH)
  ab40 = read.table(tab_fn, h=T, stringsAsFactors=F)
  env = schic_init_env()
  env = schic_init_perlim_cell_groups(env)
  env = schic_init_tss_intervs(env)
  rna_e9 = schic_rna_gen_abs_rna(env, bin_size=bin.size, 'e9_orig_bs500f')
  
  emb_max = apply(rna_e9$emb_bin_max, 1, max)
  emb_e = apply(rna_e9$emb_bin_max, 1, mean)
  
  esc_max = apply(rna_e9$esc_bin_max, 1, max)
  esc_e = apply(rna_e9$esc_bin_max, 1, mean)
  
  ery_max = apply(rna_e9$ery_bin_max, 1, max)
  ery_e = apply(rna_e9$ery_bin_max, 1, mean)
  
  emb_e = emb_e[rownames(ab40)]
  names(emb_e) = rownames(ab40)
  emb_e[is.na(emb_e)] = 2e-6
  
  esc_e = esc_e[rownames(ab40)]
  names(esc_e) = rownames(ab40)
  esc_e[is.na(esc_e)] = 2e-6
  
  ery_e = ery_e[rownames(ab40)]
  names(ery_e) = rownames(ab40)
  ery_e[is.na(ery_e)] = 2e-6
  
  
  all.exp = data.frame(emb=emb_e, ery=ery_e, esc=esc_e)
  rownames(all.exp) = names(emb_e)
  stopifnot(all(names(emb_e) == names(ery_e)))
  stopifnot(all(names(emb_e) == names(esc_e)))
  return(all.exp)
}


my_schic_ab_on_tsss = function(tab_fn, tag="ab", ylim=c(0,1), fig.dir=get.fig.dir(), bin.size=4e4)
{
	cur.exp.per.bin = get.exp.per.bin(tab_fn, bin.size=bin.size)
        ab40 = read.table(tab_fn, h=T, stringsAsFactors=F)
	emb_e = cur.exp.per.bin$emb
	esc_e = cur.exp.per.bin$esc
	ery_e = cur.exp.per.bin$ery
	names(emb_e) = rownames(cur.exp.per.bin)
	names(esc_e) = rownames(cur.exp.per.bin)
	names(ery_e) = rownames(cur.exp.per.bin)

	f_emb_up = (log2(emb_e) - log2(esc_e[names(emb_e)]))> 2
	f_emb_down = (log2(esc_e) - log2(emb_e[names(esc_e)])) >  2
	f_emb_cons = !f_emb_up & !f_emb_down
	f_emb_cons[is.na(f_emb_cons)] = F 
	f_emb_up[is.na(f_emb_up)] = F
	f_emb_down[is.na(f_emb_down)] = F

#can stratify genes by 
	png(sprintf("%s/emb_esc_%s_const.png", fig.dir, tag), w=600, h=400)
	boxplot(split(ab40[names(emb_e)[f_emb_cons],"emb"], cut(log2(emb_e[f_emb_cons]), -19:-9)), boxwex=0.3, col="blue", at=0.1+(1:10), las=2, ylim=ylim)
	boxplot(split(ab40[names(esc_e)[f_emb_cons],"esc"], cut(log2(esc_e[f_emb_cons]), -19:-9)), boxwex=0.3, col="darkgreen", at=0.5+(1:10),add=T, xaxt='n', yaxt='n')
	dev.off()

	png(sprintf("%s/emb_esc_%s_up.png", fig.dir, tag), w=600, h=400)
	boxplot(split(ab40[names(emb_e)[f_emb_up],"emb"], cut(log2(emb_e[f_emb_up]), -19:-9)), boxwex=0.3, col="blue", at=0.1+(1:10), las=2, ylim=ylim)
	boxplot(split(ab40[names(emb_e)[f_emb_up],"esc"], cut(log2(emb_e[f_emb_up]), -19:-9)), boxwex=0.3, col="darkgreen", at=0.5+(1:10),add=T, xaxt='n', yaxt='n')
	dev.off()

	png(sprintf("%s/emb_esc_%s_down.png",fig.dir, tag), w=600, h=400)
	boxplot(split(ab40[names(esc_e)[f_emb_down],"emb"], cut(log2(esc_e[f_emb_down]), -19:-9)), boxwex=0.3, col="blue", at=0.1+(1:10), las=2, ylim=ylim)
	boxplot(split(ab40[names(esc_e)[f_emb_down],"esc"], cut(log2(esc_e[f_emb_down]), -19:-9)), boxwex=0.3, col="darkgreen", at=0.5+(1:10),add=T, xaxt='n', yaxt='n')
	dev.off()

	f_ery_up = (log2(ery_e) - log2(emb_e[names(ery_e)]))> 2
	f_ery_down = (log2(emb_e) - log2(ery_e[names(emb_e)])) >  2
	f_ery_cons = !f_ery_up & !f_ery_down
	f_ery_cons[is.na(f_ery_cons)] = F 
	f_ery_up[is.na(f_ery_up)] = F
	f_ery_down[is.na(f_ery_down)] = F

#can stratify genes by 
	png(sprintf("%s/ery_emb_%s_const.png", fig.dir, tag), w=600, h=400)
	boxplot(split(ab40[names(ery_e)[f_ery_cons],"ery"], cut(log2(ery_e[f_ery_cons]), -19:-9)), boxwex=0.3, col="red", at=0.1+(1:10), las=2, ylim=ylim)
	boxplot(split(ab40[names(emb_e)[f_ery_cons],"emb"], cut(log2(emb_e[f_ery_cons]), -19:-9)), boxwex=0.3, col="blue", at=0.5+(1:10),add=T, xaxt='n', yaxt='n')
	dev.off()

	png(sprintf("%s/ery_emb_%s_up.png", fig.dir, tag), w=600, h=400)
	boxplot(split(ab40[names(ery_e)[f_ery_up],"ery"], cut(log2(ery_e[f_ery_up]), -19:-9)), boxwex=0.3, col="red", at=0.1+(1:10), las=2, ylim=ylim)
	boxplot(split(ab40[names(ery_e)[f_ery_up],"emb"], cut(log2(ery_e[f_ery_up]), -19:-9)), boxwex=0.3, col="blue", at=0.5+(1:10),add=T, xaxt='n', yaxt='n')
	dev.off()


	png(sprintf("%s/ery_emb_%s_down.png", fig.dir, tag), w=600, h=400)
	boxplot(split(ab40[names(emb_e)[f_ery_down],"ery"], cut(log2(emb_e[f_ery_down]), -19:-9)), boxwex=0.3, col="red", at=0.1+(1:10), las=2, ylim=ylim)
	boxplot(split(ab40[names(emb_e)[f_ery_down],"emb"], cut(log2(emb_e[f_ery_down]), -19:-9)), boxwex=0.3, col="blue", at=0.5+(1:10),add=T, xaxt='n', yaxt='n')
	dev.off()
}

