library("misha")
library("dplyr")
library("reshape2")
library("metacell")
library("tgconfig")
library("data.table")
library("parallel")
library("zoo")
library("Rgraphviz")
library("Matrix")


initial.hic.clustering <- function() {
  data.dir = get.data.file.dir()
  tgconfig::register_params(file.path(data.dir, 'schic_config.yaml'), package='schic')
  load.decay.metrics()
  
  set.seed(19)
  
  force_comp = F
  #force_comp = T
  
  env = schic_init_env()
  env = schic_init_perlim_cell_groups(env)
  env = schic_init_tss_intervs(env)
  message("will init bin cov")
  env = schic_init_bin_cov(env, force_comp=force_comp)
  message("will mask bin cov")
  env = schic_repli_karyo_spike_to_na(env)
  #schic_init_bin_cov(force_comp=force_comp,hres=T)
  #recompute 
  message("will cluster genome by rough repli")
  rep_gw = schic_repli_genome_clust(env, label="base", force_comp=force_comp) #add ab interv
  message("will cluster cells by rough repli")
  rep_gw_c_cl = schic_repli_cell_gw_clust(env, rep_gw, label="all_rep", force_comp=F)
  #rep_gw_c_cl = schic_repli_cell_gw_clust(env, rep_gw, label="all_rep", force_comp=T)
  
  env_emb_es = list(rep_gw=rep_gw, rep_gw_c_cl=rep_gw_c_cl)
  
  #compute cell ab
  message("will init cell,bin to strict A/B contacts")
  bin_cell_ab = schic_init_cell_bin_ab(env, rep_gw, "base", force_comp=force_comp)
  
  env_emb_es$bin_cell_ab = bin_cell_ab
  #ab per bin per ery/emb/escbin
  message("will marginalize A/B contacts on 3 repli cell clusters")
  ab_on_cell_clusts = schic_ab_bins_on_cell_clusts(env, bin_cell_ab, rep_gw_c_cl$cnv_k3) 
  env_emb_es$ab_on_cell_clusts = ab_on_cell_clusts
  
  message("will generate A/B locus cluster using the repl clust of emb, ery and ES ")
  ab_loc_clust = schic_ab_clust_loc_by_ab_on_cl(env, ab_on_cell_clusts, label="base", K=12, T_var=0.05, plot=T)
  
  env_emb_es$ab_loc_clust = ab_loc_clust
  
  message("computing ab down sampe ab stat on emb alone")
  ab_dsamp = schic_ab_downsamp(env, bin_cell_ab, label="base")
  
  message("will project AB on locus cluster for each cell")
  cell_on_ab_loc_clust = schic_ab_cell_on_loc_clust(env, ab_loc_clust=ab_loc_clust, tab_ds=ab_dsamp)
  
  env_emb_es$cell_on_ab_loc_clust = cell_on_ab_loc_clust
  
  message("will project all cells on 2D using PCA of the AB races")
  a = cell_on_ab_loc_clust$cell_rats
  an = t(t(a)/colMeans(a))
  
  
  nan.per.col = colSums(is.na(an))
  nan.cols = which(nan.per.col > 0)
  stopifnot(length(nan.cols) == 1)
  an = an[,-nan.cols]
  a2 = a[,-nan.cols]
  plot_cell_clust_using_ab(env, an, label='base')
  plot_cell_clust_using_ab(env, a2, label='base2')
  
  ab_cell_pca = prcomp(t(an)) #pc on cluster total per cell
  
  
  #define ery cells - write it
  message("Ery cells are defined by the ery-early and ery-late loci")
  
  first.index = min(setdiff(1:12, nan.cols))
  last.index = max(setdiff(1:12, nan.cols))
  # the offset may differ in different randomizations...
  ery.offset = 0.1
  ery_nms = rownames(cell_on_ab_loc_clust$cell_rats)[which(cell_on_ab_loc_clust$cell_rats[,last.index]-
  							cell_on_ab_loc_clust$cell_rats[,first.index] > ery.offset)]
  
  fig.dir = file.path(get.fig.dir(), 'initial_clustering')
  png(file.path(fig.dir, "ery_on_ab_cls.png"), w=600, h=600)
  
  plot(cell_on_ab_loc_clust$cell_rats[,last.index], 
  			cell_on_ab_loc_clust$cell_rats[,first.index], pch=21, 
  			bg="gray",xlab="Ery specific A", ylab="Ery specific B")
  points(cell_on_ab_loc_clust$cell_rats[ery_nms, last.index], 
  			cell_on_ab_loc_clust$cell_rats[ery_nms, first.index], cex=1,pch=21, bg="darkred")
  abline(a=-ery.offset, b=1)
  dev.off()
  
  if (!file.exists(ERYS.PATH)) {
    write.table(data.frame(ery_nms=ery_nms), quote=F, file=ERYS.PATH)
  } else {
    ery_nms = get.erys()
  }
  
  env$emb_cells_no_ery = setdiff(env$emb_cells, ery_nms)
  env$ery_cells = ery_nms
  
  cell.cycle.phase = c(env$decay_metrics$group, env$esc_decay_metrics$group)
  names(cell.cycle.phase) = c(rownames(env$decay_metrics), rownames(env$esc_decay_metrics))
  png(file.path(fig.dir, "rep_cell_pca_cell_cycle.png"), w=500, h=500)
  plot(ab_cell_pca$rotation[,1], ab_cell_pca$rotation[,2], 
       cex=0.8,pch=21,lwd=0.5, 
       bg=CELL.CYCLE.COLS[cell.cycle.phase[rownames(ab_cell_pca$rotation)]],
       xlab="PC 1", ylab="PC 2")
  dev.off()
  
  c.cls = c(rep(1, length(env$esc_cells)), rep(2, length(env$emb_cells_no_ery)), rep(3, length(env$ery_cells))) 
  names(c.cls) = c(env$esc_cells, env$emb_cells_no_ery, env$ery_cells)
  png(file.path(fig.dir, "rep_cell_pca_cls.png"), w=500, h=500)
  plot(ab_cell_pca$rotation[,1], ab_cell_pca$rotation[,2], 
       cex=0.8,pch=21,lwd=0.5, 
       bg=c("darkgreen","blue","red")[c.cls[rownames(ab_cell_pca$rotation)]],
       xlab="PC 1", ylab="PC 2")
  dev.off()
  
  message("will generate for each locus its repli and A/B state in emb, ery and ESC")
  strict_a = rownames(rep_gw$intervs_apx_ab)[rep_gw$intervs_apx_ab$ab_tor=="a"]
  strict_b = rownames(rep_gw$intervs_apx_ab)[rep_gw$intervs_apx_ab$ab_tor=="b"]
  cell_repA = colSums(env$bin_cell_cov_ds[strict_a,])
  cell_repB = colSums(env$bin_cell_cov_ds[strict_b,])
  
  cell_AB_cov = log2(cell_repA/cell_repB)
  G1_cells = intersect(names(which(cell_AB_cov < quantile(cell_AB_cov, 0.25))), 
  							colnames(env$bin_cell_cov_ds))
  esc_S_cells = names(which(rep_gw_c_cl$cnv_k3==1))
  emb_S_cells = names(which(rep_gw_c_cl$cnv_k3==2))
  ery_S_cells = names(which(rep_gw_c_cl$cnv_k3==3))
  esc_G1_cells = intersect(G1_cells, env$esc_cells)
  emb_G1_cells = intersect(G1_cells, env$emb_cells_no_ery)
  ery_G1_cells = intersect(G1_cells, env$ery_cells)
  
  emb_loc_SG1 = log2(rowMeans(env$bin_cell_cov_ds[, emb_S_cells])/
  					rowMeans(env$bin_cell_cov_ds[, emb_G1_cells]))
  esc_loc_SG1 = log2(rowMeans(env$bin_cell_cov_ds[, esc_S_cells])/
  					rowMeans(env$bin_cell_cov_ds[, esc_G1_cells]))
  ery_loc_SG1 = log2(rowMeans(env$bin_cell_cov_ds[, ery_S_cells])/
  					rowMeans(env$bin_cell_cov_ds[, ery_G1_cells]))
  ery_emb_loc_SG1 = log2(rowMeans(env$bin_cell_cov_ds[, c(ery_S_cells, emb_S_cells)]) /
  					rowMeans(env$bin_cell_cov_ds[, c(ery_G1_cells, emb_G1_cells)]))
  
  all_S_cells = c(ery_S_cells, esc_S_cells, emb_S_cells)
  all_loc_SG1 = log2(rowMeans(env$bin_cell_cov_ds[, all_S_cells])/
  					rowMeans(env$bin_cell_cov_ds[, G1_cells]))
  fil_S_cells = intersect(all_S_cells, colnames(ab_dsamp$atab_ds))
  all_loc_AB = rowSums(ab_dsamp$atab_ds[,fil_S_cells]) / (rowSums(ab_dsamp$atab_ds[,fil_S_cells]) + rowSums(ab_dsamp$btab_ds[,fil_S_cells]))
  g1_cell_cov = rowMeans(env$bin_cell_cov_ds[, intersect(G1_cells, colnames(env$bin_cell_cov_ds))])
  low_cov_bins = names(g1_cell_cov)[g1_cell_cov < 3]
  
  emb_S_cells = intersect(colnames(ab_dsamp$atab_ds), emb_S_cells)
  esc_S_cells = intersect(colnames(ab_dsamp$atab_ds), esc_S_cells)
  ery_S_cells = intersect(colnames(ab_dsamp$atab_ds), ery_S_cells)
  emb_loc_strictA = rowSums(ab_dsamp$atab_ds[,emb_S_cells])
  esc_loc_strictA = rowSums(ab_dsamp$atab_ds[,esc_S_cells])
  ery_loc_strictA = rowSums(ab_dsamp$atab_ds[,ery_S_cells])
  emb_loc_strictB = rowSums(ab_dsamp$btab_ds[,emb_S_cells])
  esc_loc_strictB = rowSums(ab_dsamp$btab_ds[,esc_S_cells])
  ery_loc_strictB = rowSums(ab_dsamp$btab_ds[,ery_S_cells])
  
  emb_loc_AB = emb_loc_strictA/(emb_loc_strictB+emb_loc_strictA)
  esc_loc_AB = esc_loc_strictA/(esc_loc_strictB+esc_loc_strictA)
  ery_loc_AB = ery_loc_strictA/(ery_loc_strictB+ery_loc_strictA)
  
  if (!file.exists(EMB.ERY.CHROM.DATA.PATH)) {
    save(ery_loc_AB, emb_loc_AB, ery_loc_SG1, emb_loc_SG1, file=EMB.ERY.CHROM.DATA.PATH)
  }
  if (!file.exists(EMB.ESC.CHROM.DATA.PATH)) {
    save(esc_loc_AB, emb_loc_AB, esc_loc_SG1, emb_loc_SG1, file=EMB.ESC.CHROM.DATA.PATH)
  }
  if (!file.exists(EMB.G1.CELLS.PATH)) {
    save(emb_G1_cells, ery_G1_cells, esc_G1_cells, file=EMB.G1.CELLS.PATH)
  }
  
  common_bins = setdiff(intersect(names(all_loc_AB), names(all_loc_SG1)), low_cov_bins)
  png(file.path(fig.dir, "SG1_vs_AB.png"))
  plot(all_loc_SG1[common_bins], all_loc_AB[common_bins], xlim=c(-1, 1), ylim=c(0, 1), col='darkblue', cex=0.5, pch=19)
  grid(col='grey', lwd=1)
  dev.off()
  
  strict_bins = c(strict_a, strict_b)
  png(file.path(fig.dir, "esc_vs_emb_ery_cov.png"), w=400, h=400)
  plot(esc_loc_SG1, ery_emb_loc_SG1, cex=0.3, pch=19, xlim=c(-1,1), ylim=c(-1,1), col="gray")
  abline(a=0,b=1,lty=2)
  points(esc_loc_SG1[strict_bins], ery_emb_loc_SG1[strict_bins], cex=0.5, pch=21,bg="black")
  dev.off()
  
  
  # Now code for the e10 data
  env.e10 = schic_init_env()
  e10_decay_metrics$cell = rownames(e10_decay_metrics)
  env.e10$decay_metrics = e10_decay_metrics
  env.e10$emb_cells = rownames(e10_decay_metrics)
  env.e10 = schic_init_perlim_cell_groups(env.e10)
  env.e10 = schic_init_tss_intervs(env.e10)
  env.e10 = schic_init_bin_cov(env.e10, force_comp=force_comp, 
                               fname='bin_cell_cov_e10', selected.cells=rownames(e10_decay_metrics))
  env.e10 = schic_repli_karyo_spike_to_na(env.e10)
  bin_cell_ab_e10 = schic_init_cell_bin_ab(env10, rep_gw, "e10", force_comp=force_comp, cells=rownames(e10_decay_metrics))

  # code for trans
  if (!file.exists(BIN.CELL.AB.TRANS.PATH)) {
    bin_cell_ab_trans = schic_init_cell_bin_ab(env, rep_gw, "trans", force_comp=force_comp, cis=F)
    save(bin_cell_ab_trans, file=BIN.CELL.AB.TRANS.PATH)
  }
  
  # Now code for hres
  set.seed(19)
  env.hres = schic_init_env()
  env.hres = schic_init_perlim_cell_groups(env.hres)
  env.hres = schic_init_tss_intervs(env.hres)
  env.hres = schic_init_bin_cov(env.hres, force_comp=force_comp, hres=T, fname='bin_cov_40k')

  bin_cell_cov_hres = env.hres$bin_cell_cov_hres
  bin_cell_cov_hres_ds = env.hres$bin_cell_cov_hres_ds

  bin_cell_ab_hres = schic_init_cell_bin_ab(env, rep_gw, "40k", force_comp=F, cov_binsize=4e4, is_hres=T)

  # repli_scores_c_clusters
  emb_loc_SG1_hres = log2(rowMeans(bin_cell_cov_hres_ds[, emb_S_cells]) / rowMeans(bin_cell_cov_hres_ds[, esc_G1_cells]))
  esc_loc_SG1_hres = log2(rowMeans(bin_cell_cov_hres_ds[, esc_S_cells]) / rowMeans(bin_cell_cov_hres_ds[, emb_G1_cells]))
  ery_loc_SG1_hres = log2(rowMeans(bin_cell_cov_hres_ds[, ery_S_cells]) / rowMeans(bin_cell_cov_hres_ds[, ery_G1_cells]))
  repli_scores_c_clusters = data.frame(emb=emb_loc_SG1_hres, ery=ery_loc_SG1_hres, esc=esc_loc_SG1_hres)
  if (!file.exists(REPLI.SCORES.C.CLUSTERS.PATH)) {
    write.table(repli_scores_c_clusters, file=REPLI.SCORES.C.CLUSTERS.PATH)
  }

  # ab_scores_c_clusters
  ab.form.hres = total.ab.format(names(c.cls), bin_cell_ab_hres)
  chrom.ab.hres = compute.chrom.ab(c.cls, ab.form.hres)[[1]]
  ab_scores_c_clusters = data.frame(emb=chrom.ab.hres[,2], ery=chrom.ab.hres[,3], esc=chrom.ab.hres[,1])
  if (!file.exists(AB.SCORES.C.CLUSTERS.PATH)) {
    write.table(ab_scores_c_clusters, file=AB.SCORES.C.CLUSTERS.PATH)
  }
}
