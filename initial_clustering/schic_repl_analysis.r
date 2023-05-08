
schic_init_bin_cov = function(env, force_comp=T, hres=F, selected.cells=NULL, fname='bin_cell_cov')
{
	data.dir = get.data.file.dir()
	dat_fn = file.path(data.dir, paste0(fname, ".Rda"))
	if(!force_comp & file.exists(dat_fn)) {
		if(hres) {
			load(dat_fn)
			env$bin_cell_cov_hres = bin_cell_cov_hres
			env$bin_cell_cov_hres_ds = bin_cell_cov_hres_ds
		} else {
			load(dat_fn)
			env$bin_cell_cov = bin_cell_cov
			env$bin_cell_cov_ds = bin_cell_cov_ds
		}
	} else {
		min_cov_ds = get_param_strict("min_cov_ds", "schic")
		if(hres) {
			cov_binsize = get_param_strict("cov_binsize_hres", "schic")
		} else {
			cov_binsize = get_param_strict("cov_binsize", "schic")
		}

		intervs_all = gintervals.all()
		intervs_all$f2 = "all"

		split_chrom = function(x) {
			chr = as.character(x[1])
			max = as.numeric(x[3])
			bins = ceiling(max/cov_binsize)
			return(data.frame(chrom=rep(chr,times = bins),
								  start = seq(1, max, cov_binsize),
								  end = cov_binsize+  seq(1, max, cov_binsize)))
		}

		intervs_bin = do.call(rbind,apply(intervs_all, 1, split_chrom))
		intervs_bin = gintervals.force_range(intervs_bin)
		intervs_bin$f1 = paste(intervs_bin$chrom, intervs_bin$start, sep="_")
		rownames(intervs_bin) = intervs_bin$f1

		if (is.null(selected.cells)) {
		  selected.cells = c(env$emb_cells, env$esc_cells)
		}
		cell_bin_all = schic_bin_mat_multi(selected.cells,
								intervs_bin, intervs_all,
								"f1", "f2", min_dist=0, cis=T)

		a = cell_bin_all %>% dcast(f1+f2 ~ cell, value.var='cnt', fill=0)
		cov= a[,c(-1,-2)]
		rownames(cov) = a$f1
		cov_ds = scm_downsamp(cov[,colSums(cov)>min_cov_ds], min_cov_ds)
		if(hres) {
			bin_cell_cov_hres = cov
			bin_cell_cov_hres_ds = cov_ds
			save(bin_cell_cov_hres, bin_cell_cov_hres_ds, file=dat_fn)
			env$bin_cell_cov_hres = cov
			env$bin_cell_cov_hres_ds = cov_ds
		} else {
			bin_cell_cov = cov
			bin_cell_cov_ds = cov_ds
			save(bin_cell_cov, bin_cell_cov_ds, file=dat_fn)
			env$bin_cell_cov = cov
			env$bin_cell_cov_ds = cov_ds
		}
	}
	return(env)
}

schic_repli_karyo_spike_to_na = function(env)
{
	chrom_map = sub("chr","", rownames(env$bin_cell_cov_ds))
	chrom_map = sub("_\\d+","", chrom_map)

	chrom_marg = tgs_matrix_tapply(t(env$bin_cell_cov_ds), chrom_map, sum)
	marg_n = chrom_marg/apply(chrom_marg, 1, median)
	env$bin_cell_cov_ds_mask_karyo = env$bin_cell_cov_ds
	env$bin_cell_cov_ds_mask_karyo[which(marg_n[chrom_map,]>1.3)]=NA

	env$bin_cell_cov_ds_mask_karyo = env$bin_cell_cov_ds_mask_karyo
	return(env)
}

schic_repli_genome_clust = function(env, label="base", force_comp=F)
{
	data.dir = get.data.file.dir()
	fig.dir = file.path(get.fig.dir(), 'initial_clustering')
	dir.create(fig.dir, showWarnings=F)
	fn = sprintf(file.path(data.dir, "%s.gw_cnvs.Rda"), label)
	if(!force_comp & file.exists(fn)) {
		load(fn)
		return(gw)
	}
#chrom marg (tapply)
	cov_ds = as.matrix(env$bin_cell_cov_ds_mask_karyo)

	f_cov = rowSums(env$bin_cell_cov_ds)>15000 & !grepl("chrX", rownames(cov_ds))

	loc_cnv_cor= tgs_cor(as.matrix(t(cov_ds[f_cov,])),pairwise.complete.obs=T)
	diag(loc_cnv_cor) = 0
	hc = hclust(tgs_dist(loc_cnv_cor),"ward.D2")

	png(file.path(fig.dir, "gw_cnv.png"), w=1200,h=1300)
	layout(matrix(1:3, ncol=1), h=c(1,10,1))
	par(mar=c(0,2,2,2))
	chrom_num = as.integer(as.factor(do.call(rbind, strsplit(rownames(loc_cnv_cor),"_"))[,1]))
	plot(chrom_num[hc$order], type="p", xaxs='i', cex=0.5, pch=19)
	par(mar=c(2,2,0,2))
	shades = colorRampPalette(c("darkblue","white","darkred"))(1000)
	loc_cnv_cor = pmax(pmin(loc_cnv_cor, 0.2),-0.2)
	image.plot(loc_cnv_cor[hc$order, hc$order], col=shades,zlim=c(-0.2,0.2))
	dev.off()
	nms = colnames(cov_ds)

	cl4=cutree(hc, 4)	
	cor_cl_k4 = tgs_matrix_tapply(
							tgs_matrix_tapply(loc_cnv_cor, cl4, mean), cl4, mean)
	a_id = ((which.min(cor_cl_k4)-1) %% 4)+1
	b_id = floor((which.min(cor_cl_k4)-1)/ 4)+1
	f = cl4==a_id | cl4==b_id
	bins = do.call(rbind, strsplit(names(cl4[f]),"_"))
	intervs_ab = data.frame(chrom=bins[,1], start = as.numeric(bins[,2]))
	intervs_ab$end = intervs_ab$start + 2e+5 - 1
	intervs_ab$ab_tor = ifelse(cl4[f]==a_id, "a", "b")
	intervs_ab = gintervals.force_range(intervs_ab)
	rownames(intervs_ab) = names(cl4[f])
	gw = list(gw_hc = hc, intervs_apx_ab = intervs_ab)

	png(file.path(fig.dir, "gw_cnv_clusters_legend.png"), w=1200,h=1300)
	cells.ord = rownames(loc_cnv_cor[hc$order, hc$order])
	plot(rep(1, length(cells.ord)), type='h', col=cl4[cells.ord], xaxt='s', xaxs='i')
	dev.off()

	save(gw, file=fn)
	return(gw)
}

schic_repli_cell_gw_clust = function(env, gw_loc, label="base", two_phase_norm=F, exclude=NULL, force_comp=F, draw_cluster_assignment=T)
{
	data.dir = get.data.file.dir()
	fig.dir = file.path(get.fig.dir(), 'initial_clustering')
	fn = sprintf(file.path(data.dir, "%s.gw_c_cor_cnvs.Rda"), label)
	fn_bins = sprintf(file.path(data.dir, "%s.early_late_bins.Rda"), label)
	if(!force_comp & file.exists(fn)) {
		load(fn)
		return(gw_c_cor)
	}

#	el_loci = c(gw_loc$early_loc, gw_loc$late_loc)
	strict_A = rownames(gw_loc$intervs_apx_ab)[gw_loc$intervs_apx_ab$ab_to=="a"]
	strict_B = rownames(gw_loc$intervs_apx_ab)[gw_loc$intervs_apx_ab$ab_to=="b"]
	early = colSums(env$bin_cell_cov_ds[strict_A,])
	late = colSums(env$bin_cell_cov_ds[strict_B,])

	cell_el_ratio = log2(early/late)
	early_late_bins = list(strict_A=strict_A, strict_B=strict_B)
	save(early_late_bins, file=fn_bins)

	cnms = names(cell_el_ratio)
	png(file.path(fig.dir, "EL_ratios.png"), w=400, h=400)
	plot(density(cell_el_ratio[intersect(cnms, env$esc_cells)]), lwd=4, col="darkseagreen")
	lines(density(cell_el_ratio[intersect(cnms, env$emb_cells)]), lwd=4, col="darkblue")
	abline(v=1.8, lwd=2, lty=2)
	dev.off()

	all_repli_cells = names(which(cell_el_ratio > 1.8))
	all_repli_cells = intersect(all_repli_cells, colnames(env$bin_cell_cov_ds))
	message("retaining ", length(all_repli_cells), " repl cells for clustering")

	f_cov = rowSums(env$bin_cell_cov_ds)>15000 &
							 !grepl("chrX", rownames(env$bin_cell_cov_ds)) 
#							 !rownames(env$bin_cell_cov_ds) %in% el_loci

	mat = env$bin_cell_cov_ds[f_cov, all_repli_cells]
#	mat = mat[,order(early/late)]
	mat = mat[,order(cell_el_ratio[colnames(mat)])]
	cell_cor= tgs_cor(as.matrix(mat))
	cell_cor_trend = t(apply(cell_cor, 1, rollmedian, 101, fill="extend"))
	cell_cor_n2 = cell_cor - cell_cor_trend

	cell_cor_n2 = cell_cor_n2+t(cell_cor_n2)

	diag(cell_cor_n2) = 0
#	hc = hclust(tgs_dist(cell_cor_n2), "ward.D2")
	hc = hclust(as.dist(1-cell_cor_n2), "ward.D2")
	k3 = cutree(hc, 3)

	nms = rownames(cell_cor_n2)

	if (draw_cluster_assignment) {
	  png(sprintf(file.path(fig.dir, "%s.gw_cell_cor_with_assignment.png"), label), w=4800, h=6000)
	  layout(matrix(1:2,ncol=1), h=c(2,8))
	  par(mar= c(0,2,2,2))
	  esc.cls = which.max(sapply(1:3, function(i) sum(names(k3)[k3 == i] %in% env$esc_cells)))
	  ery.cls = which.min(table(k3))
	  emb.cls = (1:3)[-c(esc.cls, ery.cls)]
	  cls.cols = c('darkgreen', 'blue', 'red')[invPerm(c(esc.cls, emb.cls, ery.cls))][k3]
	  names(cls.cols) = names(k3)

          cells.ord = rownames(cell_cor_n2)[hc$order]
	  cells.ord = c(cells.ord[k3[cells.ord] == esc.cls], cells.ord[k3[cells.ord] == emb.cls], cells.ord[k3[cells.ord] == ery.cls])
	  plot(1:length(k3), rep(1, length(k3)), col=cls.cols[cells.ord], type='h', xaxs='i', ylab='', xlab='', lwd=10, xaxt='n', yaxt='n')
	  shades = colorRampPalette(c("brown", "darkblue", "white", "red", "yellow"))(1000)
	  cell_cor_n2 = pmin(pmax(cell_cor_n2, -0.4),0.4)
	  image.plot(cell_cor_n2[cells.ord, cells.ord], col=shades,zlim=c(-0.4,0.4), legend.width=5)
	  dev.off()
	} else {
	  png(sprintf(file.path(fig.dir, "%s.gw_cell_cor.png"), label), w=1600, h=2000)
	  layout(matrix(1:3,ncol=1), h=c(1,1,8))
	  par(mar= c(0,2,2,2))
	  early_n = log2(early/mean(early))
	  late_n = log2(late/mean(late))
	  plot(seq(0,1,len=ncol(cell_cor_n2)),early_n[nms[hc$order]], cex=0.3, pch=19,xaxs='i',col="red")
	  points(seq(0,1,len=ncol(cell_cor_n2)),late_n[nms[hc$order]], cex=0.3, pch=19,col="blue")
	  par(mar= c(0,2,0,2))
	  grp = rep(0, ncol(cell_cor_n2))
	  names(grp) = nms
	  nms1 = intersect(nms, rownames(env$decay_metrics))
	  nms2 = intersect(nms, rownames(env$esc_decay_metrics))
	  grp[nms1] = env$decay_metrics[nms1, "group"]+6
	  grp[nms2] = env$esc_decay_metrics[nms2, "group"]
	  grp_col = c(colorRampPalette(c("lightgreen", "darkgreen"))(6),
	  					colorRampPalette(c("lightblue", "darkblue"))(6))

	  plot(grp[nms[hc$order]], cex=0.3, pch=19, col=grp_col[grp[hc$order]], xaxs='i')

	  shades = colorRampPalette(c("brown","darkblue","white","red","yellow"))(1000)
	  cell_cor_n2 = pmin(pmax(cell_cor_n2, -0.4),0.4)
	  image.plot(cell_cor_n2[hc$order, hc$order], col=shades,zlim=c(-0.4,0.4), legend.width=5)
	  dev.off()
	}

	loc_cnv_on_k3 = t(tgs_matrix_tapply(env$bin_cell_cov_ds[,names(k3)], k3, sum))
	f_cov = apply(loc_cnv_on_k3, 1, min) > 150
	cnv_prof3 = t(t(loc_cnv_on_k3[f_cov,])/colSums(loc_cnv_on_k3))
	rownames(cnv_prof3) = rownames(env$bin_cell_cov_ds)[f_cov]

	cnv_prof3 = log2(cnv_prof3/apply(cnv_prof3, 1, mean))
	cnv_prof3 = cnv_prof3[apply(abs(cnv_prof3),1,max)>0.2,]

	shades = colorRampPalette(c("darkblue","blue", "white","red", "darkred"))(1000)
	K = 12
	km = kmeans(cnv_prof3, K)	
	cl_map = 1:K
	cl_map[order(km$centers[,3])] = 1:K
#	cl_map[order(km$centers[,2]-km$centers[,5])] = 1:K

	km$cluster = cl_map[km$cluster]
	names(km$cluster) = rownames(cnv_prof3)

	loc_km12 = km$cluster

	png(sprintf(file.path(fig.dir, "%s.k12_loc_cor.png"), label), w=500,h=1000)
	image.plot(t(cnv_prof3[order(loc_km12),]), col=shades)
	dev.off()

	cnv_m = as.matrix(t(env$bin_cell_cov_ds[names(loc_km12),names(k3)])) 
	cell_on_cnv_k12 = t(tgs_matrix_tapply(cnv_m, loc_km12, sum))
	rownames(cell_on_cnv_k12) = names(k3)

	gw_c_cor = list(hc=hc, 
							cnv_k3 = k3, 
							cors=cell_cor_n2, 
							loc_on_c3_k12 = 
							loc_km12, 
							cell_on_cnv_k12 = cell_on_cnv_k12, 
							cnv_prof3 = cnv_prof3)
	save(gw_c_cor, file=fn)
	return(gw_c_cor)
}

schic_repli_loc_clust_diff_tss = function(env, loc_clust, tab_fn)
{
	bins = loc_clust$clust
	bins$end = bins$start + 2e+5
	bins$id = rownames(loc_clust$clust)
	bins = bins[,c("chrom","start","end","id", "clust")]
	bins = gintervals.force_range(bins)

	gene_bin = gintervals.neighbors(env$tss_intervs_25k, bins, maxdist=0)

	x = data.frame(tapply(gene_bin$gene, gene_bin$id, paste, sep=" "))
	colnames(x) = "genes"
	x$id = rownames(x)
	y = bins %>% left_join(x)

	y$genes = as.character(y$genes)
	write.table(as.data.frame(y[order(y$clust),]), file=tab_fn, quote=F, sep="\t")

	loc_clust$tss_clust = gene_bin
	loc_clust$tss_on_bin = y
	return(loc_clust)
}
