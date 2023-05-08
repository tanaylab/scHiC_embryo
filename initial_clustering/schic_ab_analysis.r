
#this gets as an input a genome wide cnv based locus clusters object as 
#computed by schic_repli_genome_clust, and aggragate contacts
#in cis and tarns from each genomic bin to the bona-fide strict A/B as
#determined by the genome clust.
schic_init_cell_bin_ab = function(env, cnv_gw, label, force_comp=F, min_dist=1e6, max_dist=1e10, cis=T, cells=NULL, cov_binsize=NULL, is_hres=F)
{
	data.dir = get.data.file.dir()
	fn = sprintf(file.path(data.dir, "%s.bin_ab.Rda"), label)
	if(!force_comp & file.exists(fn)) {
		load(fn)
		if (is_hres) {
		  return(bin_cell_ab_hres)
		} else {
		  return(bin_cell_ab)
		}
	} else {
		if (is.null(cov_binsize)) {
		  cov_binsize = get_param_strict("cov_binsize", "schic")
		}

		intervs_ab = cnv_gw$intervs_apx_ab
#		intervs_ab = sch_interv_tad_ab[,c("chrom","start","end", "ab_trans")]
#		intervs_ab$ab_trans = ifelse(intervs_ab$ab_trans > 0.5, 1, 0)

		split_chrom = function(x) {
			chr = as.character(x[1])
			max = as.numeric(x[3])
			bins = ceiling(max/cov_binsize)
			return(data.frame(chrom=rep(chr,times = bins),
								  start = seq(1, max, cov_binsize),
								  end = cov_binsize+  seq(1, max, cov_binsize)))
		}

		intervs_bin = do.call(rbind,apply(gintervals.all(), 1, split_chrom))
		intervs_bin = gintervals.force_range(intervs_bin)
		intervs_bin$f1 = paste(intervs_bin$chrom, intervs_bin$start, sep="_")
		rownames(intervs_bin) = intervs_bin$f1

#		stat =  schic_bin_mat_multi(c(emb_cells, esc_cells), 
		if (is.null(cells)) {
		  cells = c(env$emb_cells, env$esc_cells)
		}
		stat =  schic_bin_mat_multi(cells,
								intervs_bin, intervs_ab,
								"f1", "ab_tor", min_dist=min_dist, cis=cis, max_dist=max_dist)

		if (is_hres) {
		  bin_cell_ab = stat %>% dcast(f1+cell ~ ab_tor, value.var='cnt', fill=0)
		  save(bin_cell_ab, file=fn)
	          return(bin_cell_ab)
		} else {
		  bin_cell_ab_hres = stat %>% dcast(f1+cell ~ ab_tor, value.var='cnt', fill=0)
		  save(bin_cell_ab_hres, file=fn)
	          return(bin_cell_ab_hres)
		}
	}
}

schic_init_cell_bin_ab_direct = function(cells, intervs_ab, cov_binsize=4e4, min_dist=1e6) {
  split_chrom = function(x) {
  	chr = as.character(x[1])
  	max = as.numeric(x[3])
  	bins = ceiling(max/cov_binsize)
  	return(data.frame(chrom=rep(chr,times = bins),
  						  start = seq(1, max, cov_binsize),
  						  end = cov_binsize+  seq(1, max, cov_binsize)))
  }
  
  intervs_bin = do.call(rbind,apply(gintervals.all(), 1, split_chrom))
  intervs_bin = gintervals.force_range(intervs_bin)
  intervs_bin$f1 = paste(intervs_bin$chrom, intervs_bin$start, sep="_")
  rownames(intervs_bin) = intervs_bin$f1
  
  stat =  schic_bin_mat_multi(cells,
								intervs_bin, intervs_ab,
								"f1", "ab_tor", min_dist=min_dist, cis=T)

		bin_cell_ab = stat %>% dcast(f1+cell ~ ab_tor, value.var='cnt', fill=0)
	return(bin_cell_ab)
}

schic_ab_bins_on_cell_clusts = function(env, bin_cell_ab, cell_clust, label="base", force_recomp=F)
{
	data.dir = get.data.file.dir()
	fn = sprintf(file.path(data.dir, "%s.ab_on_c_cl.Rda"), label)
	if(!force_recomp & file.exists(fn)) {
		load(fn)
		return(ab_on_c_cl)
	}
#chrom marg (tapply)
   atab = bin_cell_ab  %>% dcast(cell ~ f1, value.var='a')
   btab = bin_cell_ab  %>% dcast(cell ~ f1, value.var='b')
	cell_nms = atab$cell
	atab = as.matrix(atab[,-1])
	btab = as.matrix(btab[,-1])
	rownames(atab) = cell_nms
	rownames(btab) = cell_nms
	
	atab[is.na(atab)] = 0
	btab[is.na(btab)] = 0

	f = colSums(atab)>100
	atab = atab[,f]
	btab = btab[,f]

	cell_tot_a = rowSums(atab)
	cell_tot_b = rowSums(btab)
	cell_rat  = log2(cell_tot_a/cell_tot_b)
	cell_rat_n = 2**(cell_rat - mean(cell_rat))
	atab = atab/cell_rat_n

	a_clust = tgs_matrix_tapply(t(atab[names(cell_clust),]), cell_clust, sum)
	b_clust = tgs_matrix_tapply(t(btab[names(cell_clust),]), cell_clust, sum)

	ds_nms = intersect(names(cell_clust), colnames(env$bin_cell_cov_ds))
	cnv_clust = t(tgs_matrix_tapply(env$bin_cell_cov_ds_mask_karyo[,ds_nms], cell_clust[ds_nms], function(x) mean(x, na.rm=T)))
	cnv_clust_lr = log2(cnv_clust/rowMeans(env$bin_cell_cov_ds))

	rownames(cnv_clust_lr) = rownames(env$bin_cell_cov_ds_mask_karyo)

	colnames(a_clust) = colnames(atab)
	colnames(b_clust) = colnames(btab)

	a_rat = t(a_clust/(ifelse((a_clust+b_clust)>20, a_clust+b_clust, NA)))

#clust and plotd
	a_rat_n = a_rat[!is.na(rowSums(a_rat)),]
	a_rat = a_rat[!is.na(rowSums(a_rat)),]
	a_rat_n = a_rat_n - rowMeans(a_rat_n)

	cnv_clust_lr = cnv_clust_lr[rownames(a_rat_n),]
	ab_on_c_cl = list(a_rat = a_rat, a_rat_n = a_rat_n,
				a_cl = a_clust, b_cl=b_clust, 
				cnv_clust_lr = cnv_clust_lr,
				cell_a_scale=cell_rat_n)
				

	save(ab_on_c_cl, file=fn)
	return(ab_on_c_cl)
}

schic_ab_clust_loc_by_ab_on_cl = function(env, ab_on_cl, label, K=12, T_var=0.1, plot=F)
{
	fig.dir = file.path(get.fig.dir(), 'initial_clustering')
	f_var = apply(abs(ab_on_cl$a_rat_n),1,max)>T_var
	a_rat_n = as.matrix(ab_on_cl$a_rat_n[f_var,])
	cnv_clust_lr = as.matrix(ab_on_cl$cnv_clust_lr[f_var,])
	cnv_clust_lr = cnv_clust_lr - rowMeans(cnv_clust_lr)
	browser()

#	km = kmeans(a_rat_n,K)
	km = kmeans(cbind(a_rat_n,cnv_clust_lr),K)

	cl_map = 1:K
	cl_map[order(km$centers[,3])] = 1:K
#	cl_map[order(km$centers[,2]-km$centers[,5])] = 1:K

	km$cluster = cl_map[km$cluster]
	names(km$cluster) = rownames(ab_on_cl$a_rat_n)[f_var]

	km = km$cluster
	
#plot heatmaps (hclust by cnv in each clust)
	ord_nms = c()
	for(i in 1:K) {
		#clust by cnv
		k_nms = names(which(km==i))
		cnvs = cnv_clust_lr[k_nms,]
		abs = a_rat_n[k_nms,]
		hc = hclust(tgs_dist(cbind(abs,cnvs)), "ward.D2")
		ord_nms = c(ord_nms, k_nms[hc$order])
	}

	toty = cumsum(table(km))/length(km)
#find genes
	ab_rat = a_rat_n[ord_nms,]
	cnv_rat = cnv_clust_lr[ord_nms,]

	if(plot) {
		png(sprintf(file.path(fig.dir, "%s_ab_on_cl.png"), label),h=1200, w=800)
		layout(matrix(1:2, nrow=1))
		par(mar=c(2,2,2,2))
		ab_rat = pmin(pmax(ab_rat, -0.5), 0.5)
		shades = colorRampPalette(c("darkblue","blue", "white","red", "darkred"))(1000)
		image.plot(t(ab_rat), col=shades, zlim=c(-0.5,0.5),xaxt='n', yaxt='n')
		for(i in 1:K) {
			abline(h=toty[i])
		}
		par(mar=c(2,0.2,2,2))
		shades = colorRampPalette(c("darkblue", "white","darkred"))(1000)
		cnv_rat = pmin(pmax(cnv_rat, -0.6), 0.6)
		image.plot(t(cnv_rat), col=shades, zlim=c(-0.6, 0.6),xaxt='n', yaxt='n')
		for(i in 1:K) {
			abline(h=toty[i])
		}
		dev.off()

		ab_rat_for_plot = ab_rat
		cnv_rat_for_plot = pmin(pmax(cnv_rat, -0.5), 0.5)
		shades = colorRampPalette(c("darkblue","blue", "white","red", "darkred"))(1000)
		fig_path = sprintf(file.path(fig.dir, "%s_ab_on_cl_pheatmap.png"), label)
		bin_cls = km[ord_nms]
		#png(fig_path ,h=3000, w=2000)
		png(fig_path ,h=650, w=440)
                pheatmap(cbind(ab_rat_for_plot, cnv_rat_for_plot), col=shades, cluster_rows=F, cluster_cols=F, breaks = seq(-0.5,0.5,l=1000),
                         show_rownames=F, show_colnames=F, annotation_names_row=F, annotation_names_col=F, annotation_legend=F,
                         #filename=fig_path, width=1600, height=2400,
	                 annotation_row=as.data.frame(bin_cls), annotation_colors=list(bin_cls=brewer.pal(K, 'Set3')))
	        dev.off()
	}
	bins = do.call(rbind, strsplit(names(km),"_"))
	clust = data.frame(chrom = bins[,1],
						start = as.numeric(bins[,2]),
						clust = km)
	rownames(clust) = names(km)
	ab_loc_clust = list(clust = clust)
	data.dir = get.data.file.dir()
	ab_loc_clust = schic_repli_loc_clust_diff_tss(env, ab_loc_clust, file.path(data.dir, sprintf("%s.ab_clust_tss_K%d.txt", label, K)))
	return(ab_loc_clust)

}

schic_ab_downsamp = function(env, bin_cell_ab, label="base", force_recomp=F, focus_cell=NULL)
{
	data.dir = get.data.file.dir()
	fn = file.path(data.dir, sprintf("%s.dsamp_ab.Rda", label))
	if(!force_recomp & file.exists(fn)) {
		load(fn)
		return(ab_tab_ds)
	} 
#chrom marg (tapply)
   atab = bin_cell_ab  %>% dcast(cell ~ f1, value.var='a', fill=0)
   btab = bin_cell_ab  %>% dcast(cell ~ f1, value.var='b', fill=0)
	if(!is.null(focus_cell)) {
		atab = atab[atab$cell %in% focus_cell,]
		btab = btab[btab$cell %in% focus_cell,]
	}

	tota = rowSums(atab[,-1])
	totb = rowSums(btab[,-1])
	n_a = round(quantile(tota,0.2))
	n_b = round(quantile(totb,0.2))
	f_cov = tota > n_a & totb > n_b
	atab_cov = t(atab[f_cov,-1])
	btab_cov = t(btab[f_cov,-1])
	
	atab_ds = flex_n_downsamp(atab_cov, n_a)
	btab_ds = flex_n_downsamp(btab_cov, n_b)
	colnames(atab_ds) = atab[f_cov,"cell"]
	colnames(btab_ds) = btab[f_cov,"cell"]

	ab_tab_ds = list(atab_ds = atab_ds, btab_ds = btab_ds)
	save(ab_tab_ds, file=fn)
	return(ab_tab_ds)
}

schic_ab_cell_on_loc_clust = function(env, ab_loc_clust, cell_bin_ab=NULL, tab_ds=NULL)
{
	if(!is.null(tab_ds)) {
		atab = tab_ds$atab_ds
		btab = tab_ds$btab_ds
		lclst = ab_loc_clust$clust[rownames(atab),"clust"]
		lclst[is.na(lclst)] = max(lclst,na.rm=T)+1
		a_tot = t(tgs_matrix_tapply(as.matrix(t(atab)), lclst, sum))
		b_tot = t(tgs_matrix_tapply(as.matrix(t(btab)), lclst, sum))
		cell_ratios = a_tot/(a_tot+b_tot)
		rownames(cell_ratios) = colnames(atab)
	} else if(!is.null(cell_bin_ab)) {
	   atab = cell_bin_ab  %>% dcast(cell ~ f1, value.var='a')
   	btab = cell_bin_ab  %>% dcast(cell ~ f1, value.var='b')
		atab[is.na(atab)]=0
		btab[is.na(btab)]=0
		lclst = ab_loc_clust$clust[colnames(atab)[-1],"clust"]
		lclst[is.na(lclst)] = max(lclst,na.rm=T)+1
		a_tot = t(tgs_matrix_tapply(as.matrix(atab[,-1]), lclst, sum))
		b_tot = t(tgs_matrix_tapply(as.matrix(btab[,-1]), lclst, sum))
		cell_ratios = a_tot/(a_tot+b_tot)
		rownames(cell_ratios) = atab$cell
	} else {
		stop("Must provide either ab_loc_clust or downsampled tab tab_ds")
	}
	return(list(cell_rats = cell_ratios))	
}

plot_cell_clust_using_ab <- function(env, cell_rats, label='base') {
	fig.dir = file.path(get.fig.dir(), 'initial_clustering')

	orig_cell_rats = cell_rats
	cell_rats = cell_rats[,-13]
	shades = colorRampPalette(c("darkblue","blue", "white","red", "darkred"))(1000)


	cell_rats = t(t(cell_rats) - colMeans(cell_rats, na.rm=T))
	cell_rats = t(t(cell_rats) / apply(cell_rats, 2, function(x) sd(x, na.rm=T)))
	cell_ab_cor = tgs_cor(t(cell_rats), pairwise.complete.obs=T)
	hc = hclust(as.dist(1-cell_ab_cor), "ward.D2")
        #shades = colorRampPalette(c("darkblue", "white","darkred"))(1000)
	png(file.path(fig.dir, sprintf("%s.cell_clust_using_ab.png",label)), w=500, h=1000)
	image(t(cell_rats[hc$ord,]), col=shades)
	dev.off()

	diag(cell_ab_cor) = 0
	k3 = cutree(hc, 3)
	#esc.cls = which.max(sapply(1:3, function(i) sum(names(k3)[k3 == i] %in% env$esc_cells)))
	#ery.cls = which.min(table(k3))
	#emb.cls = (1:3)[-c(esc.cls, ery.cls)]
	#cls.cols = c('darkseagreen', 'darkblue', 'red')[invPerm(c(esc.cls, emb.cls, ery.cls))][k3]
	cells.ord = rownames(cell_ab_cor)[hc$ord]
	esc.cls = which.max(sapply(1:3, function(i) sum(names(k3)[k3 == i] %in% env$esc_cells)))
	ery.cls = which.min(table(k3))
	emb.cls = (1:3)[-c(esc.cls, ery.cls)]
	cells.ord = c(cells.ord[k3[cells.ord] == esc.cls], cells.ord[k3[cells.ord] == emb.cls], cells.ord[k3[cells.ord] == ery.cls])
	k3 = k3[cells.ord]
	cls.cols = c('darkgreen', 'blue', 'red')[invPerm(c(esc.cls, emb.cls, ery.cls))][k3]
	#cls.cols = c('darkgreen', 'blue', 'red')[invPerm(c(3, 2, 1))][k3]
	names(cls.cols) = names(k3)

	png(file.path(fig.dir, sprintf("%s.cell_cor_using_ab_orig.png",label)), w=1000,h=1125)
	layout(matrix(1:2,ncol=1), h=c(1, 8))
	#layout(matrix(1:3,ncol=1), h=c(1, 1, 8))
	#par(mar=c(0,2,2,2))
	#common.feats = intersect(colnames(env$esc_decay_metrics), colnames(env$decay_metrics))
	#decay_metrics = rbind(env$esc_decay_metrics[,common.feats], env$decay_metrics[,common.feats])
	#plot(decay_metrics[cells.ord, 'f_mitotic_band'], xaxs='i', ylab='', xlab='', xaxt='n', yaxt='n', pch=19, cex=1.5)
	par(mar=c(0,2,2,2))
	plot(1:length(k3), rep(1, length(k3)), col=cls.cols[cells.ord], type='h', xaxs='i', ylab='', xlab='', lwd=5, xaxt='n', yaxt='n', ylim=c(0,1))

	#cell_ab_cor = pmin(pmax(cell_ab_cor, -0.4),0.4)
	#image(cell_ab_cor[cells.ord, cells.ord], col=shades,zlim=c(-0.4,0.4))
	par(mar=c(0,2,2,2))
	image.plot(cell_ab_cor[cells.ord, cells.ord], col=shades)
	dev.off()

	png(file.path(fig.dir, sprintf("%s.cell_clust_using_ab_raw.png",label)), w=500, h=1000)
	image.plot(orig_cell_rats[cells.ord, -12], col=shades)
	dev.off()

	png(file.path(fig.dir, sprintf("%s.cell_clust_using_ab_raw_heatmap.png",label)), w=500, h=1000)
        pheatmap(t(orig_cell_rats[cells.ord, -12]), col=shades, cluster_rows=F, cluster_cols=F,
                 show_rownames=F, show_colnames=F, annotation_names_row=F, annotation_names_col=F, annotation_legend=F,
	         annotation_col=as.data.frame(k3), annotation_colors=list(k3=c('darkgreen', 'blue', 'red')[invPerm(c(esc.cls, emb.cls, ery.cls))]))
	dev.off()
	return(k3)
}
