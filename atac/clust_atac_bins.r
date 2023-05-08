
clust_atac_bins = function(mm, k_n=100, T_vclst_min=-16)
{
	a_legc_n = mm$a_legc - rowMeans(mm$a_legc)

	km = tglkmeans::TGL_kmeans(a_legc_n, k=k_n, id_column=F,seed=19)

	rownames(km$centers) = 1:k_n

	message("done kmeans")
	f_big_km = table(km$cluster)>99

	clst_cent_abs = tgs_matrix_tapply(t(mm$a_legc), km$cluster, mean) 

	colnames(clst_cent_abs) = colnames(mm$a_legc)

	clst_min = apply(clst_cent_abs,1,function(x) sort(x)[4])
	clst_max = apply(clst_cent_abs,1,max)
	fvar_clst = (clst_max-clst_min)>0.7 & clst_min < T_vclst_min & f_big_km
	vclst_nms = which(fvar_clst)
	vclst_cent_abs = clst_cent_abs[fvar_clst,]

	peak_cls_col = mm$mclst_color[apply(clst_cent_abs,1,which.max)]

#km = kmeans of peaks to clusters
#clst_cent_abs - the mean (of logs) of atac values per peak cluster and clmc
#vclst_cent_abs - only the variable clusters heatmap
#peak_cls_col - we are not really using it (the "cell type" of the peak cluster)
#vclst_nms - the cluster numbers of those that are "variable"
	atac_clsts = list(km=km, 
							clst_cent_abs = clst_cent_abs, 
							vclst_cent_abs = vclst_cent_abs, 
							paek_cls_col = peak_cls_col,
							vclst_nms = vclst_nms)
	save(atac_clsts, file=file.path(ATAC.DATA.DIR, 'atac_clsts.Rda'))
	return(atac_clsts)
}

plot_atac_clusts = function(mm, ac, k_sup=14)
{
	fold_shades = colorRampPalette(c("darkblue","white","darkred"))(1000)
	wfold_shades = colorRampPalette(c("darkblue","blue","lightblue","white","pink","red","darkred"))(1000)
	abs_shades = colorRampPalette(c("white", "lightblue", "blue", "yellow", "red", "darkred"))(1000)

	f_big_km = table(ac$km$cluster)>99

	ol = order_clmc_by_tfs(mm, ac$vclst_cent_abs)
	ord_clmc = ol$ord_clmc
	hc_clmc = ol$hc_clmc
	tfs_top_cor = ol$tfs_top_cor

	hc_vclst = hclust(as.dist(1-tgs_cor(t(ac$vclst_cent_abs))),"ward.D2")
	vclst_ord = hc_vclst$order

#	vclst_ord = order(apply(vclst_cent_abs[,ord_clmc],1,which.max))

	png(file.path(ATAC.FIGURE.DIR, 'rna_types.png'), w=2600, h=400)
	barplot(rep(1,ncol(mm$a_legc)), col=mm$mclst_color[ord_clmc], border=NA)
	dev.off()

	hc_pclst = hclust(tgs_dist(ac$km$center[f_big_km,]), "ward.D2")

	pheatmap::pheatmap(t(ac$km$center[f_big_km,][hc_pclst$order, ord_clmc]),
			fontsize=5, breaks=seq(-1,1,l=999), 
			col=wfold_shades, filename=file.path(ATAC.FIGURE.DIR, 'atac_km.png'), 
			w=12,h=12, 
			cluster_rows=F, cluster_cols=F)
	png(file.path(ATAC.FIGURE.DIR, 'km_sizes.png'), w=1500,h=200)
	barplot(table(ac$km$cluster)[f_big_km][hc_pclst$order],las=2)
	dev.off()
	ashades =colorRampPalette(c("white","gray","lightblue", "blue", "darkblue","brown","yellow"))(1000)
	pheatmap::pheatmap(t(ac$clst_cent_abs[f_big_km,][hc_pclst$order,  ord_clmc]), 
				fontsize=5, 
				col=ashades,
				filename=file.path(ATAC.FIGURE.DIR, 'atac_km_abs.png'), w=12,h=12, 
				cluster_rows=F, cluster_col=F)

	pheatmap::pheatmap(t(ac$vclst_cent_abs[vclst_ord,  ord_clmc]), 
				fontsize=5, 
				col=ashades,
				filename=file.path(ATAC.FIGURE.DIR, 'atac_vkm_abs.png'), w=12,h=12, 
				cluster_rows=F, cluster_col=F)

	f_var_peak = ac$km$cluster %in% ac$vclst_nms
	legc_act_v = mm$a_legc[f_var_peak,]
	map_cls = 1:length(ac$vclst_nms)
	names(map_cls) = ac$vclst_nms[vclst_ord]
	peak_ord = order(map_cls[as.character(ac$km$cluster[f_var_peak])])
	pheatmap::pheatmap(t(legc_act_v[peak_ord, ord_clmc]), 
			col=abs_shades, breaks=c(-16.6,seq(-16,-13.8,l=999),-11), 
			fontsize=6, w=12,h=12, 
			filename=file.path(ATAC.FIGURE.DIR, 'atac_all_var_km.png'), 
			cluster_rows=F, cluster_cols=F)

	sup_clst = cutree(hc_vclst, k_sup)

#ordering the rna metacel clusters in a "Reasonable" way
#vclst_ord - visualization order of the peak clust orders
#hc_vclst - hierarchical clustering of the vclsts peaks
#sup_clst - grouping vclst into fewer super clusters
#f_var_peaks - factor defining the peaks that are part of variable clusters 
#peak_ord - ordering of the peak themselves
#tfs_top_Cor - the tfs that are most correlated (in RNA) with each peak cluster
	acn = list(ord_clmc=ord_clmc, 
							hc_clmc = hc_clmc,
							vclst_ord = vclst_ord,
							hc_vclst = hc_vclst,
							sup_clst = sup_clst,
							f_var_peak = f_var_peak,
							peak_ord = peak_ord,
							tfs_top_cor = tfs_top_cor)

	acn = project_interv_mm9(mm, ac, acn)
	acn = project_ab_on_atac_interv(mm, ac, acn)

	atac_clsts_annot = acn
	
	save(acn, file=file.path(ATAC.DATA.DIR, 'atac_clsts_annot.Rda'))
	return(atac_clsts_annot)
}

plot_atac_supclst_tfs = function(mm, ac, acn,  
										sup_cl, 
										foc_tfs=NULL, width=1200) 
{
	tfs_legc = mm$r_legc[acn$all_tfs,]
	legc_act_v = mm$a_legc[acn$f_var_peak,]

	tfs_legc = mm$r_legc[names(acn$tfs_top_cor),]
	tfs_lfp = tfs_legc - rowMeans(tfs_legc)

	pos_shades = colorRampPalette(c("white","white", "white","pink","red","darkred"))(1000)
	abs_shades = colorRampPalette(c("white", "lightblue", "blue", "yellow", "red", "darkred"))(1000)

	message("plot cl ")	
	message("cls ", paste(sup_cl, collapse=" "))
	f = ac$km$cluster[acn$f_var_peak][acn$peak_ord] %in% sup_cl

	y_of_cl = tapply(seq(0,1, l=sum(f)), ac$km$cluster[acn$f_var_peak][acn$peak_ord][f], mean)

	message("tot peaks ", sum(f))

	if(is.null(foc_tfs)) {
		foc_tfs = names(acn$tfs_top_cor)[which(acn$tfs_top_cor %in% sup_cl)]

		message("tfs ", paste(foc_tfs, collapse=" "))

		if(length(foc_tfs) <= 2) {
			top_tfs = c()
			vca = ac$vclst_cent_abs
			rownames(ac$vclst_cent_abs)
			for(cl in sup_cl) {
				add_tfs = names(tail(sort(apply(tfs_legc, 1, cor, ac$vclst_cent_abs[cl,])),2))
				foc_tfs = c(foc_tfs, add_tfs)
			}
			foc_tfs = unique(foc_tfs)
		}
		hc_tfs = hclust(tgs_dist(tfs_lfp[foc_tfs, acn$ord_clmc]),"ward.D2")
		foc_tfs = foc_tfs[hc_tfs$order]
	}

	png(file.path(ATAC.FIGURE.DIR, 'atac_clfig5.png'), 
						w = width, h=50*(4*length(sup_cl)+1+length(foc_tfs))+50)
	layout(matrix(1:3,nrow=3), h=c(4*length(sup_cl),1,length(foc_tfs)))
	par(mar = c(0,6,2,12)) 
	image(t(legc_act_v[acn$peak_ord, acn$ord_clmc][f,]), 
						col=abs_shades, breaks=c(-16.6,seq(-16,-13.8,l=999),-11), 
						xaxt='n', yaxt='n')
	mtext(sup_cl, side = 4, las=2, at=y_of_cl[as.character(sup_cl)], cex=3, line=0.5)
	par(mar=c(0,1.4,0,7.3))
	barplot(rep(1,ncol(mm$a_legc)), col=mm$mclst_col[acn$ord_clmc], 
					border=NA, xaxt='n', yaxt='n')
	par(mar=c(4,6,0,12))
	image(t(as.matrix(tfs_lfp[foc_tfs, acn$ord_clmc])), col=pos_shades, 
							breaks=c(-10,seq(-2,2,l=999),10), xaxt='n', yaxt='n')
	mtext(foc_tfs, side=4, las=2, at=seq(0,1,l=length(foc_tfs)), line=0.5, cex=2.6)
	dev.off()
		
}

order_clmc_by_tfs = function(mm, vclst_cent_abs)
{
	all_tfs = intersect(read.table(file.path(ATAC.DATA.DIR, 'all_tfs.txt'))$x, rownames(mm$r_legc))

	tfs_legc = mm$r_legc[all_tfs,]
	tfs_lfp = tfs_legc-rowMeans(tfs_legc)

	xcor = tgs_cor(t(tfs_legc), t(vclst_cent_abs))

	top_lf_tfs = apply(tfs_lfp,2,function(x) names(tail(sort(x),3)))
	top_tfs = apply(xcor,2,function(x) names(tail(sort(x),3)))
	top_tfs = unique(c(top_tfs[1,], top_tfs[2,],top_tfs[3,],top_lf_tfs[3,]))
#	top_tfs = unique(c(top_tfs[1,],top_tfs[2,], top_tfs[3,]))
	hot_xcor = xcor[top_tfs,]

	tfs_lfp = tfs_lfp[top_tfs,]
#	hc = hclust(tgs_dist(tgs_cor(tfs_lfp)),"ward.D2")
	hc = hclust(tgs_dist(t(tfs_lfp)),"ward.D2")
	ord_tf = order(apply(tfs_lfp[,hc$order],1,which.max))

	pos_shades = colorRampPalette(c("white","white", "white","pink","red","darkred", "yellow"))(1000)

	pheatmap::pheatmap(tfs_lfp[ord_tf, hc$order], cluster_cols=F, cluster_rows=F, 
								col=pos_shades, 
								file=file.path(ATAC.FIGURE.DIR, 'tfs_r_on_clmc.png'), 
								w=18,h=18)

	vclst_ord = order(apply(vclst_cent_abs[,hc$order],1,which.max))
	pheatmap::pheatmap(hot_xcor[ord_tf,vclst_ord], 
					cluster_rows=F, cluster_cols=F, 
					breaks=c(-1,seq(0,1,l=1000)), col=pos_shades, 
					file=file.path(ATAC.FIGURE.DIR, 'tfs_on_plcust.png'), w=18, h=14)

	tfs_top_cor = rownames(vclst_cent_abs)[apply(hot_xcor, 1, which.max)]
	names(tfs_top_cor) = top_tfs
	return(list(tfs_top_cor=tfs_top_cor, hc_clmc = hc, ord_clmc = hc$order))
}

main_tfs = c("Pou5f1", "Utf1", "Eomes", "Mesp1", "Elf5", "Elf4", "Satb2","Creb3l3", "Ets1", "Gata1","Hand2", "Pitx1", "Rfx4", "Sox1", "Hoxd4","Nr2f1", "Cdx2", "Foxc1", "Prrx2", "Pax7","Eya1","Foxa1", "Foxj1", "Trp63", "Grhl2", "T")

