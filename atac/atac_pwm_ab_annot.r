
library("misha")
library("metacell")
library("ggplot2")
library("ggrepel")

project_interv_mm9 = function(mm, ac, acn)
{
	gsetroot(SCHIC.MISHA.PATH)
	intervs = mm$atac_intervs
	# TODO: where did this come from?
	imm9 = gintervals.liftover(intervs, chain=file.path(ATAC.DATA.DIR, 'mm10ToMm9.over.chain.fixed1'))

	imm9_300 = imm9
	x=round((imm9_300$start+imm9_300$end)/2)
	imm9_300$start = x-150
	imm9_300$end = x+150
	imm9_300 = gintervals.force_range(imm9_300)
	acn$interv_mm9 = imm9_300
	return(acn)
}

project_ab_on_atac_interv = function(mm, ac, acn)
{
	km = ac$km
	imm9_300 = acn$interv_mm9
	# TODO: where did this come from?
	ab_per_bin40 = read.table(file.path(ATAC.DATA.DIR, 'ab_per_bin40.txt'), h=T)

	ab_annots_mm9 = gintervals.neighbors(imm9_300, ab_per_bin40)
	ab_annots = mm$atac_intervs
	ab_annots$C21 = rep(NA, nrow(ab_annots))
	ab_annots$C21[ab_annots_mm9$intervalID] = ab_annots_mm9$C21
	ab_annots$C22 = rep(NA, nrow(ab_annots))
	ab_annots$C22[ab_annots_mm9$intervalID] = ab_annots_mm9$C22
	ab_annots$C23 = rep(NA, nrow(ab_annots))
	ab_annots$C23[ab_annots_mm9$intervalID] = ab_annots_mm9$C23
	ab_annots$ery = rep(NA, nrow(ab_annots))
	ab_annots$ery[ab_annots_mm9$intervalID] = ab_annots_mm9$ery
	ab_annots$esc = rep(NA, nrow(ab_annots))
	ab_annots$esc[ab_annots_mm9$intervalID] = ab_annots_mm9$esc
	ab_annots$emb = rep(NA, nrow(ab_annots))
	ab_annots$emb[ab_annots_mm9$intervalID] = ab_annots_mm9$emb

	ab_annots$km_cls = km$cluster

	acn$ab_annots = ab_annots
	return(acn)
}

add_tad_to_atac_interv = function(mm, ac, acn)
{

	# TODO: where did this come from?
	tads = as.data.frame(fread(file.path(ATAC.DATA.DIR, 'sch_tads.txt'), sep="\t", stringsAsFactors=F, h=T))

	imm9_tad = gintervals.neighbors(acn$interv_mm9, tads, maxdist=0)

	acn$ab_annots$tad_id = rep(NA, nrow(acn$ab_annots))
	acn$ab_annots$tad_id[imm9_tad$intervalID] = imm9_tad$id
	acn$ab_annots$tad_len = rep(NA, nrow(acn$ab_annots))
	acn$ab_annots$tad_len[imm9_tad$intervalID] = imm9_tad$len
	return(acn)
}

plot_vclst_ab_scatters = function(mm, ac, acn)
{
	
	embesc = tapply((acn$ab_annots$emb-acn$ab_annots$esc) > 0.1, 
											acn$ab_annots$km_cls, mean, na.rm=T)

	escemb = tapply((acn$ab_annots$emb-acn$ab_annots$esc) < -0.1, 
											acn$ab_annots$km_cls, mean, na.rm=T)

	eryemb = tapply((acn$ab_annots$ery-acn$ab_annots$emb) > 0.1, 
											acn$ab_annots$km_cls, mean, na.rm=T)

	embery = tapply((acn$ab_annots$ery-acn$ab_annots$emb) < -0.1, 
											acn$ab_annots$km_cls, mean, na.rm=T)

	ectomeso= tapply((acn$ab_annots$C21-acn$ab_annots$C23) > 0.1, 
											acn$ab_annots$km_cls, mean, na.rm=T)

	mesoecto = tapply((acn$ab_annots$C21-acn$ab_annots$C23) < -0.1, 
											acn$ab_annots$km_cls, mean, na.rm=T)

	nms = ac$vclst_nms

	df= data.frame(emb=embesc[nms], esc=escemb[nms], col=ifelse(escemb[nms]>0.08,"black","gray"), lab=nms )
	lab_theme = element_text(size=24)
	force=30
	f = embesc[nms]>0.04 | escemb[nms]>0.04
   p = ggplot(df, aes(emb, esc)) +
            geom_point(col=df$col, size=5) +
            labs(x="A: Embryo > ESC") +
            labs(y="A: ESC > Embryo") +
            theme(axis.text.x = lab_theme, axis.text.y = lab_theme,
                  axis.title.x = lab_theme, axis.title.y = lab_theme) +
            geom_text_repel(data = df[f,],
                     size=5,force=force, max.iter=1e+4,aes(label=lab))
	pdf(file.path(ATAC.FIGURE.DIR, 'AB_embesc.pdf'), h=10, w=10)
	print(p)
	dev.off()
	
	df= data.frame(emb=embery[nms], ery=eryemb[nms], col=ifelse(eryemb[nms]>0.25,"red","gray"), lab=nms )
	lab_theme = element_text(size=24)
	force=30
	f = embery[nms]>0.08 | eryemb[nms]>0.08
   p = ggplot(df, aes(emb, ery)) +
            geom_point(col=df$col, size=5) +
            labs(x="A: Embryo > Ery") +
            labs(y="A: Ery > Embryo") +
            theme(axis.text.x = lab_theme, axis.text.y = lab_theme,
                  axis.title.x = lab_theme, axis.title.y = lab_theme) +
            geom_text_repel(data = df[f,],
                     size=5,force=force, max.iter=1e+4,aes(label=lab))
	pdf(file.path(ATAC.FIGURE.DIR, 'AB_embery.pdf'), h=10, w=10)
	print(p)
	dev.off()

	df= data.frame(ecto=ectomeso[nms], meso=mesoecto[nms], col=ifelse(ectomeso[nms]>0.1,"darkgreen", ifelse(mesoecto[nms]>0.1,"orange","gray")), lab=nms )
	lab_theme = element_text(size=24)
	force=30
	f = ectomeso[nms]>0.04 | mesoecto[nms]>0.04
   p = ggplot(df, aes(ecto, meso)) +
            #geom_point(col=df$col, size=7) +
            geom_point(col=df$col, size=5) +
            labs(x="A: Ectoderm > Mesoderm") +
            labs(y="A: Mesoderm > Ectoderm") +
            theme(axis.text.x = lab_theme, axis.text.y = lab_theme,
                  axis.title.x = lab_theme, axis.title.y = lab_theme) +
            geom_text_repel(data = df[f,],
                     size=5,force=force, max.iter=1e+4,aes(label=lab))
            #geom_text_repel(data = df[f,],
            #         size=8,force=force, max.iter=1e+4,aes(label=lab))
	pdf(file.path(ATAC.FIGURE.DIR, 'AB_ectomeso.pdf'), h=10, w=10)
	print(p)
	dev.off()


}

plot_peak_proxim_mat = function(mm, ac, acn)
{
	flags = acn$ab_annots[,c(1,2,3,5,12,13,14)]
	a = gintervals.neighbors(flags, flags, maxdist=1e+6,mindist=1e+3, maxneighbors=1000)

	f = a$dist < 200000 & a[,6]==a[,13]
	prox = table(a[f,5], a[f,12])
	cnt = table(a[f,5])
	e_prox = (cnt %*% t(cnt))/sum(cnt)
	prox_n = (1+prox)/(1+e_prox)
	clst_cent_norm = ac$clst_cent_abs-rowMeans(ac$clst_cent_abs)

	shades = colorRampPalette(c("black","white","darkred"))(1000)
	vshades = colorRampPalette(c("white","white", "darkblue", "brown","yellow"))(1000)
	
#split abs/var and order it
	vnms = ac$vclst_nms
   f_big_km = table(ac$km$cluster)>99
	cnms = setdiff(colnames(prox_n)[f_big_km],ac$vclst_nms)

	prox_n_v = prox_n[vnms, vnms]
	v_hc = hclust(tgs_dist(tgs_cor(log2(prox_n_v))),"ward.D2")
	prox_n_c = prox_n[cnms, cnms]
	c_hc = hclust(tgs_dist(tgs_cor(log2(prox_n_c))),"ward.D2")

	prox_ord = c(rev(cnms[c_hc$order]), vnms[v_hc$order])

	pheatmap::pheatmap(log2(prox_n)[prox_ord, prox_ord], 
							cluster_rows=F, cluster_cols=F, col=shades, 
							breaks=seq(-2,2,l=1001), file=file.path(ATAC.FIGURE.DIR, 'prox_clst.png'), w=20,h=20)

	
	clmc_clst = cutree(acn$hc_clmc, 17)
	clmc_clst_map = rank(tapply(1:length(acn$ord_clmc), 
												clmc_clst[acn$ord_clmc], mean))
	clmc_clst = clmc_clst_map[clmc_clst]
	clmc_clst_cent_norm = t(tgs_matrix_tapply(clst_cent_norm, clmc_clst, mean))
	clmc_clst_cent = t(tgs_matrix_tapply(ac$clst_cent_abs, clmc_clst, mean))

	abs_shades = colorRampPalette(c("white", "white", "lightblue","blue","darkblue", "brown"))(1000)
	
	# I (Nimrod) added the next lines:
	rownames(clmc_clst_cent_norm) = 1:nrow(clmc_clst_cent_norm)
	rownames(clmc_clst_cent) = 1:nrow(clmc_clst_cent)

	pheatmap::pheatmap(clmc_clst_cent_norm[prox_ord, ], 
							cluster_rows=F, cluster_cols=F, col=vshades, 
							breaks=seq(-1,3,l=1001), 
							file=file.path(ATAC.FIGURE.DIR, 'norm_agg_by_prox_ord.png'), w=20,h=20)

	pheatmap::pheatmap(clmc_clst_cent[prox_ord, ], 
							cluster_rows=F, cluster_cols=F, col=vshades, 
							breaks=seq(-17, -13, l=1001), 
							file=file.path(ATAC.FIGURE.DIR, 'agg_by_prox_ord.png'), w=20,h=20)

	pheatmap::pheatmap(clst_cent_norm[prox_ord, acn$ord_clmc], 
							cluster_rows=F, cluster_cols=F, col=vshades, 
							breaks=seq(-1,3,l=1001), 
							file=file.path(ATAC.FIGURE.DIR, 'norm_by_prox_ord.png'), w=20,h=20)

	pheatmap::pheatmap(ac$clst_cent_abs[prox_ord, acn$ord_clmc], 
							cluster_rows=F, cluster_cols=F, col=abs_shades, 
							breaks=c(seq(-17, -14.5, l=1000),-10), 
							file=file.path(ATAC.FIGURE.DIR, 'abs_by_prox_ord.png'), w=20,h=20)

}
