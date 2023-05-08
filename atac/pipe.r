library(metacell)

#this takes time
#it write the multi model object 
prepare_reik_multimodl()

load(file.path(ATAC.DATA.DIR, "multi_reik.Rda"), v=T); 
mm = multi_model

ac = clust_atac_bins(mm, k_n=120, T_vclst_min=-16)

rownames(mm$r_legc) = rownames(mm$r_all_legc)

acn = plot_atac_clusts(mm, ac, k_sup=14)

acn = project_interv_mm9(mm, ac, acn)
acn = project_ab_on_atac_interv(mm, ac, acn)
acn = add_tad_to_atac_interv(mm, ac, acn)

#fig 5A lower panel (cropped)
plot_atac_supclst_tfs(mm, ac, acn, sup_cl = c(38,37,27,9,8,5,117,99,95,98,75,69,73,76),
				   foc_tfs=NULL, width=1200) 

tfs = c("Pou5f1", "Utf1", "Eomes", "Mesp1", "Elf5", "Elf4", "Satb2", "Creb3l3", "Ets1", "Gata1", "Hand2","Pitx1","Rfx4", "Sox1", "Hoxd4", "Nr2f1", "Cdx2", "Foxc1", "Prrx2", "Pax7","Eya1", "Foxa1", "Foxj1", "Trp63", "Grhl2", "T")
r_lfp = mm$r_legc-rowMeans(mm$r_legc)

shades = colorRampPalette(c("white","white", "pink","darkred", "black", "black"))(1000)
pheatmap::pheatmap(r_lfp[tfs, acn$ord_clmc], 
							cluster_rows=F, cluster_cols=F, col=shades,
							filename = file.path(ATAC.FIGURE.DIR, 'fig5a_top.png'), width=12, height=4)

#fig 5B,C,D
plot_vclst_ab_scatters(mm, ac, acn)

#fig 5E
gsetroot(ATAC.MISHA.PATH)
plot_peak_proxim_mat(mm, ac, acn)

