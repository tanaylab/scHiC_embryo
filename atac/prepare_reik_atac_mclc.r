
library("misha")

prepare_reik_multimodl = function()
{

	# assume that there's a track called wt_reik.marginal with the coverage to detect peaks
	gsetroot(ATAC.MISHA.PATH)
	options(gmax.data.size=5e+7)
	atac_hits = gscreen("wt_reik.marginal>300")

	hits_v = gextract(c("wt_reik.marginal"), intervals=atac_hits, iterator=20)
	hits_v = hits_v[order(-hits_v$wt_reik.marginal),]
	hits_peak = hits_v[!duplicated(hits_v$intervalID),] 
	hits_peak$start = hits_peak$start - 140
	hits_peak$end = hits_peak$end + 140
	hits_peak = gintervals.force_range(hits_peak)
	hits_peak = gintervals.canonic(hits_peak)
	hits_v = gextract(c("wt_reik.marginal"), intervals=hits_peak, iterator=20)
	hits_v = hits_v[order(-hits_v$wt_reik.marginal),]
	hits_peak = hits_v[!duplicated(hits_v$intervalID),] 
	hits_peak$start = hits_peak$start - 140
	hits_peak$end = hits_peak$end + 140

	hits_peak = hits_peak[order(hits_peak$intervalID),]

	scdb_init(EXPRESSION.MC.DIR, force_reinit=T)
	mc = scdb_mc("reik_multiome_rna_f")
	gset = scdb_gset("reik_multiome_rna_f_feats")

	load(file.path(ATAC.DATA.DIR, "rna_md.Rmd"), v=T)
	md$cell_type = as.character(md$cell_type)
	md$metacell= as.character(md$metacell)
	md$color= as.character(md$color)
	rownames(md) = md$metacell

	good_mcs = md$metacell[md$cell_type != "Doublet" & md$cell_type != "Mixture"]

	feat_legc = log2(1e-5+mc@e_gc[names(gset@gene_set),good_mcs])

	hc = hclust(tgs_dist(t(feat_legc)),"ward.D2")
	g_hc = hclust(tgs_dist(tgs_cor(t(feat_legc))),"ward.D2")

	k_clust = 300
	mclust = cutree(hc, k_clust)

	mclst_types_n = table(mclust,md[names(mclust),"cell_type"])
	mclst_type = colnames(mclst_types_n)[apply(mclst_types_n,1,which.max)]
	col_key = unique(md[,c("cell_type","color")])
	rownames(col_key) = col_key$cell_type
	mclst_color = col_key[mclst_type, "color"]


	add_mc_minic_peak_v = function(c_i)
	{
		message("processing clust ", c_i)
		mcs = names(which(mclust==c_i))
		tnms = paste("wt_reik.mc", mcs, sep="")
		tnms = intersect(tnms, gtrack.ls("wt_reik.mc"))
		peak_vs = gextract(tnms, intervals=hits_peak, iterator=hits_peak); 
		peak_vs = peak_vs[order(peak_vs$intervalID),]
		if(length(mcs)==1) {
			return(peak_vs[,4])
		} else {
			return(rowSums(peak_vs[,-c(1:3,ncol(peak_vs))], na.rm=T))
		}
	}
	if(1) {
	mc_cores =8 
	doMC::registerDoMC(mc_cores)
	peak_tot = do.call(cbind, parallel::mclapply(1:k_clust, add_mc_minic_peak_v, mc.cores = mc_cores))
	peak_tot[is.na(peak_tot)]=0
	colnames(peak_tot) = 1:k_clust
	n_calib = 10*quantile(colSums(peak_tot),0.1)
	f3 = colSums(peak_tot)>quantile(colSums(peak_tot),0.05)
	peak_tot_n = t(t(peak_tot[,f3])/colSums(peak_tot[,f3]))*n_calib
	a_legc = log2(1e-5+t(t(peak_tot[,f3])/colSums(peak_tot[,f3])))
	} 

	r_all_legc = log2(1e-5+mc@e_gc[, good_mcs])

	r_legc = t(tgs_matrix_tapply(r_all_legc, mclust, mean))
	r_legc = r_legc[,f3]

#atac_intervs - list of peaks intervals
#peak_tot - the total umi per peak per mc cluster
#a_legc, r_legc - log atac and rna values per peak/mclc
#r_all_legc - the original rna over mc
#mclst - the definition of metacell clusters
#mclst_color,type - consensus celltype /color per mccl
#filt_cov_clmc3 - mclc we filtered as low atac count
	multi_model = list(atac_intervs = hits_peak, 
						  peak_tot = peak_tot,
						  a_legc = a_legc,
						  r_legc = r_legc, 
						  r_all_legc = r_all_legc,
						  mclust = mclust, 
						  mclst_color = mclst_color[f3],
						  mclst_type = mclst_type[f3],
						  filt_cov_clmc3 = f3) 

	save(multi_model, file=file.path(ATAC.DATA.DIR, "multi_reik.Rda"))
}


