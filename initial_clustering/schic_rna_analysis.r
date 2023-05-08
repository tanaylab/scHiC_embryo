
schic_rna_gen_abs_rna = function(env, bin_size=2e+5, emb_mc = "e9_orig_bs500f")
{
	library("metacell")
	scdb_init(EXPRESSION.MC.DIR, force_reinit=T)
	mc_emb = scdb_mc(emb_mc)
	mc_esc = scdb_mc("eseb_wt0_bs500f")

	emb_multi_nms = rownames(mc_emb@e_gc)
	emb_nms = strsplit(emb_multi_nms, ";")
	emb_n = unlist(lapply(emb_nms, length))
	emb_v = rep(emb_multi_nms,times=emb_n)
	names(emb_v) = unlist(emb_nms)

	esc_multi_nms = rownames(mc_esc@e_gc)
	esc_nms = strsplit(esc_multi_nms, ";")
	esc_n = unlist(lapply(esc_nms, length))
	esc_v = rep(esc_multi_nms,times=esc_n)
	names(esc_v) = unlist(esc_nms)

	ref_multi_nms = rownames(env$tss_intervs)
	ref_nms = strsplit(ref_multi_nms, ";")
	ref_n = unlist(lapply(ref_nms, length))
	ref_v = rep(ref_multi_nms,times=ref_n)
	names(ref_v) = unlist(ref_nms)

	nm_chrom = env$tss_intervs[ref_v,"chrom"]
	names(nm_chrom) = names(ref_v)
	nm_start = (env$tss_intervs[ref_v,"start"]+env$tss_intervs[ref_v,"end"])/2
	names(nm_start) = names(ref_v)

	emb_nms_chrom = nm_chrom[names(emb_v)]
	emb_nms_start = nm_start[names(emb_v)]
	emb_multi_nms_start = tapply(emb_nms_start, emb_v, mean)[emb_multi_nms]
	emb_multi_nms_chrom = tapply(as.character(emb_nms_chrom), emb_v, function(x) x[1])[emb_multi_nms]
	emb_bin = paste(emb_multi_nms_chrom, sprintf("%d", 1+floor(emb_multi_nms_start/bin_size)*bin_size), sep="_")

	emb_bin_max = apply(mc_emb@e_gc,2,function(x) tapply(x, emb_bin, max))
	
	esc_nms_chrom = nm_chrom[names(esc_v)]
	esc_nms_start = nm_start[names(esc_v)]
	esc_multi_nms_start = tapply(esc_nms_start, esc_v, mean)[esc_multi_nms]
	esc_multi_nms_chrom = tapply(as.character(esc_nms_chrom), esc_v, function(x) x[1])[esc_multi_nms]
	esc_bin = paste(esc_multi_nms_chrom, sprintf("%d", 1+floor(esc_multi_nms_start/bin_size)*bin_size), sep="_")

	esc_bin_max = apply(mc_esc@e_gc,2,function(x) tapply(x, esc_bin, max))

	f_is_ery = log2(mc_emb@mc_fp["Hba-x",])>4
	ery_bin_max = emb_bin_max[,f_is_ery]
	emb_bin_max = emb_bin_max[,!f_is_ery]

	rna = list(esc_bin_max = esc_bin_max,
                   emb_bin_max = emb_bin_max, 
                   ery_bin_max = ery_bin_max)
	return(rna)
}
