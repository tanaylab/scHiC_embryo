
schic_init_env = function()
{
	decay_metrics = sch_decay_metrics
	esc_decay_metrics = esc_sch_decay_metrics
	data.dir = get.data.file.dir()
	# precomputed TADs
	tads = as.data.frame(fread(file.path(data.dir, "sch_tads.txt"), sep="\t", stringsAsFactors=F, h=T))
	interv_tad_ab = as.data.frame(fread(file.path(data.dir, "/sch_interv_tad_ab.txt"), sep="\t", stringsAsFactors=F, h=T))
	load(file.path(data.dir, "good_cells.Rda"))

	rownames(decay_metrics) = decay_metrics$cell
	rownames(esc_decay_metrics) = esc_decay_metrics$cell
	rownames(tads) = tads$id

	env = list(decay_metrics = decay_metrics,
				  esc_decay_metrics = esc_decay_metrics,
				  tads = tads,
				  interv_tad_ab = interv_tad_ab,
			     emb_cells = good_cells)
	return(env)
}

schic_init_perlim_cell_groups = function(env)
{
	# use early_f quantiles from Yaniv's paper
	early_f_q95  <- 0.534
	early_f_q99  <- 0.565
	early_f_q995 <- 0.572
	emb_cells = env$emb_cells
	# Group g1 cells as either locked are control based on early_f fraction
	env$g1_cells = env$decay_metrics[emb_cells,] %>%
		  filter(group==2 & early_f < early_f_q95) %>%
		  pull(cell)
	env$locked_cells = env$decay_metrics[emb_cells,] %>%
				filter(group==2 & early_f > early_f_q995) %>%
				pull(cell)
	env$repl_cells = env$decay_metrics[setdiff(env$emb_cells,env$locked_cells),] %>%
				 filter(early_f > early_f_q995) %>%
				 pull(cell)

	env$g1_cells = setdiff(env$g1_cells, env$locked_cells)

	env$esc_cells = env$esc_decay_metrics %>%
				 pull(cell)
	env$esc_repl_cells = env$esc_decay_metrics %>%
				 filter(early_f > early_f_q995) %>%
				 pull(cell)
	return(env)
}

schic_init_tss_intervs = function(env)
{
	tss = gintervals.load("intervs.global.tss")
	tss_intervs = data.frame(chrom=tss$chrom, 
						start = tss$start-25e+3, 
						end = tss$start+25e+3, 
						gene = tss$geneSymbol)
	tss_intervs_r = gintervals.force_range(intervals=tss_intervs)
	tss_intervs = tss_intervs_r[!duplicated(tss_intervs$gene),]
	rownames(tss_intervs) = tss_intervs$gene

	tss_intervs = data.frame(chrom=tss$chrom, 
						start = tss$start-100e+3, 
						end = tss$start+100e+3, 
						gene = tss$geneSymbol)
	tss_intervs_r = gintervals.force_range(intervals=tss_intervs)
	tss_intervs = tss_intervs_r[!duplicated(tss_intervs$gene),]
	rownames(tss_intervs) = tss_intervs$gene
	env$tss_intervs_25k = tss_intervs
	env$tss_intervs = tss_intervs
	return(env)
}
