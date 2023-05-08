
schic_bin_mat = function(cell_id, intervs1, intervs2, 
									f1="id", f2="id", min_dist=1e+4, cis=T, max_dist=NULL)
{
	message("binning ", cell_id)
	if(is.null(which(f1==colnames(intervs1)))) {
		stop("Missing factor field ", f1, " in intervs1 when tryig to bin schic")
	}
	if(is.null(which(f2==colnames(intervs2)))) {
		stop("Missing factor field ", f2, " in intervs1 when tryig to bin schic")
	}

	colnames(intervs1)[colnames(intervs1)==f1] = "f1"
	colnames(intervs2)[colnames(intervs2)==f2] = "f2"

	d = gextract(cell_id, gintervals.2d.all())

	if(cis) {
		d= d %>% filter(chrom1==chrom2 & abs(start1-start2)>min_dist)
		if (!is.null(max_dist)) {
		  d= d %>% filter(chrom1==chrom2 & abs(start1-start2)<max_dist)
		}
	} else {
		d = d %>% filter(chrom1!=chrom2)
	}

	pt1 = d %>% mutate(chrom=chrom1, start=start1, end=start1+1) %>% select(chrom, start, end)
	pt2 = d %>% mutate(chrom=chrom2, start=start2, end=start2+1) %>% select(chrom, start, end)
	pt1$i = 1:nrow(pt1)
	pt2$i = 1:nrow(pt2)

	id1 = gintervals.neighbors(pt1, intervs1, maxdist=0)[,c("i","f1")]
	colnames(id1) = c("i", f1)

	id2 = gintervals.neighbors(pt2, intervs2, maxdist=0)[,c("i","f2")]
	colnames(id2) = c("i", f2)

	bins = merge(id1, id2, by=c("i"))
	
	b = bins %>% filter(f2 !='NA') %>% select(f1, f2) %>% group_by_all %>% summarise(cnt=n()) %>% mutate(cell=cell_id)

	return(b)
}

schic_bin_mat_multi = function(cells, intervs1, intervs2, 
												f1="id", f2="id",min_dist, cis=T, max_dist=NULL)
{
	if(is.null(which(f1==colnames(intervs1)))) {
		stop("Missing factor field ", f1, " in intervs1 when tryig to bin schic")
	}
	if(is.null(which(f2==colnames(intervs2)))) {
		stop("Missing factor field ", f2, " in intervs1 when tryig to bin schic")
	}
   a = do.call(rbind, mclapply(cells, schic_bin_mat, intervs1=intervs1, intervs2=intervs2, f1=f1, f2=f2, min_dist=min_dist, cis=cis, max_dist=max_dist, mc.cores=25))

	return(a)
}

flex_n_downsamp = function(umis, n)
{
	umis = umis[,colSums(umis)>= n]
	m = nrow(umis)
	.downsamp_one=function(v, replace = F) {
		n = v[1]
		v = v[-1]
		a = tabulate(sample(rep(1:length(v),times=v),replace=replace,size=n),nbins=m)
		return (a)
	}
	max_bin = tgconfig::get_param("mc_cores", "schic")
	doMC::registerDoMC(max_bin)

	max_bin = min(max_bin, ceiling(ncol(umis)/500))

	if(max_bin*10000 < ncol(umis)) {
		max_bin =  round(ncol(umis))/10000
	}
	cell_quant = ceiling(ncol(umis)/max_bin)
	seed = 19
	if(length(n) == 1) {
		n = rep(n,times=ncol(umis))
	}	
	sub_dsamp = function(x) {
		set.seed(seed)
		i = 1+(x-1)*cell_quant
		j = min(x*cell_quant, ncol(umis))
	   ret = Matrix(apply(rbind(n[i:j],umis[,i:j]), 2, .downsamp_one))
	   rownames(ret) = rownames(umis)
		return(as(ret,"dgCMatrix"))
	}
	res <- plyr::alply(1:max_bin, 1, sub_dsamp, .parallel=TRUE)
	umis_ds = do.call(cbind, res)
	return(umis_ds)
}
