
# before running this schicluster should be installed in python,
# and that python should be on the path
# if you already have the output file.output.path file, you can just load it
# and skip running the rest
export.data.for.schicluster <- function() {
  # we'll start with the embryo clusters
  repl.data.dir = get.data.file.dir()
  embryo.cls.path = file.path(repl.data.dir, 'embryo_cls')
  embryo.cls.ret = get.or.create(embryo.cls.path, cluster.embryo.cells)
  model.params = embryo.cls.ret$model.params
  clusters = model.to.clusters(model.params)
  selected.cells = names(clusters)

  #load.data()
  #selected.cells = colnames(bin_cell_cov)
  #decay.metrics = get.common.decay.metrics()
  #selected.cells = colnames(bin_cell_cov)[grep('elad', colnames(bin_cell_cov))]
  #selected.cells = sample(selected.cells, 100, replace=F)

  schicluster.input.dir = file.path(SCHICLUSTER.DIR, 'input/')
  schicluster.imputed.dir = file.path(SCHICLUSTER.DIR, 'imputed/')
  schicluster.embedded.dir = file.path(SCHICLUSTER.DIR, 'embedded/')

  # use 1Mbp resolution, as they do in their paper
  #for (i in seq_along(selected.cells)) {
  schicluster.res = 1e6
  #for (i in 1:2) {
  mclapply(seq_along(selected.cells), function(i) {
    print(i)
    cur.cell = selected.cells[i]
    #all.output.path.prefixes = sprintf('%s_%s_.*npz', cur.cell, paste0('chr', 1:19))
    all.output.path.prefixes = sprintf('%s_%s_', cur.cell, paste0('chr', 1:19))
    cur.files.in.dir = list.files(schicluster.imputed.dir)
    if (all(sapply(all.output.path.prefixes, function(cur.prefix) length(grep(cur.prefix, cur.files.in.dir))) > 0)) {
      return(NULL)
    }
    vtrack.name = paste0(cur.cell, '_area')
    gvtrack.create(vtrack.name, cur.cell, "area")
    contacts = gextract(vtrack.name, intervals=gintervals.2d.all(), iterator=c(schicluster.res, schicluster.res), colnames='count')
    cis.contacts = filter(contacts, chrom1 == chrom2)
    for (cur.chrom in paste0('chr', 1:19)) {
      cur.cis.contacts = filter(cis.contacts, chrom1 == cur.chrom & start1 <= start2 & count > 0)
      start.indices = trimws(format(floor(cur.cis.contacts$start1 / 1e6), scientific=F))
      end.indices = trimws(format(floor(cur.cis.contacts$start2 / 1e6), scientific=F))
      cur.df = data.frame(start.indices, end.indices, count=cur.cis.contacts$count, stringsAsFactors=F)
      cur.df[cur.df$start.indices == cur.df$end.indices, 'count'] = cur.df[cur.df$start.indices == cur.df$end.indices, 'count'] / 2
      cur.file.path = file.path(schicluster.input.dir, sprintf('%s_%s.txt', cur.cell, cur.chrom))
      write.table(cur.df, file=cur.file.path, sep='\t', quote=F, row.names=F, col.names=F) 
      cur.command = sprintf('hicluster impute-cell --indir %s --outdir %s --cell %s --chrom %s --res 1000000 --chrom_file %s', 
                            schicluster.input.dir, schicluster.imputed.dir, cur.cell, cur.chrom, file.path(SCHICLUSTER.DIR, 'chrom_sizes.txt'))
      system(cur.command)
    }
  }, mc.cores=25)
  for (cur.chrom in paste0('chr', 1:19)) {
    print(paste0('cur chrom ', cur.chrom))
    chrom.fnames = file.path(schicluster.imputed.dir, sapply(sprintf('%s_%s_.*hdf5', selected.cells, cur.chrom), function(cur.re) {
    #chrom.fnames = file.path(schicluster.imputed.dir, sapply(sprintf('%s_%s_.*npz', selected.cells, cur.chrom), function(cur.re) {
      cur.fname = grep(cur.re, list.files(schicluster.imputed.dir), value=T)
      #if (length(cur.fname) != 1) browser()
      stopifnot(length(cur.fname) == 1)
      return(cur.fname)
    }))
    imputed.file.list.path = file.path(SCHICLUSTER.DIR, paste0(cur.chrom, '_imputed_file_list'))
    write.table(chrom.fnames, quote=F, row.names=F, col.names=F, file=imputed.file.list.path)

    #cur.output.path = file.path(schicluster.embedded.dir, paste0('pad1_std1_rp0.5_sqrtvc_', cur.chrom))
    cur.output.path = file.path(schicluster.embedded.dir, cur.chrom)
    cur.command = sprintf('hicluster embed-concatcell-chr --cell_list %s --outprefix %s --res 1000000',
                          imputed.file.list.path, cur.output.path)
    
    if (!file.exists(file.path(cur.output.path))) {
      system(cur.command)
    }
    #system(cur.command)
  }
  #all.output.paths = file.path(schicluster.embedded.dir, paste0('pad1_std1_rp0.5_sqrtvc_chr', 1:19, '.npz'))
  all.output.paths = file.path(schicluster.embedded.dir, paste0('chr', 1:19, '.svd50.npy'))
  embedded.file.list.path = file.path(SCHICLUSTER.DIR, 'embedded_file_list')
  final.output.path = file.path(SCHICLUSTER.DIR, 'output')
  write.table(all.output.paths, quote=F, row.names=F, col.names=F, file=embedded.file.list.path)
  final.command = sprintf('hicluster embed-mergechr --embed_list %s --outprefix %s',
                          embedded.file.list.path, final.output.path)
  #if (!file.exists(final.output.path)) {
  #  system(final.command)
  #}
  system(final.command)


  schicluster.output = anndata::read_hdf(paste0(final.output.path, '.svd50.hdf5'), 'data')
  umap.config = umap.defaults
  umap.config$random_state = 42
  umap.ret = umap(data.matrix(schicluster.output), umap.config)
  #cls.col = MODEL.CLUSTER.COLORS[-2]
  fig.dir = file.path(get.fig.dir(), 'repl_model')
  png(file.path(fig.dir, 'schicluster_umap_model_cells.png'))
  plot(umap.ret$layout[,1], umap.ret$layout[,2], bg=MODEL.CLUSTER.COLORS[clusters], pch=21, cex=1.2)
  dev.off()

}


# This runs the schic topic model
# we executed the pre and post processing in R 3.5,
# the the actual topic model in R 3.6, because the topic model didn't support R 3.5.
export.data.for.schic.topic.model <- function() {
  # we'll start with the embryo clusters
  repl.data.dir = get.data.file.dir()
  embryo.cls.path = file.path(repl.data.dir, 'embryo_cls')
  embryo.cls.ret = get.or.create(embryo.cls.path, cluster.embryo.cells)
  model.params = embryo.cls.ret$model.params
  clusters = model.to.clusters(model.params)
  erys = get.erys()
  #selected.cells = setdiff(colnames(bin_cell_cov)[grep('elad', colnames(bin_cell_cov))], erys)
  #selected.cells = names(clusters)
  selected.cells = colnames(bin_cell_cov)

  #load.data()
  #decay.metrics = get.common.decay.metrics()
  #selected.cells = colnames(bin_cell_cov)[grep('elad', colnames(bin_cell_cov))]
  #selected.cells = sample(selected.cells, 100, replace=F)

  schic.topic.dir = SCHIC.TOPIC.MODEL.DIR
  chrom.sizes.df = read.table(file.path(SCHICLUSTER.DIR, 'chrom_sizes.txt'), stringsAsFactors=F)
  chrom.sizes = chrom.sizes.df[,2]
  names(chrom.sizes) = chrom.sizes.df[,1]

  # use 0.5Mbp resolution, as they do in their paper
  topic.res = 5e5
  all.bins = c()
  for (i in seq_along(chrom.sizes)) {
    bin.starts = seq(1, chrom.sizes[i], by=topic.res)
    cur.bin.names = sprintf('%s_%s', names(chrom.sizes)[i], trimws(format(bin.starts, scientific=F)))
    cur.bin.indices = (length(all.bins) + 1):(length(all.bins) + length(cur.bin.names))
    names(cur.bin.indices) = cur.bin.names
    all.bins = c(all.bins, cur.bin.indices)
  }

  mclapply(seq_along(selected.cells), function(i) {
    print(i)
    cur.cell = selected.cells[i]
    output.fname = file.path(schic.topic.dir, cur.cell)
    if (file.exists(output.fname)) {
      return(NULL)
    }

    vtrack.name = paste0(cur.cell, '_area')
    gvtrack.create(vtrack.name, cur.cell, "area")
    contacts = gextract(vtrack.name, intervals=gintervals.2d.all(), iterator=c(topic.res, topic.res), colnames='count')
    cis.contacts = filter(contacts, chrom1 == chrom2 & start1 <= start2 & chrom1 != 'chrX' & chrom1 != 'chrY' & chrom1 != 'chrM' & count > 0)
    bin.names1 = sprintf('%s_%s', cis.contacts$chrom1, trimws(format(cis.contacts$start1 + 1, scientific=F)))
    bin.names2 = sprintf('%s_%s', cis.contacts$chrom2, trimws(format(cis.contacts$start2 + 1, scientific=F)))
    cis.contacts$indices1 = all.bins[bin.names1]
    cis.contacts$indices2 = all.bins[bin.names2]
    # 19 is for 20 bins of 5e5 = 10Mbp
    cis.contacts.fil = cis.contacts[(cis.contacts$indices2 - cis.contacts$indices1) <= 19,]
    cur.df = data.frame(bin1=cis.contacts.fil$indices1, bin2=cis.contacts.fil$indices2, count=cis.contacts.fil$count)
    cur.df[cur.df$bin1 == cur.df$bin2, 'count'] = cur.df[cur.df$bin1 == cur.df$bin2, 'count'] / 2
    write.table(cur.df, quote=F, row.names=F, col.names=F, sep='\t', file=output.fname)
  }, mc.cores=25)
  all.cell.df = do.call(rbind, lapply(seq_along(selected.cells), function(i) {
    print(i)
    cur.cell = selected.cells[i]
    output.fname = file.path(schic.topic.dir, cur.cell)
    cur.df = read.table(output.fname, stringsAsFactors=F)
    colnames(cur.df) = c('bin1', 'bin2', 'count')
    cur.df$cell = i
    cur.df = cur.df[,c(4, 1:3)]
    return(cur.df)
  }))
  locus.pair = sprintf('%s_%s', all.cell.df$bin1, all.cell.df$bin2)
  all.cell.df$locus.pair = as.numeric(as.factor(locus.pair))
  df.to.write = all.cell.df[,c('cell', 'locus.pair', 'count')]
  #write.table(df.to.write, quote=F, col.names=F, row.names=F, sep='\t', file=file.path(schic.topic.dir, 'all_concat.tsv'))
  #write.table(df.to.write, quote=F, col.names=F, row.names=F, sep='\t', file=file.path(schic.topic.dir, 'all_concat_all_cells.tsv'))
  write.table(df.to.write, quote=F, col.names=F, row.names=F, sep='\t', file=file.path(schic.topic.dir, 'all_concat_embryo_cells.tsv'))

  # The next bit is run in R-3.6 for installation reasons
  #mat_file = file.path(schic.topic.dir, 'all_concat.tsv')
  #mat_file = file.path(schic.topic.dir, 'all_concat_all_cells.tsv')
  mat_file = file.path(schic.topic.dir, 'all_concat_embryo_cells.tsv')
  df = fread(
      mat_file,
      col.names = c("cell.idx", "lp.idx", "count"),
      colClasses = c("integer", "integer", "integer"))
  df$lp.idx = df$lp.idx + 1
  mat = sparseMatrix(i = df$lp.idx, j = df$cell.idx, x = df$count)
  colnames(mat) = 1:ncol(mat)
  rownames(mat) = sprintf('chr1:%s-%s', 100 * (1:nrow(mat)), 100 * (1:nrow(mat)) + 1)
  cisTopicObject <- createcisTopicObject(mat, project.name="full_cisTopic", keepCountsMatrix=FALSE)
  topic_list = 30
  cisTopicObject <- runCGSModels(cisTopicObject, topic=topic_list, seed=999, nCores=length(topic_list))
  alpha = 0
  tmp = (cisTopicObject@models[[1]])
  topic_cell_mat <- apply(tmp$document_sums, 2, function(x) {(x + alpha)/sum(x + alpha)})
  #save(topic_cell_mat, file=file.path(schic.topic.dir, 'topic_output'))
  #save(topic_cell_mat, file=file.path(schic.topic.dir, 'topic_output_all_cells'))
  save(topic_cell_mat, file=file.path(schic.topic.dir, 'topic_output_embryo_cells'))

  # and now back to R-3.5
  #load(file.path(schic.topic.dir, 'topic_output_all_cells'))
  #load(file.path(schic.topic.dir, 'topic_output_embryo_cells'))
  load(file.path(schic.topic.dir, 'topic_output'))

  umap.config = umap.defaults
  umap.config$random_state = 42
  topic.umap.ret = umap(t(topic_cell_mat), umap.config)
  #plot(topic.umap.ret$layout[,1], topic.umap.ret$layout[,2], col=clusters)

  fig.dir = get.fig.dir()
  png(file.path(fig.dir, 'topic_umap_model_cells.png'))
  plot(topic.umap.ret$layout[,1], topic.umap.ret$layout[,2], bg=MODEL.CLUSTER.COLORS[clusters], pch=21, cex=1.2)
  dev.off()



}
