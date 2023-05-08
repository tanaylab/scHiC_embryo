
get.pooled.mean.a.trace <- function(ab.scores, bins, should.reverse, diff.range=seq(-4e5, 4e5, 4e4)) {
  a.trace.per.bin = do.call(rbind, lapply(1:length(bins), function(i) {
    bin = bins[i]
    splitted = strsplit(bin, '_')
    chrom = splitted[[1]][1]
    coord = as.numeric(splitted[[1]][2])
    bins.in.range = paste0(chrom, '_', trimws(format(coord + diff.range, scientific=F)))
    a.scores = ab.scores[bins.in.range]
    if (should.reverse[i]) {
      a.scores = rev(a.scores)
    }
    return(a.scores)
  }))
  colnames(a.trace.per.bin) = diff.range
  return(colMeans(a.trace.per.bin, na.rm=T))
}

get.pooled.shaman.score <- function(track.conts, intervs, ab.scores.all.clusters, window.size=2e6, res=1e5, is.bin.exact=F) {

  intervs$orig_mid = as.numeric(intervs$orig_start) + as.numeric(intervs$seq_len) * 2e5 / 2
  intervs$tss = intervs$orig_mid
  track.conts$chrom1 = as.character(track.conts$chrom1)
  track.conts$chrom2 = as.character(track.conts$chrom2)
  possible.bins1 = seq(-window.size / 2, window.size / 2, res)
  possible.bins2 = seq(0, window.size / 2, res)
  all.shamans = mclapply(1:nrow(intervs), function(i) {
    print(paste0('i is ', i))
    print(paste0('out of ', nrow(intervs)))
    cur.interv = intervs[i,]
    cur.start = as.numeric(cur.interv$orig_start)
    cur.end = as.numeric(cur.interv$orig_start) + as.numeric(cur.interv$seq_len) * 2e5
    cur.chrom = as.character(cur.interv$chrom)
    if (is.bin.exact) {
      sel.bin.name = as.character(cur.interv$bin.name)
    } else {
      rel.bin.names = intersect(get.all.bins(cur.chrom, cur.start, cur.end, 4e4)[[1]], rownames(ab.scores.all.clusters))
      sel.bin.name = rel.bin.names[which.max(ab.scores.all.clusters[rel.bin.names, 'ery'] - ab.scores.all.clusters[rel.bin.names, 'emb'])]
    }
    coord = as.numeric(strsplit(sel.bin.name, '_')[[1]][2])
    start1.shifted = track.conts$start1 - coord
    start2.shifted = track.conts$start2 - coord
    # start1.new is the col, start2.new is the row
    start1.new = start1.shifted + start2.shifted
    start2.new = start1.shifted - start2.shifted
    if (cur.interv$should_reverse) {
      start1.new = -start1.new
    }
    track.conts$start1.new = start1.new
    track.conts$start2.new = start2.new
    conts.in.window = filter(track.conts, chrom1 == cur.interv$chrom & chrom2 == cur.interv$chrom & 
                      between(start1.new, -window.size / 2, window.size / 2) & 
		      between(start2.new, 0, window.size / 2))

    if (nrow(conts.in.window) == 0) {
      ret.mat = matrix(NA, nrow=length(possible.bins2), ncol=length(possible.bins1))
      rownames(ret.mat) = possible.bins2
      colnames(ret.mat) = possible.bins1
      return(ret.mat)
    }

    conts.in.window$bin1 = factor(floor(conts.in.window$start1.new / res) * res, levels=as.numeric(possible.bins1))
    conts.in.window$bin2 = factor(floor(conts.in.window$start2.new / res) * res, levels=as.numeric(possible.bins2))
    conts.in.window$var = 'mean_shaman'
    mean.shaman = conts.in.window %>% group_by(bin1, bin2) %>% summarise(mean_shaman=mean(score)) %>% ungroup() %>% complete(bin1, bin2)
    mean.shaman.mat = as.matrix(spread(mean.shaman, bin1, mean_shaman)) # for some reason this line cannot be parallelized...
    rownames(mean.shaman.mat) = mean.shaman.mat[,1]
    mean.shaman.mat = mean.shaman.mat[,-1]
    mean.shaman.mat.numeric = matrix(as.numeric(mean.shaman.mat), nrow=nrow(mean.shaman.mat))
    rownames(mean.shaman.mat.numeric) = rownames(mean.shaman.mat)
    colnames(mean.shaman.mat.numeric) = colnames(mean.shaman.mat)
    return(mean.shaman.mat.numeric)
  }, mc.cores=1)
  mat.sum = Reduce(`+`, lapply(all.shamans, function(x) {x[is.na(x)] = 0;x}))
  mat.na.sum = Reduce(`+`, lapply(all.shamans, function(x) {!is.na(x)}))
  shaman.mean = mat.sum / mat.na.sum
  return(shaman.mean)
}

plot.pooled.shaman <- function(emb.ery.tracks.conts, ery.clust, fig.dir=get.fig.dir(), repl.data.dir=get.data.file.dir()) {
  
  ab.scores = ab_scores_c_clusters
  bin.clust = ery.clust$clusters
  should.reverse = ery.clust$should_reverse
  bins.per.clust = split(names(bin.clust), bin.clust)
  data.dir = file.path(repl.data.dir, 'pooled_shamans')
  dir.create(data.dir, showWarnings=F)
  for (j in 2:1) {
    cur.conts = emb.ery.tracks.conts[[j]]
    cur.track.name = list('emb', 'ery')[[j]]
    line.col = list('black', 'red')[[j]]
    cur.ab.scores = ab.scores[, cur.track.name]
    names(cur.ab.scores) = rownames(ab.scores)
    lapply(seq_along(bins.per.clust), function(i) {
      cur.bins = bins.per.clust[[i]]
      coord = sapply(strsplit(cur.bins, '_'), function(x) as.numeric(x[2]))
      chrom = sapply(strsplit(cur.bins, '_'), function(x) x[1])
      intervs = data.frame(start=coord, end=coord, bin.name=cur.bins, seq_len=0, orig_start=coord, chrom=chrom, should_reverse=should.reverse[cur.bins])
      cur.file.path = sprintf('%s_%d.png', cur.track.name, i)
      pooled.shaman.score = get.or.create(file.path(data.dir, cur.file.path), function() get.pooled.shaman.score(cur.conts, intervs, NULL, window.size=2e6, res=2e4, is.bin.exact=T))
      pooled.a.score = get.pooled.mean.a.trace(cur.ab.scores, cur.bins, should.reverse[cur.bins], diff.range=seq(-1e6, 1e6, 4e4))

      png(file.path(fig.dir, cur.file.path))
      par(mar=c(2, 2, 2, 2))
      shades = shaman_score_pal()
      layout.mat = matrix(1, ncol=4, nrow=4)
      layout.mat = rbind(layout.mat, rep(2, 4))
      layout(layout.mat)
      image(t(pooled.shaman.score), col=shades, xaxt='n', yaxt='n')
      plot(as.numeric(names(pooled.a.score)), pooled.a.score, type='l', lwd=5, ylim=c(0, 1), col=line.col, xaxs='i', yaxs='i')
      dev.off()

    })
  }
}
