
get.g1.cells <- function(decay.metrics) {
  return(rownames(decay.metrics)[decay.metrics$group == 2])
}

get.g2.cells <- function(decay.metrics) {
  return(rownames(decay.metrics)[decay.metrics$group == 4])
}

get.g1.non.ery.cells <- function(decay.metrics) {
  g1.cells = get.g1.cells(decay.metrics)
  return(setdiff(g1.cells, get.erys()))
}

get.g2.non.ery.cells <- function(decay.metrics) {
  g2.cells = get.g2.cells(decay.metrics)
  return(setdiff(g2.cells, get.erys()))
}

compute.chrom.cov <- function(clusters, total.cov, decay.metrics) {
  g1.cells.fil = intersect(colnames(total.cov), get.g1.non.ery.cells(decay.metrics))
  g1.cov = as.matrix(total.cov[, g1.cells.fil])
  bin.probs = rowSums(g1.cov) / sum(g1.cov)
  total.cov = as.matrix(total.cov[,names(clusters)])
  per.bin.cov = tgs_matrix_tapply(total.cov, clusters, sum)
  colnames(per.bin.cov) = rownames(total.cov)

  per.bin.cov.norm = t(t(per.bin.cov) / bin.probs)
  per.bin.frac = per.bin.cov.norm / rowSums(per.bin.cov.norm)
  without.per.bin.frac = t(do.call(cbind, lapply(1:nrow(per.bin.frac), function(clust) {
    cov.without.clust = colSums(per.bin.cov[-clust, , drop=F]) 
    per.bin.cov.without.norm = t(t(cov.without.clust) / bin.probs)
    return(per.bin.cov.without.norm / sum(per.bin.cov.without.norm))
  })))

  stopifnot(all(colnames(per.bin.frac) == names(bin.probs)))
  stopifnot(all(colnames(without.per.bin.frac) == names(bin.probs)))
  stopifnot(all(!is.na(per.bin.frac)))
  stopifnot(all(!is.na(without.per.bin.frac)))
  return(list(t(per.bin.frac), t(without.per.bin.frac)))
}

compute.chrom.ab <- function(clusters, ab.form, min.cov.discard=300, min.cov.per.clust=100, min.cov.cells=NULL) {
  atab.fil = ab.form$atab[names(clusters),]
  btab.fil = ab.form$btab[names(clusters),]
  if (is.null(min.cov.cells)) {
    deeply.covered.bins = colSums(atab.fil + btab.fil) > min.cov.discard
  } else {
    deeply.covered.bins = colSums(atab.fil[min.cov.cells,] + btab.fil[min.cov.cells,]) > min.cov.discard
  }
  atab.fil = atab.fil[, deeply.covered.bins]
  btab.fil = btab.fil[, deeply.covered.bins]

  clust.total.a = tgs_matrix_tapply(t(atab.fil), clusters, sum)
  colnames(clust.total.a) = colnames(atab.fil)
  clust.total.b = tgs_matrix_tapply(t(btab.fil), clusters, sum)
  colnames(clust.total.b) = colnames(btab.fil)
  a.scores = (clust.total.a) / (clust.total.a + clust.total.b)
  a.scores[clust.total.a + clust.total.b < min.cov.per.clust] = NA

  a.without.scores = do.call(rbind, lapply(1:nrow(clust.total.a), function(clust) {
    a.without.score = colSums(clust.total.a[-clust, , drop=F]) / colSums(clust.total.a[-clust, , drop=F] + clust.total.b[-clust, , drop=F])
    return(a.without.score)
  }))
  a.without.scores[is.na(a.scores)] = NA
  return(list(t(a.scores), t(a.without.scores)))
}

get.num.umis <- function(cells) {
  umis = do.call(rbind, mclapply(1:length(cells), function(i) {
    if (i %% 100 == 1) {
      print(paste('i is', i))
    }
    cell = cells[i]
    vtrack.name = paste0(cell, '_area')
    gvtrack.create(vtrack.name, cell, "area")
    contacts = gextract(vtrack.name, intervals=gintervals.2d.all(), iterator=c(1e10, 1e10))
    num.contacts = sum(contacts[,vtrack.name])
    num.trans.contacts = sum(contacts[contacts$chrom1 != contacts$chrom2, vtrack.name])
    return(c(num.contacts, num.trans.contacts) / 2)
  }, mc.cores=25))
  rownames(umis) = cells
  colnames(umis) = c('num_contacts', 'num_trans_contacts')
  return(umis)
}

total.ab.format <- function(cells, total.ab) {
  atab = (total.ab %>% dcast(cell ~ f1, value.var='a'))
  rownames(atab) = atab[,1]
  btab = (total.ab %>% dcast(cell ~ f1, value.var='b'))
  rownames(btab) = btab[,1]
  
  atab.fil = atab[cells, -1]
  btab.fil = btab[cells, -1]
  atab.fil[is.na(atab.fil)] = 0
  btab.fil[is.na(btab.fil)] = 0
  return(list(atab=atab.fil, btab=btab.fil))
}

get.all.bins <- function(chrom, start, end, bin.size=2e5) {
  chr.start = floor(start / bin.size) * bin.size + 1
  chr.end = floor(end / bin.size) * bin.size + 1
  bin.coords = seq(chr.start, chr.end, bin.size)
  bin.names = paste0(chrom, '_', trimws(format(bin.coords, scientific=F)))
  return(list(bin.names, bin.coords))
}

get.chrom.sizes <- function(chrom.size.path=CHROM.SIZE.PATH) {
  return(read.table(chrom.size.path, row.names=1)[as.character(1:19), 1])
}

get.gene.metadata.jumps <- function(jump=2e5) {
  tss = gintervals.load("intervs.global.tss")
  tss$bin = jump * floor((tss$start - 1) / jump) + 1
  gdata = tss %>% select("geneSymbol", 'chrom', 'start', 'bin')
}

genes.to.intervs <- function(de.genes, gene.md, is.unique.bins=T, delta=2e6) {
    md.filtered = gene.md[gene.md$geneSymbol %in% de.genes, ]
    #md.filtered = gene.md[!is.na(match.gene.names(de.genes, gene.md$geneSymbol)),]
    if (is.unique.bins) {
      is.dup = duplicated(paste0(md.filtered$chrom, '_', md.filtered$bin))
      md.filtered = md.filtered[!is.dup,]
    }
    is.dup = duplicated(md.filtered$geneSymbol)
    md.filtered = md.filtered[!is.dup,]
    intervs = data.frame(geneSymbol=as.character(md.filtered$geneSymbol), start=md.filtered$start - delta,
                         chrom=md.filtered$chrom, end=md.filtered$start + delta, tss=md.filtered$start)
    return(intervs)
}

get.ab.up.bins <- function(chrom.ab.ret, cluster, compared.to, threshold=0.2) {
  ab.diff = get.ab.differential.bins(chrom.ab.ret, cluster, compared.to, threshold)
  return(ab.diff$bins[ab.diff$sign > 0])
}

group.bins <- function(bins, bin.size=2e5) {
  bin.groups = NULL
  for (i in 1:19) {
    chr.name = paste0('chr', i)
    chr.bins = get.bins.sorted(bins, chr.name)
    if (length(chr.bins) == 0) {
      next
    }
    bin.nums = sapply(strsplit(chr.bins, '_'), function(x) as.numeric(x[[2]]))
    bin.indices = (bin.nums - 1) / bin.size
    seq = seqToIntervals(bin.indices)
    for (j in 1:nrow(seq)) {
      seq.len = seq[j,2] - seq[j,1] + 1
      if (is.null(bin.groups)) {
        bin.groups = data.frame(chrom=chr.name, start=bin.nums[which(bin.indices == seq[j, 1])],
                                                end=bin.nums[which(bin.indices == seq[j, 2])], seq_len=seq.len, stringsAsFactors=F)
      } else {
        bin.groups = rbind(bin.groups, c(chr.name, bin.nums[which(bin.indices == seq[j, 1])],
                                                bin.nums[which(bin.indices == seq[j, 2])], seq.len))
      }
    }
  }
  return(bin.groups)
}

get.4c.trace <- function(all.tracks.conts, interv, tss.dist=1e4, filter.chrom=T) {
  chrom = interv$chrom
  coord = interv$tss
  track.conts.fil = lapply(all.tracks.conts, function(track.conts) {
    if (filter.chrom) {
      track.conts = track.conts[track.conts$chrom1 == chrom & track.conts$chrom2 == chrom,]
    }
    track.conts$dist1 = abs(track.conts$start1 - coord)
    track.conts = track.conts[track.conts$dist1 < tss.dist,]
    return(track.conts)
  })
  df.4c = do.call(cbind, lapply(2:length(all.tracks.conts), function(i) {
    if (nrow(track.conts.fil[[1]]) == 0 | nrow(track.conts.fil[[i]]) == 0) {
      return(data.frame(numeric(0), numeric(0)))
    }
    knn.ret = get.knnx(data.frame(track.conts.fil[[i]]$start1, track.conts.fil[[i]]$start2), 
                       data.frame(track.conts.fil[[1]]$start1, track.conts.fil[[1]]$start2), k=1)
    return(data.frame(track.conts.fil[[1]]$score, track.conts.fil[[i]][knn.ret$nn.index, 'score']))
  }))
  df.4c = df.4c[,c(1, seq(2, ncol(df.4c), 2))]
  colnames(df.4c) = paste0('score', 1:length(all.tracks.conts))
  cont.name = paste0(chrom, '_', trimws(format(track.conts.fil[[1]]$start2, scientific=F)))
  df.4c = df.4c[!duplicated(cont.name),]
  cont.name = cont.name[!duplicated(cont.name)]
  rownames(df.4c) = cont.name

  return(df.4c)
}

get.4c.trace.mult.dist <- function(all.tracks.conts, interv, max.dist=1e6, filter.chrom=T, use.max.dist=F) {
  get.4c.trace.fil <- function(all.tracks.conts, interv, tss.dist, dist.range, filter.chrom=T) {
    gene.4c = get.4c.trace(all.tracks.conts, interv, tss.dist, filter.chrom=filter.chrom)
    coords = as.numeric(sapply(strsplit(rownames(gene.4c), '_'), function(x) x[2]))
    gene.4c$dist = coords - interv$tss
    min.dist = dist.range[1]
    max.dist = dist.range[2]
    gene.4c.fil = gene.4c[between(abs(gene.4c$dist), min.dist, max.dist),]
    return(gene.4c.fil)
  }

  if (use.max.dist) {
    gene.4c.fil = rbind(get.4c.trace.fil(all.tracks.conts, interv, 3e3, c(2e4, 1e5), filter.chrom=filter.chrom),
                        get.4c.trace.fil(all.tracks.conts, interv, 1e4, c(1e5, 5e5), filter.chrom=filter.chrom),
                        get.4c.trace.fil(all.tracks.conts, interv, 3e4, c(5e5, max.dist), filter.chrom=filter.chrom))
  } else {
    gene.4c.fil = rbind(get.4c.trace.fil(all.tracks.conts, interv, 3e3, c(2e4, 1e5), filter.chrom=filter.chrom),
                        get.4c.trace.fil(all.tracks.conts, interv, 1e4, c(1e5, 5e5), filter.chrom=filter.chrom))
  }
                      
  return(gene.4c.fil)
}
