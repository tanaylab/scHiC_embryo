
########## Simulation code #############
simulate.data <- function(ncells, nbins, nclusters=5, repl.duration=1, num.bin.clusters=3) {
  cluster.size = ncells / nclusters
  max.copy.num = 2
  total.cov = matrix(NA, nrow=nbins, ncol=ncells)
  rownames(total.cov) = 1:nbins
  colnames(total.cov) = 1:ncells
  cluster.assignment = unlist(lapply(1:nclusters, function(i) rep(i, cluster.size)))

  bin.clusters = matrix(NA, nrow=nbins, ncol=nclusters)
  num.bins.in.bin.clust = (nbins / 2) / num.bin.clusters
  stopifnot(num.bins.in.bin.clust %% 1 == 0)
  for (i in 1:num.bin.clusters) {
    bin.clusters[(1 + num.bins.in.bin.clust * (i - 1)):(num.bins.in.bin.clust * i), ] = i
  }
  
  step.size = (nbins / 2) / (num.bin.clusters - 1) / nclusters
  stopifnot(step.size %% 1 == 0)
  for (i in 1:(num.bin.clusters - 1)) {
    bin.clust.start = (nbins / 2) + ((nbins / 2) / (num.bin.clusters - 1)) * (i - 1)
    for (j in 1:nclusters) {
      bin.clusters[(bin.clust.start + (j - 1) * step.size + 1):(bin.clust.start + j * step.size), j] = i
      bin.clusters[(bin.clust.start + (j - 1) * step.size + 1):(bin.clust.start + j * step.size), -j] = i + 1
    }
  }

  s.scores = rep(max.min.rescale(1:cluster.size, 1.2, 1.8), nclusters)

  reads.per.bin.and.cell = 10
  reads.per.cell = reads.per.bin.and.cell * nbins
  var.mult = matrix(1 / reads.per.cell, nrow=nbins, ncol=ncells)
  for (i in 1:ncells) {
    all.bin.exps = c()
    for (j in 1:nbins) {
      bin.exp = get.bin.expectation(s.scores[i], bin.clusters[j, cluster.assignment[i]], max.copy.num=max.copy.num, 
                                    num.bin.clusters=num.bin.clusters, repl.duration=repl.duration)
      all.bin.exps = c(all.bin.exps, bin.exp)
    }
    
    cell.cov = rmultinom(1, reads.per.cell, all.bin.exps / sum(all.bin.exps))
    total.cov[, i] = cell.cov
  }
  unnorm.total.cov = total.cov

  total.cov = total.cov / reads.per.cell

  model.params = list()
  model.params$mixture.fractions = rep(1 / nclusters, nclusters)
  model.params$max.copy.num = max.copy.num
  model.params$repl.duration = repl.duration
  model.params$var.mult = var.mult
  model.params$num.bin.clusters = num.bin.clusters
  model.params$bin.probs = rep(1 / nbins, nbins)

  return(list(total.cov=total.cov, model.params=model.params, unnorm.total.cov=unnorm.total.cov, s.scores=s.scores, 
    bin.clusters=bin.clusters, cluster.assignment=cluster.assignment))
}

run.simulation <- function(ncells=120, nbins=240, nclusters=5, repl.duration=1, num.bin.clusters=3) {
  sim.data = simulate.data(ncells, nbins, nclusters=nclusters, repl.duration=repl.duration, num.bin.clusters=num.bin.clusters)
  sim.total.cov = sim.data$total.cov
  model.params = sim.data$model.params

  cells.ord = colnames(sim.total.cov)[order(sim.data$s.scores + rnorm(length(sim.data$s.scores), mean=0, sd=0.1))]
  #init.clusts = init.clustering(sim.total.cov, cells.ord, nclusters, feat.min=-Inf, feat.max=Inf, ratio.min=-Inf, rollmean.k=10)
  init.clusts = sample(1:nclusters, ncells, replace=T)
  e.z = get.init.e.z(colnames(sim.total.cov), nclusters, init.clusts)
  mixture.fractions = colMeans(e.z)
  model.params$e.z = e.z
  model.params$mixture.fractions = mixture.fractions
  
  init.s.scores = max.min.rescale(1:ncells, 1.2, 1.8)
  names(init.s.scores) = cells.ord
  init.s.scores = init.s.scores[colnames(sim.total.cov)]
  model.params$s.scores = init.s.scores
  model.params$lambda = 20

  clust.ret = run.clustering.iteration.ds(sim.total.cov, model.params, start.with.e.step=F, ncores=25)
  return(list(clust.ret=clust.ret, s.scores=sim.data$s.scores, bin.clusters=sim.data$bin.clusters, clusters=sim.data$cluster.assignment))
}

########## Stability code #############
test.embryo.stability <- function(embryo.model.params, num.folds, use.init.func=T) {
  set.seed(42)
  stability.data.dir = STABILITY.DIR

  cells = colnames(embryo.model.params$total.cov)
  ncells = length(cells)
  cell.perm = sample.int(ncells, ncells)
  cell.per.fold = split(cell.perm, 1:num.folds)
  s.scores = embryo.model.params$s.scores
  names(s.scores) = colnames(embryo.model.params$total.cov)
  clusters = model.to.clusters(embryo.model.params)
  rownames(embryo.model.params$bin.clusters) = rownames(embryo.model.params$total.cov)
  
  all.model.params = lapply(1:num.folds, function(i) {
    fold.path = file.path(stability.data.dir, paste('fold', i, use.init.func, sep='_'))
    
    if (!file.exists(fold.path)) {
      if (num.folds == 1) {
        selected.cells = cells
      } else {
        selected.cells = cells[-cell.per.fold[[i]]]
      }
      fold.model.params = cluster.embryo.cells(selected.cells=selected.cells, use.init.func=use.init.func)[[1]]
      save(fold.model.params, file=fold.path)
    } else {
      load(fold.path)
    }
    return(fold.model.params)
  })

  all.s.cors = sapply(1:num.folds, function(i) {
    cur.model.params = all.model.params[[i]]
    cur.s.scores = cur.model.params$s.scores
    return(cor(cur.s.scores, s.scores[colnames(cur.model.params$total.cov)]))
  })

  num.diff.bins = sapply(1:num.folds, function(i) {
    cur.model.params = all.model.params[[i]]
    embryo.model.bins = rownames(embryo.model.params$total.cov)
    cur.model.bins = rownames(cur.model.params$total.cov)
    common.bins = intersect(embryo.model.bins, cur.model.bins)
    rownames(cur.model.params$bin.clusters) = rownames(cur.model.params$total.cov)
    return(sum(embryo.model.params$bin.clusters[common.bins, ] != cur.model.params$bin.clusters[common.bins, ]))
  })

  return(list(all.model.params=all.model.params, all.s.cors=all.s.cors, num.diff.bins=num.diff.bins))
}

########## Cross validation code #############

cross.val.embryo <- function(num.bin.clusters, repl.duration, nclusters=1, lambda=0, bins.to.predict=NULL) {
  decay.metrics = sch_decay_metrics
  ret = get.embryo.cells.model.params()
  model.params = ret$model.params
  unnorm.total.cov = ret$unnorm.total.cov
  cells.ord = colnames(unnorm.total.cov)[order(decay.metrics[colnames(unnorm.total.cov), 'near_f'])]

  model.params$num.bin.clusters = num.bin.clusters
  model.params$repl.duration = repl.duration
  model.params$lambda = lambda

  return(do.model.selection.specific(unnorm.total.cov, model.params, num.bin.clusters, repl.duration, 
         nclusters, lambda, init.clust.func=embryo.init.clusters.func, cells.ord=cells.ord, bins.to.predict=bins.to.predict))
}

do.model.selection <- function(unnorm.total.cov, orig.model.params, bins.groups=NULL, max.num.bin.clusters = 6) {
  # bins.groups should be a list whose length is the number of bin groups, and each 
  # value in that list is a vector with T/F for each bin saying whether it belongs to the group
  if (is.null(bins.groups)) {
    bins.groups = lapply(1:19, function(i) grepl(paste0('chr', i, '_'), rownames(unnorm.total.cov)))
  }

  all.mean.errors = matrix(NA, nrow=max.num.bin.clusters, ncol=max.num.bin.clusters)
  for (num.bin.clusters in 3:max.num.bin.clusters) {
    for (repl.duration in 1:(num.bin.clusters - 1)) {
      orig.model.params$num.bin.clusters = num.bin.clusters
      orig.model.params$repl.duration = repl.duration
      mean.error = mean(cluster.cross.val(unnorm.total.cov, orig.model.params, bins.groups))
      all.mean.errors[num.bin.clusters, repl.duration] = mean.error
    }
  }
  return(all.mean.errors)
}

do.model.selection.specific <- function(unnorm.total.cov, orig.model.params, num.bin.clusters, repl.duration, nclusters, lambda, bins.groups=NULL,
                                        init.clust.func=NULL, cells.ord=NULL, bins.to.predict=NULL) {
  set.seed(42)
  if (any(is.null(bins.groups))) {
    bins.groups = lapply(1:19, function(i) grepl(paste0('chr', i, '_'), rownames(unnorm.total.cov)))
  }
  orig.model.params$num.bin.clusters = num.bin.clusters
  orig.model.params$repl.duration = repl.duration
  orig.model.params$lambda = lambda 
  return(cluster.cross.val(unnorm.total.cov, orig.model.params, nclusters, bins.groups, init.clust.func, cells.ord, bins.to.predict=bins.to.predict))
}

cluster.cross.val <- function(unnorm.total.cov, orig.model.params, nclusters,
                    bins.groups, init.clust.func, cells.ord, loc.prob.vector=NULL, bins.to.predict=NULL) {

  set.seed(42)
  cross.val.data.dir = CROSS.VAL.DIR

  if (is.null(loc.prob.vector)) {
    loc.prob.vector = orig.model.params$bin.probs
  }

  #if (is.null(bins.to.predict)) {
  #  bins.to.predict = 1:nrow(unnorm.total.cov)
  #}

  locs.per.fold = list()
  all.selected.locs = rownames(unnorm.total.cov)
  cell.unnorm.total.cov = colSums(unnorm.total.cov)
  all.locs.total.cov = t(t(unnorm.total.cov) / cell.unnorm.total.cov) 
  all.locs.var.mult = t(1 / (cell.unnorm.total.cov %*% t(rep(1, nrow(unnorm.total.cov)))))

  num.cells.folds = 10
  cell.perm = sample.int(ncol(unnorm.total.cov), ncol(unnorm.total.cov))
  cell.per.fold = split(cell.perm, 1:num.cells.folds)

  num.loc.folds = length(bins.groups)
  model.pred = all.locs.total.cov
  model.pred[,] = NA
  expected.pred = all.locs.total.cov
  expected.pred[,] = NA
  prediction.file.path = file.path(cross.val.data.dir, paste('pred', orig.model.params$num.bin.clusters, orig.model.params$repl.duration, 
    nclusters, orig.model.params$lambda, sep='_'))
  if (file.exists(prediction.file.path)) {
    load(prediction.file.path)
    # loads model.pred, expected.pred
    return(list(model.pred=model.pred, expected.pred=expected.pred,
    cell.error=calc.cell.pred.error(model.pred, expected.pred, bins.to.predict), bin.error=calc.bin.pred.error(model.pred, expected.pred)))
  }

  for (i in 1:num.loc.folds) {
    library(parallel)
    model.params = orig.model.params
    is.in.chrom = bins.groups[[i]]
    chrom.indices = which(is.in.chrom)
    non.chrom.indices = 1:nrow(all.locs.total.cov)
    non.chrom.indices = non.chrom.indices[-chrom.indices]
    current.fold.locs = all.selected.locs[chrom.indices]
    locs.per.fold[[i]] = current.fold.locs
    total.cov = unnorm.total.cov[non.chrom.indices,]
    cell.total.cov = colSums(total.cov)
    total.cov = t(t(total.cov) / cell.total.cov) 
    unnorm.cur.loc.probs = loc.prob.vector[-chrom.indices]
    cur.loc.prob.vector = unnorm.cur.loc.probs / sum(unnorm.cur.loc.probs)
    var.mult = t(1 / (cell.total.cov %*% t(rep(1, nrow(total.cov)))))
    model.params$var.mult = var.mult
    model.params$bin.probs = cur.loc.prob.vector
    init.clusters = init.clust.func(total.cov, cells.ord, nclusters)
    model.params$e.z = get.init.e.z(colnames(total.cov), nclusters, init.clusters)
    model.params$mixture.fractions = colSums(model.params$e.z) / sum(model.params$e.z)

    iter.file.path = file.path(cross.val.data.dir, paste('iter', model.params$num.bin.clusters, model.params$repl.duration, 
                                                                 nclusters, model.params$lambda, i, sep='_'))
    if (file.exists(iter.file.path)) {
      load(iter.file.path, verbose=T)
      # loads iter.ret
    } else {
      iter.ret = run.clustering.iteration.ds(total.cov, model.params, FALSE)
      save(iter.ret, file=iter.file.path)
    }
    model.params = iter.ret[[1]]
    cur.e.z = iter.ret[[2]]

    ncells = ncol(total.cov)
    num.bin.clusters = model.params$num.bin.clusters
    nclusters = ncol(cur.e.z)

    cov.sum = matrix(NA, nrow=ncells, ncol=nclusters)
    exp.cache = matrix(NA, nrow=ncells, ncol=num.bin.clusters)
    for (i in 1:ncells) {
      for (l in 1:num.bin.clusters) {
        exp.cache[i, l] = get.bin.expectation(model.params$s.scores[i], l, model.params$max.copy.num, 
                               model.params$num.bin.clusters, model.params$repl.duration)
      }
    }

    # I estimate the cov sum based on the loc probs of the whole data, but for only the bins not in the hidden chromosome. 
    # I then add another factor which estimates the contribution of the missing chromsome to the cov.sum.
    for (k in 1:nclusters) {
      cov.sum[, k] = Reduce("+", lapply(1:length(unnorm.cur.loc.probs), function(j) unnorm.cur.loc.probs[j] * exp.cache[, model.params$bin.clusters[j, k]]))
    }
    cov.sum = cov.sum + sum(loc.prob.vector[chrom.indices]) * cov.sum

    all.folds.error = c()
    all.folds.cor = c()
    for (cell.fold in 1:num.cells.folds) {
      cells.to.leave.out = cell.per.fold[[cell.fold]]
      # The +1 thing is a hack so that the input won't contain NA
      bin.clusters.init = matrix(model.params$num.bin.clusters + 1, nrow=nrow(unnorm.total.cov), ncol=nclusters)
      bin.clusters.init[non.chrom.indices,] = model.params$bin.clusters
      bin.clusters = all.bins.m.step.ds(all.locs.total.cov[,-cells.to.leave.out], 
                                        s.scores=model.params$s.scores[-cells.to.leave.out], e.z=model.params$e.z[-cells.to.leave.out, , drop=F], 
                                        num.bin.clusters=model.params$num.bin.clusters, 
                                        max.copy.num=model.params$max.copy.num, 
                                        repl.duration=model.params$repl.duration, 
                                        bin.probs=loc.prob.vector, 
                                        var.mult=all.locs.var.mult[,-cells.to.leave.out], 
                                        cov.sum.input=cov.sum[-cells.to.leave.out, , drop=F], bin.clusters=bin.clusters.init,
                                        bins.to.change=chrom.indices)
  
      pred.total.cov = all.locs.total.cov
      pred.total.cov[,] = NA
      for (cell.to.predict in cells.to.leave.out) {
        bin.exps = loc.prob.vector * get.bin.expectation(model.params$s.scores[cell.to.predict], bin.clusters, 
                                       model.params$max.copy.num, model.params$num.bin.clusters, model.params$repl.duration)
        bin.exps = t(t(bin.exps) / colSums(bin.exps))
        pred.total.cov[, cell.to.predict] = rowSums(t(cur.e.z[cell.to.predict,] * t(bin.exps)))
      }

      if (is.null(bins.to.predict)) {
        chrom.bins.to.predict = chrom.indices
      } else {
        chrom.bins.to.predict = which(is.in.chrom & (rownames(unnorm.total.cov) %in% bins.to.predict))
      }

      expected.total.cov = all.locs.total.cov[chrom.bins.to.predict, cells.to.leave.out]
      pred.total.cov.filtered = pred.total.cov[chrom.bins.to.predict, cells.to.leave.out]

      model.pred[chrom.bins.to.predict, cells.to.leave.out] = pred.total.cov.filtered
      expected.pred[chrom.bins.to.predict, cells.to.leave.out] = expected.total.cov 
    }    
  }
  
  save(model.pred, expected.pred, file=prediction.file.path)
  
  return(list(model.pred=model.pred, expected.pred=expected.pred,
    cell.error=calc.cell.pred.error(model.pred, expected.pred, bins.to.predict), bin.error=calc.bin.pred.error(model.pred, expected.pred)))
}

get.replication.cross.val.errors <- function(max.num.bin.clusters = 6, bins.to.predict=NULL) {
  res.mat = matrix(NA, nrow=max.num.bin.clusters, ncol=max.num.bin.clusters)
  for (num.bin.clusters in 3:max.num.bin.clusters) {
    for (repl.duration in 1:(num.bin.clusters - 1)) {
      res.mat[num.bin.clusters, repl.duration] = cross.val.embryo(num.bin.clusters, repl.duration, bins.to.predict=bins.to.predict)$cell.error
    }
  }
  return(res.mat)
}

get.replication.cross.val.errors.specific <- function(num.bin.clusters, repl.duration, bins.to.predict=NULL) {
  return(cross.val.embryo(num.bin.clusters, repl.duration, bins.to.predict=bins.to.predict)$cell.error)
}

get.clustering.cross.val.errors <- function(num.bin.clusters, repl.duration, nclust.options, lambda.options) {
  res.mat = matrix(NA, nrow=length(nclust.options), ncol=length(lambda.options))
  for (nclust.index in 1:length(nclust.options)) {
    nclust = nclust.options[nclust.index]
    for (lambda.index in 1:length(lambda.options)) {
      lambda = lambda.options[lambda.index]
      res.mat[nclust.index, lambda.index] = cross.val.embryo(num.bin.clusters, repl.duration, nclust, lambda)$cell.error
    }
  }
  return(res.mat)
}

calc.cell.pred.error <- function(model.pred, expected.pred, bins.to.predict=NULL) {
  mean(sapply(1:ncol(model.pred), function(i) {
    cell.pred = model.pred[, i]
    cell.exp = expected.pred[, i]
    predicted.bins = which(!(is.na(cell.pred) | is.na(cell.exp)))

    if (!is.null(bins.to.predict)) {
      predicted.bins = intersect(predicted.bins, bins.to.predict)
    }
    
    return(cor(cell.pred[predicted.bins], cell.exp[predicted.bins]))
  }))
}

calc.bin.pred.error <- function(model.pred, expected.pred) {
  not.na = !(apply(model.pred, 1, function(x) any(is.na(x))) | apply(expected.pred, 1, function(x) any(is.na(x))))
  mean(sapply(which(not.na), function(i) {
    bin.pred = model.pred[i, ]
    bin.exp = expected.pred[i, ]
    return(cor(bin.pred, bin.exp))
  }))
}


########## Validation plots code #############
plot.stability <- function(model.params, stability.ret, fig.dir=get.fig.dir(), fig.name='sampling_s_stability') {
	if (is.null(model.params$cols)) {
          model.params$cols = MODEL.CLUSTER.COLORS[1:length(model.params$mixture.fractions)]
	}
	print(stability.ret$all.s.cors)
	print(mean(stability.ret$all.s.cors))

	# this number has no meaning, because I currently do no cluster matching
	print(stability.ret$num.diff.bins)
	print(mean(stability.ret$num.diff.bins))

	s.scores = model.params$s.scores
  	names(s.scores) = colnames(model.params$total.cov)

	png(file.path(fig.dir, paste0(fig.name, '.png')))
	plot(s.scores[colnames(stability.ret$all.model.params[[1]]$total.cov)], 
		stability.ret$all.model.params[[1]]$s.scores, pch=19, cex=1.5)
	dev.off()

	new.cols = c()
	orig.cls = model.to.clusters(model.params)
	new.cls = model.to.clusters(stability.ret$all.model.params[[1]])
	stopifnot(length(table(orig.cls)) == length(table(new.cls)))
	cls.table = table(orig.cls[names(new.cls)], new.cls)

	png(file.path(fig.dir, paste0(fig.name, '_orig_col.png')))
	plot(s.scores[colnames(stability.ret$all.model.params[[1]]$total.cov)], 
		stability.ret$all.model.params[[1]]$s.scores, pch=19, cex=1.5, col=model.params$cols[orig.cls[names(new.cls)]])
	dev.off()

	for (i in 1:length(table(new.cls))) {
	  new.cols = c(new.cols, model.params$cols[which.max(cls.table[,i])])
	}
	png(file.path(fig.dir, paste0(fig.name, '_new_col.png')))
	plot(s.scores[colnames(stability.ret$all.model.params[[1]]$total.cov)], 
		stability.ret$all.model.params[[1]]$s.scores, pch=19, cex=1.5, col=new.cols[new.cls])
	dev.off()
}

plot.lambda.cross.val.errors <- function(lambda.cross.val.errors, fig.dir=get.fig.dir()) {
  cols = brewer.pal(8, 'Dark2')
  png(file.path(fig.dir, 'lambda_cross_val.png'), w=1000, h=1000)
  ylim = range(lambda.cross.val.errors, na.rm=T)
  for (i in 1:nrow(lambda.cross.val.errors)) {
    if (i == 1) {
      plot(as.numeric(colnames(lambda.cross.val.errors)), lambda.cross.val.errors[i,], col=cols[i], type='l', ylim=ylim, lwd=15)
    } else {
      lines(as.numeric(colnames(lambda.cross.val.errors)), lambda.cross.val.errors[i,], col=cols[i], lwd=15)
    }
  }
  legend('bottomright', legend=rownames(lambda.cross.val.errors), col=cols[1:(nrow(lambda.cross.val.errors))], lty=1, title='# clusters', cex=3, lwd=15)
  grid(col='black')
  dev.off()
}

plot.repli.cross.val.errors <- function(repli.cross.val.errors, repli.cross.val.errors.specific, fig.dir=get.fig.dir()) {
        cols = brewer.pal(8, 'Dark2')
	png(file.path(fig.dir, 'repli_cross_val.png'), w=1000, h=1000)
	ylim = range(repli.cross.val.errors, na.rm=T)
	for (i in 3:nrow(repli.cross.val.errors)) {
		if (i == 3) {
			plot(repli.cross.val.errors[i,], col=cols[i - 2], type='l', ylim=ylim, lwd=15)
		} else {
			lines(repli.cross.val.errors[i,], col=cols[i - 2], lwd=15)
		}
	}
	legend('bottomright', legend=3:nrow(repli.cross.val.errors), col=cols[1:(nrow(repli.cross.val.errors)-2)], lty=1, title='# replication regimes', cex=3, lwd=15)
	grid(col='black')
	dev.off()

	png(file.path(fig.dir, 'repli_cross_val_specific.png'), w=1000, h=1000)
	plot(as.numeric(names(repli.cross.val.errors.specific)), repli.cross.val.errors.specific, type='l', lwd=16)	
	points(12, repli.cross.val.errors.specific['12'], pch=19, col='red', cex=10)
	grid(col='black')
	dev.off()
}

plot.cv.num.cls <- function(num.clust.cross.val.errors, fig.dir) {
  png(file.path(fig.dir, 'cross_val_num_cls.png'), height=400, width=400)
  plot(1:nrow(num.clust.cross.val.errors), num.clust.cross.val.errors[,1], pch=19, cex=2)
  lines(1:nrow(num.clust.cross.val.errors), num.clust.cross.val.errors[,1], lwd=4)
  grid(col='black')
  dev.off()
}

plot.cross.val.cor <- function(cross.val.ret, unnorm.total.cov, fig.dir=get.fig.dir()) {
	cell.cors = sapply(1:ncol(unnorm.total.cov), function(i) {
		cor(cross.val.ret$expected.pred[,i], cross.val.ret$model.pred[,i])
	})
	png(file.path(fig.dir, 'pred_cell_vs_coverage.png'))
	plot(colSums(unnorm.total.cov), cell.cors, pch=19, cex=2)
	dev.off()
}
