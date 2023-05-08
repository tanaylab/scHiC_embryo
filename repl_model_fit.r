
######### Model optimization ############

run.clustering.iteration.ds <- function(filtered.total.cov, model.params, start.with.e.step=T, ncores=25) {
  # filtered.total.cov is assumed to be normalized (sum to 1), and already contain only sufficiently covered bins.
  # Rows are bins and cols are cells.
  start.time = Sys.time()
  orig.model.params = model.params
  mixture.fractions = model.params$mixture.fractions
  max.copy.num = model.params$max.copy.num
  repl.duration = model.params$repl.duration
  num.bin.clusters = model.params$num.bin.clusters
  lambda = model.params$lambda
  
  # assume that the input total.cov is already normalized
  total.cov = filtered.total.cov
  ncells = ncol(total.cov)
  nbins = nrow(total.cov)
  nclusters = length(mixture.fractions)
  s.scores = model.params$s.scores
  bin.clusters = model.params$bin.clusters
  e.z = model.params$e.z
  var.mult = model.params$var.mult
  bin.probs = model.params$bin.probs
  stopifnot(length(dim(var.mult)) == 2)
  cell.clusters = apply(e.z, 1, which.max)

  # Stopping parameters
  max.num.m.step.iters = model.params$max.num.m.step.iters
  diff.bin.thresh = model.params$diff.bin.thresh
  diff.cell.thresh = model.params$diff.cell.thresh

  if (is.null(diff.bin.thresh)) {
    diff.bin.thresh = nbins * nclusters / 1e3  
  }
  if (is.null(diff.cell.thresh)) {
    diff.cell.thresh = ncells / 1e2
  }
  

  if (start.with.e.step) {
    e.z = do.e.step.ds(total.cov, model.params, bin.clusters, mixture.fractions, s.scores, bin.probs=bin.probs, var.mult=var.mult, ncores=ncores)
  }

  iter = 0
  continue.optimization = T
  while (continue.optimization) {
    message('Starting optimization iteration')
    iter = iter + 1
    if (!is.null(max.num.m.step.iters) && iter == max.num.m.step.iters) {
      continue.optimization = F
    }

    old.bin.clusters = bin.clusters
    bin.clusters = all.bins.m.step.ds(total.cov, s.scores=s.scores, e.z=e.z, 
                       num.bin.clusters=num.bin.clusters, max.copy.num=max.copy.num, repl.duration=repl.duration, 
                       bin.probs=bin.probs, var.mult=var.mult, bin.clusters=bin.clusters, lambda=lambda, ncores=ncores)
    diff.in.bin.clusters = sum(bin.clusters != old.bin.clusters)
    message('difference in bin.clusters: ', diff.in.bin.clusters)
    is.bin.unchanged = (iter > 1) & (diff.in.bin.clusters <= diff.bin.thresh)
    message('Total time after bins optim: ', Sys.time() - start.time)

    e.z = do.e.step.ds(total.cov, model.params, bin.clusters, mixture.fractions, s.scores, bin.probs=bin.probs, var.mult=var.mult, ncores=ncores)

    mixture.fractions = colSums(e.z) / sum(e.z)
    s.scores = unlist(mclapply(1:ncells, function(index) single.cell.m.step.ds(index, total.cov=total.cov, 
                        bin.clusters=bin.clusters, e.z=e.z, num.bin.clusters=num.bin.clusters, 
                        max.copy.num=max.copy.num, repl.duration=repl.duration, bin.probs=bin.probs, var.mult=var.mult), mc.cores=ncores))
    
    message('Total time after s.score optim: ', Sys.time() - start.time)
    e.z = do.e.step.ds(total.cov, model.params, bin.clusters, mixture.fractions, s.scores, bin.probs=bin.probs, var.mult=var.mult, ncores=ncores)
    
    old.cell.clusters = cell.clusters
    cell.clusters = apply(e.z, 1, which.max)
    diff.in.cell.clusters = sum(cell.clusters != old.cell.clusters)
    is.clusters.unchanged = (iter > 1) & (diff.in.cell.clusters <= diff.cell.thresh)

    if (is.clusters.unchanged & is.bin.unchanged) {
      continue.optimization = F
    }
  }

  model.params$mixture.fractions = mixture.fractions
  model.params$s.scores = s.scores
  model.params$bin.clusters = bin.clusters
  model.params$e.z = e.z
  model.params$total.cov = total.cov
  return(list(model.params=model.params, e.z=e.z))
}


do.e.step.ds <- function(total.cov, model.params, bin.clusters, mixture.fractions, s.scores, bin.probs, var.mult, ncores=25) {
  max.copy.num = model.params$max.copy.num
  num.bin.clusters = model.params$num.bin.clusters
  repl.duration = model.params$repl.duration
  ncells = ncol(total.cov)
  nclusters = length(mixture.fractions)
  if (nclusters == 1) {
    return(matrix(1, ncol=1, nrow=ncells))
  }
  
  # e-step: calculate expectations
  e.z = do.call(rbind, mclapply(1:ncells, function(i) {
    bin.exps = bin.probs * get.bin.expectation(s.scores[i], bin.clusters, max.copy.num, num.bin.clusters, repl.duration)
    bin.exps = t(t(bin.exps) / colSums(bin.exps))
    log_pr_count_ij_cond_sample_i_cluster_k = colSums(dnorm(total.cov[, i], mean=bin.exps, sd=sqrt(bin.exps * var.mult[, i]), log=T))
    cell.e.z = log(mixture.fractions) + log_pr_count_ij_cond_sample_i_cluster_k
    return(cell.e.z)
  }, mc.cores=ncores))

  e.z.norm = e.z
  for (k in 1:nclusters) {
    e.z.norm[, k] = 1 / (Reduce("+", lapply(1:nclusters, function(k2) exp(e.z[, k2] - e.z[, k]))))
  }

  e.z = e.z.norm
  return(e.z)
}


single.cell.m.step.ds <- function(cell.index, total.cov, bin.clusters, e.z, num.bin.clusters, max.copy.num, repl.duration, bin.probs, var.mult, grid.step=0.01) {
  # Optimize for each sub part separately and choose the best statistics
  i = cell.index
  nclusters = ncol(e.z)
  nbins = nrow(total.cov)

  obj.func <- function(params) {
    s.score = params[1]
    all.bin.exps = matrix(NA, nrow=nbins, ncol=nclusters)
    bin.exp.cache = sapply(1:num.bin.clusters, function(l) get.bin.expectation(s.score, l, max.copy.num, 
                                          num.bin.clusters, repl.duration))

    for (k in 1:nclusters) {
      all.bin.exps[,k] = bin.probs * bin.exp.cache[bin.clusters[,k]]
    }

    phase.obj = sum(-e.z[i,] * sapply(1:nclusters, function(k) {
       bin.exp = all.bin.exps[,k] / sum(all.bin.exps[,k]) 
       return(sum(dnorm(total.cov[,i], mean=bin.exp, sd=sqrt(var.mult[,i] * bin.exp), log=T)))
        }))

    return(phase.obj)
  }

  obj.values = sapply(seq(1, 2, grid.step), obj.func)
  return(seq(1, 2, grid.step)[which.min(obj.values)])
}

all.bins.m.step.ds <- function(total.cov, s.scores, e.z, num.bin.clusters, max.copy.num, repl.duration, bin.probs, var.mult, 
                               bin.clusters=NULL, cov.sum.input=NULL, bins.to.change=NULL, lambda=0, ncores=25) {
  nclusters = ncol(e.z)
  ncells = ncol(total.cov)
  nbins = nrow(total.cov)

  exp.cache = matrix(NA, nrow=ncells, ncol=num.bin.clusters)
  for (l in 1:num.bin.clusters) {
    exp.cache[, l] = get.bin.expectation(s.scores, l, max.copy.num, num.bin.clusters, repl.duration)
  }
 
  single.bin.m.step.ds <- function(j, cov.sum, bin.clusters) {

    is.bin.cluster = all(!is.na(bin.clusters))
    if (is.bin.cluster) {
      clust.ratios = colSums(e.z) / sum(e.z)
      # Probably makes more sense to have abs here instead of square, but there results were done like this...
      prob.per.bin.clust = sapply(1:num.bin.clusters, function(l) sum(clust.ratios * (bin.clusters[j, ] - l)**2)) 
      global.bin.clust = which.min(prob.per.bin.clust)
    } else {
      lambda = 0
      global.bin.clust = 1 # unused
    }

    all.likelihoods = do.call(rbind, lapply(1:nclusters, function(k) {
      frac = bin.probs[j] * exp.cache / cov.sum[, k] 
      lik.per.bin.clust = colSums(-e.z[, k] * dnorm(total.cov[j,], mean=frac, sd=sqrt(frac * var.mult[j,]), log=T))
      #lik.with.regularization = lik.per.bin.clust + lambda * abs(1:num.bin.clusters - global.bin.clust)
      lik.with.regularization = lik.per.bin.clust + lambda * colMeans(e.z)[k] * abs(1:num.bin.clusters - global.bin.clust)
      return(lik.with.regularization)
    }))
    return(apply(all.likelihoods, 1, which.min))
  }

  old.bin.clusters = bin.clusters
  
  should.continue = T
  is.parallel = T
  parallel.thresh = (nbins * nclusters) / 1e3
  non.parallel.thresh = (nbins * nclusters) / 1e3
  while (should.continue) {

    is.bin.cluster = !all(is.null(bin.clusters))
    if (!is.bin.cluster) {
      bin.clusters = matrix(NA, nrow=nbins, ncol=nclusters)
    }

    # init cov.sum
    if (all(is.null(cov.sum.input))) {
      cov.sum = matrix(NA, nrow=ncells, ncol=nclusters)
      for (k in 1:nclusters) {
        if (is.bin.cluster) {
          cov.sum[, k] = Reduce("+", lapply(1:nbins, function(j) bin.probs[j] * exp.cache[, bin.clusters[j, k]]))
        } else {
          # Note that this assumes that the coverage increases linearly.
          cov.sum[, k] = 1 + (s.scores - 1) * (max.copy.num - 1)
        }
      }
    } else {
      cov.sum = cov.sum.input
    }

    if (all(is.null(bins.to.change))) {
      bins.to.change = 1:nbins
    }

    if (is.parallel) {
      bin.clusters[bins.to.change,] = do.call(rbind, mclapply(bins.to.change, function(j) single.bin.m.step.ds(j, cov.sum, bin.clusters), mc.cores=ncores))
    } else {
      
      for (j in bins.to.change) {
        new.bin.clusters = single.bin.m.step.ds(j, cov.sum, bin.clusters)
        if (is.bin.cluster) {
          cov.sum = cov.sum - bin.probs[j] * (exp.cache[, bin.clusters[j,]] - exp.cache[, new.bin.clusters])
        }
        bin.clusters[j,] = new.bin.clusters
      }
    }

    print(paste('diff is', sum(bin.clusters != old.bin.clusters)))
    cov.sum.input = NULL # hack to only use the sum the first time

    if (!any(is.null(old.bin.clusters)) & (is.parallel & (sum(bin.clusters != old.bin.clusters) <= parallel.thresh) |
       (!is.parallel & (sum(bin.clusters != old.bin.clusters) <= non.parallel.thresh)))) {
      if (is.parallel) {
        is.parallel = F
      } else {
        should.continue = F
      }
    } 
    old.bin.clusters = bin.clusters
  }
  
  return(bin.clusters)
  
}

model.to.clusters <- function(model.params) {
  clusters = apply(model.params$e.z, 1, which.max)
  names(clusters) = colnames(model.params$total.cov)
  return(clusters)
}

get.bin.expectation <- function(s.score, bin.cluster, max.copy.num=2, num.bin.clusters=2, repl.duration=1) {
  bin.start.repl = 1 + (bin.cluster - 1) / (num.bin.clusters + repl.duration - 1)
  bin.end.repl = 1 + (bin.cluster + repl.duration - 1) / (num.bin.clusters + repl.duration - 1)
  return(ifelse(s.score < bin.start.repl, 1, ifelse(s.score > bin.end.repl, max.copy.num,
      1 + (s.score - bin.start.repl) * ((max.copy.num - 1) / (bin.end.repl - bin.start.repl)))))
}
