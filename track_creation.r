
create.pool.track <- function(conts, track.name) {
  conts = conts[!duplicated(paste(conts$chrom1, conts$start1, conts$start2, sep='_')),]
  intervs = gintervals.2d(chroms1=conts$chrom1, 
                          starts1=conts$start1, 
                          ends1=conts$end1, 
                          chroms2=conts$chrom2,
                          starts2=conts$start2,
                          ends2=conts$end2)
  gtrack.2d.create(track.name, 'track description', intervals=intervs, values=rep(1, nrow(conts)))
}

# To create ecto_for_meso_comparison and meso_for_ecto_comparison
create.ecto.meso.pool.tracks <- function(all.cluster.assignment) {
  set.misha(SCHIC.MISHA.PATH)
  all.conts.fil = list()
  for (i in 1:2) {
    print(i)
    cells = names(all.cluster.assignment)[all.cluster.assignment == i]
    all.conts = do.call(rbind, mclapply(1:length(cells), function(j) {
      if (j %% 100 == 1) {
        print(j)
      }
      cell = cells[j]
      gextract(cell, gintervals.2d.all(), colnames='score')
    }, mc.cores=25))
    all.conts.fil[[i]] = filter(all.conts, chrom1 == chrom2 & !(chrom1 %in% c('chrX', 'chrY', 'chrM')))
  }

  smaller.cluster = which.min(c(nrow(all.conts.fil[[1]]), nrow(all.conts.fil[[2]])))
  larger.cluster.conts = all.conts.fil[[3 - smaller.cluster]]
  larger.cluster.conts.one.sided = filter(larger.cluster.conts, start1 < start2)
  larger.cluster.ds = larger.cluster.conts.one.sided[sample(1:nrow(larger.cluster.conts.one.sided), nrow(all.conts.fil[[smaller.cluster]]) / 2, replace=F),]
  larger.cluster.other.side = larger.cluster.ds[,c(4:6, 1:3, 7:ncol(larger.cluster.ds))]
  colnames(larger.cluster.other.side) = colnames(all.conts.fil[[smaller.cluster]])
  all.larger.cluster.conts = rbind(larger.cluster.ds, larger.cluster.other.side)

  set.misha(SHAMAN.MISHA.PATH)
  create.pool.track(all.conts.fil[[smaller.cluster]], 'ecto_for_meso_comparison')
  create.pool.track(all.larger.cluster.conts, 'meso_for_ecto_comparison')

}


create.cell.equal.sized.tracks <- function(track.names, new.tracks.suffix='_equal_size') {
  num.contacts = c()
  for (track.name in track.names) {
    print(track.name)
    contacts = gextract(track.name, intervals=gintervals.2d.all())[,1:6]
    good.conts = (contacts$chrom1 == contacts$chrom2) & (!(contacts$chrom1 %in% c('chrX', 'chrY', 'chrM')))
    num.contacts = c(num.contacts, sum(good.conts))
  }
  num.contacts.to.ds = min(num.contacts)

  for (track.name in track.names) {
    print(track.name)
    contacts = gextract(track.name, intervals=gintervals.2d.all())[,1:6]
    good.conts = (contacts$chrom1 == contacts$chrom2) & (!(contacts$chrom1 %in% c('chrX', 'chrY', 'chrM')))
    contacts = contacts[good.conts,]
    sampled.contact.indices = sample(1:nrow(contacts), num.contacts.to.ds)
    sampled.contacts = contacts[sampled.contact.indices,]
    gtrack.2d.create(paste0(track.name, new.tracks.suffix), 'track description', intervals=sampled.contacts, values=rep(1, nrow(sampled.contacts)))
  }
}

create.cell.pool.track <- function(cells, track.name, 
                                   ind.track.misha=SCHIC.MISHA.PATH, pool.track.misha=SHAMAN.MISHA.PATH, tmp.pool.dir=TMP.POOL.DIR) {
  num.cores = 25
  gsetroot(ind.track.misha)
  core.assignment = rep(1:num.cores, length.out=length(cells))
  tmp.file.dir = file.path(tmp.pool.dir, track.name)
  dir.create(tmp.file.dir, showWarnings=F)
  mclapply(1:num.cores, function(cur.core) {
    contact.file.path = file.path(tmp.file.dir, cur.core)
    if (file.exists(contact.file.path)) {
      next
    }
    cur.cells = cells[core.assignment == cur.core]
    all.cell.conts = NULL
    for (cell in cur.cells) {
      track.conts = gextract(cell, gintervals.2d.all())[,1:6]
      all.cell.conts = rbind(all.cell.conts, track.conts)
    }
    all.cell.conts$value = 1
    all.cell.conts = format(all.cell.conts, scientific=F)
    write.table(all.cell.conts, file=contact.file.path, sep='\t', quote=F, row.names=F)
  }, mc.cores=num.cores)
  gsetroot(pool.track.misha)
  gtrack.2d.import_contacts(track.name, track.name, file.path(tmp.file.dir, 1:num.cores))
}


create.all.cells.pooled.tracks <- function(decay.metrics, esc.decay.metrics) {
  erys = get.erys()
  emb.non.erys = setdiff(rownames(decay.metrics), erys)
  escs = esc.decay.metrics$cell
  create.cell.pool.track(erys, ERY.POOL.TRACK.NAME)
  create.cell.pool.track(emb.non.erys, EMB.POOL.TRACK.NAME)
  create.cell.pool.track(escs, ESC.POOL.TRACK.NAME)
}


# To create the C1-C3 pool tracks and their shaman tracks
create.emb.esc.ery.tracks <- function() {
  misha.db.dir = SHAMAN.MISHA.PATH
  setwd(TMP.CREATE.TRACK.DIR)
  gsetroot(misha.db.dir)
  options(shaman.mc_support=1)

  orig.track.names = c(EMB.POOL.TRACK.NAME, ESC.POOL.TRACK.NAME, ERY.POOL.TRACK.NAME)
  ds.track.names = paste0(c(EMB.POOL.TRACK.NAME, ESC.POOL.TRACK.NAME, ERY.POOL.TRACK.NAME), '_equal_size')
  create.all.cells.pooled.tracks(sch_decay_metrics, esc_sch_decay_metrics)
  create.cell.equal.sized.tracks(orig.track.names)

  for (i in 1:3) {
    shaman_shuffle_hic_track(track_db=misha.db.dir, obs_track_nm=ds.track.names[i], work_dir=TMP.CREATE.TRACK.DIR)
  }

  for (i in 1:3) {
    obs.track.name = ds.track.names[i]
    #shaman_score_hic_track(track_db=misha.db.dir, work_dir=TMP.CREATE.TRACK.DIR, score_track_nm=paste0(obs.track.name, '_scores'), obs_track_nms=obs.track.name, near_cis=1e7)
    shaman_score_hic_track(track_db=misha.db.dir, work_dir=TMP.CREATE.TRACK.DIR, score_track_nm=paste0(obs.track.name, '_matshuff_scores'), obs_track_nms=obs.track.name, near_cis=1e7)
  }
}

# To create the ecto and meso pool tracks and their shaman tracks
create.emb.proper.tracks <- function() {
  load(CLUSTERING.PATH, v=T)

  set.misha(SCHIC.MISHA.PATH)
  create.cell.pool.track(names(all.cluster.assignment)[all.cluster.assignment == 1], EMB.CLS.ANALYSIS.TRACK.NAMES[1])
  create.cell.pool.track(names(all.cluster.assignment)[all.cluster.assignment == 2], EMB.CLS.ANALYSIS.TRACK.NAMES[2])

  misha.db.dir = SHAMAN.MISHA.PATH
  setwd(TMP.CREATE.TRACK.DIR)
  gsetroot(misha.db.dir)
  options(shaman.mc_support=1)

  create.cell.equal.sized.tracks(EMB.CLS.ANALYSIS.TRACK.NAMES[1:2])
  ds.tnames = paste0(EMB.CLS.ANALYSIS.TRACK.NAMES, '_equal_size')[1:2]
  for (i in 1:2) {
    shaman_shuffle_hic_track(track_db=misha.db.dir, obs_track_nm=ds.tnames[i], work_dir=TMP.CREATE.TRACK.DIR)
  }

  for (i in 1:2) {
    obs.track.name = ds.tnames[i]
    #shaman_score_hic_track(track_db=misha.db.dir, work_dir=TMP.CREATE.TRACK.DIR, score_track_nm=paste0(obs.track.name, '_scores'), obs_track_nms=obs.track.name, near_cis=1e7)
    shaman_score_hic_track(track_db=misha.db.dir, work_dir=TMP.CREATE.TRACK.DIR, score_track_nm=paste0(obs.track.name, '_matsuff_scores'), obs_track_nms=obs.track.name, near_cis=1e7)
  }

}

# To create the ery and e10 pool tracks for their comparison
create.emb.proper.tracks <- function() {
  e10.cells = rownames(e10_decay_metrics)
  erys = get.erys()

  tnames = c('all_erys_pooled_paper_ery_comp', 'all_charlie_pooled_paper_ery_comp')
  set.misha(SCHIC.MISHA.PATH)
  create.cell.pool.track(erys, tnames[1])
  create.cell.pool.track(e10.cells, tnames[2])

  misha.db.dir = SHAMAN.MISHA.PATH
  setwd(TMP.CREATE.TRACK.DIR)
  gsetroot(misha.db.dir)
  options(shaman.mc_support=1)

  create.cell.equal.sized.tracks(tnames)
  ds.tnames = paste0(tnames, '_equal_size')[1:2]
  for (i in 1:2) {
    shaman_shuffle_hic_track(track_db=misha.db.dir, obs_track_nm=ds.tnames[i], work_dir=TMP.CREATE.TRACK.DIR)
  }

  for (i in 1:2) {
    obs.track.name = ds.tnames[i]
    #shaman_score_hic_track(track_db=misha.db.dir, work_dir=TMP.CREATE.TRACK.DIR, score_track_nm=paste0(obs.track.name, '_scores'), obs_track_nms=obs.track.name, near_cis=1e7)
    shaman_score_hic_track(track_db=misha.db.dir, work_dir=TMP.CREATE.TRACK.DIR, score_track_nm=paste0(obs.track.name, '_matsuff_scores'), obs_track_nms=obs.track.name, near_cis=1e7)
  }

}

