load.decay.metrics <- function() {
  load(DECAY.METRICS.PATH)
  sch_decay_metrics <<- sch_decay_metrics
  load(ESC.DECAY.METRICS.PATH)
  esc_sch_decay_metrics <<- esc_sch_decay_metrics
  load(E10.DECAY.METRICS.PATH)
  e10_decay_metrics <<- charlie.decay.metrics
}

load.ab.and.cov <- function(load.heavy=F) {

  load(E10.BIN.CELL.AB.PATH)
  e10_bin_cell_ab <<- bin_cell_ab
  load(E10.BIN.CELL.COV.PATH)
  e10_bin_cell_cov <<- bin_cell_cov
  load(BIN.CELL.COV.PATH)
  bin_cell_cov <<- bin_cell_cov
  bin_cell_cov_ds <<- bin_cell_cov_ds
  load(BIN.CELL.AB.PATH)
  bin_cell_ab <<- bin_cell_ab

  ab_scores_c_clusters <<- read.table(AB.SCORES.C.CLUSTERS.PATH)
  repli_scores_c_clusters <<- read.table(REPLI.SCORES.C.CLUSTERS.PATH)

  load(EMB.ERY.CHROM.DATA.PATH)
  load(EMB.ESC.CHROM.DATA.PATH)
  ery_loc_SG1 <<- ery_loc_SG1
  emb_loc_SG1 <<- emb_loc_SG1
  esc_loc_SG1 <<- esc_loc_SG1

  load(EARLY.LATE.BINS.PATH)
  early_late_bins <<- early_late_bins

  if (load.heavy) {
    load(BIN.CELL.AB.HRES.PATH)
    bin_cell_ab_hres <<- bin_cell_ab_hres
    load(BIN.CELL.COV.HRES.PATH)
    bin_cell_cov_hres_ds <<- bin_cell_cov_hres_ds
    bin_cell_cov_hres <<- bin_cell_cov_hres
  }
}

write.geo.files <- function(output.dir) {
  write.table(bin_cell_cov_hres, file=file.path(output.dir, 'genomic_bin_coverage_per_cell.txt'), sep="\t", quote=F)
  write.table(bin_cell_ab_hres, file=file.path(output.dir, 'genomic_bin_a_score_per_cell.txt'), sep="\t", quote=F) 

  write.table(e10_bin_cell_cov, file=file.path(output.dir, 'genomic_bin_coverage_per_cell_pery.txt'), sep="\t", quote=F)
  write.table(e10_bin_cell_ab, file=file.path(output.dir, 'genomic_bin_a_score_per_cell_pery.txt'), sep="\t", quote=F) 

  scdb_init(EXPRESSION.MC.DIR)
  mat = scdb_mat("e9")
  write.table(data.matrix(mat@mat), sep="\t", quote=F, file=file.path(output.dir, 'UMI_counts.txt'))

  e9.mc = get.or.create(file.path(repl.data.dir, 'atlases'), function() error)$e9.exp
  write.table(e9.mc@e_gc, sep="\t", quote=F, file=file.path(output.dir, 'metacell_expression.txt'))
}


create.figs <- function() {
  # figures 1A-C are created by the single-cell Hi-C procesing pipeline
  # (https://github.com/tanaylab/schic2)
  # EDF 1 is also created by the same pipeline, except for G-H which are taken from:
  # https://doi.org/10.1038/nature23001

  load.decay.metrics()
  load.ab.and.cov()

  esc.analysis()
  ery.analysis()
  repl.model.analysis()
  emb.proper.analysis()
  ecto.meso.a.analysis()
  epigenetic.analysis()

}

get.or.create <- function(data.path, func=NULL) {
  if (file.exists(data.path)) {
    var.name = load(data.path)
    obj = get(var.name)
    return(obj)
  } else {
    obj = func()
    save(obj, file=data.path)
    return(obj)
  }
}

set.misha <- function(misha.db.dir) {
  gsetroot(misha.db.dir)
  options(shaman.mc_support=1)
  options(gmax.data.size=1e9)
}

get.fig.dir <- function() {
  dir.create(FIGURE.DIR, showWarnings=F)
  return(FIGURE.DIR)
}

get.data.file.dir <- function() {
  dir.create(DATA.DIR, showWarnings=F)
  return(DATA.DIR)
}
