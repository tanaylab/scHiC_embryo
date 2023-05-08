
SCHIC.MISHA.PATH = '/path/to/mm9/'
SHAMAN.MISHA.PATH = '/path/to/work_mm9/'
ATAC.MISHA.PATH = '/path/to/mm10/trackdb'

DATA.DIR = '/path/to/downloaded/data/files/'
FIGURE.DIR = '/path/to/fig/directory/'

EXPRESSION.MC.DIR = '/path/to/metacell/directory/'

CROSS.VAL.DIR = file.path(DATA.DIR, 'embryo_cross_val')
STABILITY.DIR = file.path(DATA.DIR, 'embryo_stability')

ERYS.PATH = file.path(DATA.DIR, 'erys.txt')

DECAY.METRICS.PATH = file.path(DATA.DIR, 'sch_decay_metrics')
ESC.DECAY.METRICS.PATH = file.path(DATA.DIR, 'esc_sch_decay_metrics')
E10.DECAY.METRICS.PATH = file.path(DATA.DIR, 'charlie_decay_metrics')

BIN.CELL.COV.PATH = file.path(DATA.DIR, 'bin_cell_cov.Rda')
BIN.CELL.AB.PATH = file.path(DATA.DIR, 'base.bin_ab.Rda')
BIN.CELL.AB.HRES.PATH = file.path(DATA.DIR, '40k.bin_ab.Rda')
BIN.CELL.COV.HRES.PATH = file.path(DATA.DIR, 'bin_cov_40k.Rda')
BIN.CELL.AB.TRANS.PATH = file.path(DATA.DIR, 'bin_cell_ab_trans')
EMB.ERY.CHROM.DATA.PATH = file.path(DATA.DIR, 'emb_ery_chrom_data')
EMB.ESC.CHROM.DATA.PATH = file.path(DATA.DIR, 'emb_esc_chrom_data')
EMB.G1.CELLS.PATH = file.path(DATA.DIR, 'emb_g1_cells')

AB.SCORES.C.CLUSTERS.PATH = file.path(DATA.DIR, 'ab_score_C_clusters_40k')
REPLI.SCORES.C.CLUSTERS.PATH = file.path(DATA.DIR, 'repli_score_C_clusters_40k')
EXP.PER.BIN.PATH = file.path(DATA.DIR, 'exp_per_bin')

E10.BIN.CELL.AB.PATH = file.path(DATA.DIR, 'e10.bin_ab.Rda')
E10.BIN.CELL.COV.PATH = file.path(DATA.DIR, 'bin_cell_cov_e10.Rda')

EARLY.LATE.BINS.PATH = file.path(DATA.DIR, 'all_rep.early_late_bins.Rda')

EMB.POOL.TRACK.NAME = 'all_emb_pooled_paper'
ERY.POOL.TRACK.NAME = 'all_erys_pooled_paper'
ESC.POOL.TRACK.NAME = 'all_esc_pooled_paper'
EMB.CLS.ANALYSIS.TRACK.NAMES = c("cluster_ds_paper_1_contacts", "cluster_ds_paper_2_contacts", 
                                 "ery_track_paper", "esc_track_paper")

E10.POOL.TRACK.NAMES = c('all_charlie_pooled_paper_ery_comp_equal_size_matshuff_scores', 'all_erys_pooled_paper_ery_comp_equal_size_matshuff_scores')

CHROM.SIZE.PATH = file.path(DATA.DIR, 'chrom_sizes.txt')

# metacell
ATLAS.PATH = file.path(EXPRESSION.MC.DIR, 'cell_atlas.Rda')
E9.EXP.PATH = file.path(EXPRESSION.MC.DIR, 'e9_mc.Rda')
ESC.EXP.PATH = file.path(EXPRESSION.MC.DIR, 'esc_mc.Rda')

# embryo proper model
MODEL.PATH = file.path(DATA.DIR, 'embryo_cls_complete')
CLUSTERING.PATH = file.path(DATA.DIR, 'full_clustering')
CLUSTERING5.PATH = file.path(DATA.DIR, 'orig_clusters_c22_splitted')

AB.SCORES.EMB.CLUSTERS.SPLITTED.PATH = file.path(DATA.DIR, 'ab_score_emb_clusters_40k_splitted')
REPLI.SCORES.EMB.CLUSTERS.SPLITTED.PATH = file.path(DATA.DIR, 'repli_score_emb_clusters_40k_splitted')
REPLI.SCORES.EMB.CLUSTERS.LOW.RES.PATH = file.path(DATA.DIR, 'repli_score_emb_clusters')

MESO.ECTO.EXP = file.path(DATA.DIR, 'exp_meso_ecto')

MODEL.CLUSTER.COLORS = brewer.pal(5, 'Set1')[3:5]
CELL.CYCLE.COLS = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33')

# atac parameters
ATAC.DATA.DIR = file.path(DATA.DIR, 'atac')
ATAC.FIGURE.DIR = file.path(FIGURE.DIR, 'atac')

SCHICLUSTER.DIR = file.path(DATA.DIR, 'schicluster_files/')
SCHIC.TOPIC.MODEL.DIR = file.path(DATA.DIR, 'schic-topic-model')

TMP.CREATE.TRACK.DIR = ''
