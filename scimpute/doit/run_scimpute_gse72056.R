library(scImpute)

data_dir <- "./"
out_dir <- "/home/jby2/XH/scImpute/doit/mousebrain_rerun_GSE72056/"
count_path <- "/home/jby2/XH/scImpute/doit/GSE72056.csv"

labels_df <- read.csv("/home/jby2/XH/scImpute/doit/new_GSE72056_ct_.csv", stringsAsFactors = FALSE)
labels_ <- labels_df$cell_type 


if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


scimpute(count_path = count_path, infile = "csv", outfile = "csv", labeled = TRUE,
         Kcluster = NULL, labels = labels_, out_dir = out_dir, drop_thre = 0.5, ncores = 36)