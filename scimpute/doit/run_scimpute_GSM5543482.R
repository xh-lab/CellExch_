library(scImpute)

data_dir <- "./"
out_dir <- "/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/generator/"
count_path <- "/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/GSM5543482.csv"

labels_df <- read.csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/cell_type.csv", stringsAsFactors = FALSE)
labels_ <- labels_df$cell_type 


if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


scimpute(count_path = count_path, infile = "csv", outfile = "csv", labeled = TRUE,
         Kcluster = NULL, labels = labels_, out_dir = out_dir, drop_thre = 0.5, ncores = 36)