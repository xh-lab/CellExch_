import liana
# needed for visualization and toy data
import scanpy as sc
filename = "/home/jby2/XH/CellExch/mouse_brain/GSM5543482_10x_Visium_processed.h5ad"
adata = sc.read(filename)
# import all individual methods
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean
# run cellphonedb
cellphonedb(adata, groupby='cell_type', expr_prop=0.1, resource_name='mouseconsensus', verbose=True, key_added='cpdb_res',use_raw=False)
adata.uns['cpdb_res'].sort_values(by='cellphone_pvals').to_csv("/home/jby2/XH/CellExch/mouse_brain/other_tools/cpdb_result.csv")
# run singlecellsignalr
singlecellsignalr(adata, groupby='cell_type', expr_prop=0.1, resource_name='mouseconsensus', verbose=True, key_added='singlecellsignalr_res',use_raw=False)
adata.uns['singlecellsignalr_res'].sort_values(by='lrscore').to_csv("/home/jby2/XH/CellExch/mouse_brain/other_tools/singlecellsignalr_result.csv")

# run connectome
connectome(adata, groupby='cell_type', expr_prop=0.1, resource_name='mouseconsensus', verbose=True, key_added='connectome_res',use_raw=False)
adata.uns['connectome_res'].sort_values(by='scaled_weight').to_csv("/home/jby2/XH/CellExch/mouse_brain/other_tools/connectome_result.csv")

# run natmi
natmi(adata, groupby='cell_type', expr_prop=0.1, resource_name='mouseconsensus', verbose=True, key_added='natmi_res',use_raw=False)
adata.uns['natmi_res'].sort_values(by='spec_weight').to_csv("/home/jby2/XH/CellExch/mouse_brain/other_tools/natmi_result.csv")

# run logfc
logfc(adata, groupby='cell_type', expr_prop=0.1, resource_name='mouseconsensus', verbose=True, key_added='logfc_res',use_raw=False)
adata.uns['logfc_res'].sort_values(by='lr_logfc').to_csv("/home/jby2/XH/CellExch/mouse_brain/other_tools/logfc_result.csv")

# run geometric_mean
geometric_mean(adata, groupby='cell_type', expr_prop=0.1, resource_name='mouseconsensus', verbose=True, key_added='geometric_mean_res',use_raw=False)
adata.uns['geometric_mean_res'].sort_values(by='gmean_pvals').to_csv("/home/jby2/XH/CellExch/mouse_brain/other_tools/geometric_mean_result.csv")
# run cellchat
cellchat(adata, groupby='cell_type', expr_prop=0.1, resource_name='mouseconsensus', verbose=True, key_added='cellchat_res',use_raw=False)
adata.uns['cellchat_res'].sort_values(by='cellchat_pvals').to_csv("/home/jby2/XH/CellExch/mouse_brain/other_tools/cellchat_result.csv")