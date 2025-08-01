import pandas as pd
import numpy as np
import scanpy as sc

# Processing  scRNA-seq data
# single_cell1 = pd.read_csv("/home/jby2/XH/Cell_Dialog/GSE72056_1.csv", header=None, index_col=None)  # scRNA-seq data
# single_cell2 = pd.read_csv("/home/jby2/XH/Cell_Dialog/GSE72056_2.csv", header=None, index_col=None)
# single_cell2 = single_cell2.drop(0)
# single_cell = pd.concat([single_cell1, single_cell2])
# single_cell.to_csv('/home/jby2/XH/Cell_Dialog/GSE72056.csv', index=False, header=False)

single_cell = sc.read_h5ad("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/GSM5543482_10x_Visium_processed.h5ad")

# acquire gene expression matrix
gene_exp = single_cell.X.toarray()

cell_name = np.array(single_cell.obs_names)
cell_name = np.insert(cell_name, 0, '/').reshape(-1, 1)
gene_name = np.array(single_cell.var_names)

gene_exp = np.vstack((gene_name, gene_exp))
gene_exp = np.hstack((cell_name, gene_exp))
gene_exp = gene_exp.T  # row represents gene, col represents cell

cell_type = single_cell.obs['cell_type']
cell_type = np.array(cell_type)
cell_type = np.expand_dims(cell_type, axis=0) 
cell_type_padded = np.insert(cell_type, 0, '/', axis=1)

gene_exp[0] = cell_type_padded


# save gene expression matrix
tmp = pd.DataFrame(gene_exp)
tmp.to_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/GSM5543482.csv", header=False, index=False)