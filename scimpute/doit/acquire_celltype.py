import pandas as pd
import numpy as np
import scanpy as sc

# Processing  scRNA-seq data
# single_cell1 = pd.read_csv("/home/jby2/XH/Cell_Dialog/GSE72056_1.csv", header=None, index_col=None)  # scRNA-seq data
# single_cell2 = pd.read_csv("/home/jby2/XH/Cell_Dialog/GSE72056_2.csv", header=None, index_col=None)
# single_cell2 = single_cell2.drop(0)
# single_cell = pd.concat([single_cell1, single_cell2])
# single_cell.to_csv('/home/jby2/XH/Cell_Dialog/GSE72056.csv', index=False, header=False)

single_cell = sc.read_h5ad("/home/jby2/XH/Cell_Dialog/mouse_brain/GSM5543482_10x_Visium_processed.h5ad")


#cell_type = np.delete(single_cell[0], 0)
cell_type = single_cell.obs['cell_type']
cell_type = np.array(cell_type)

cell_type_df = pd.DataFrame(cell_type, columns=["cell_type"])

cell_type_df.to_csv('/home/jby2/XH/scImpute/doit/mouse_brain_cell_type.csv', index=False)
