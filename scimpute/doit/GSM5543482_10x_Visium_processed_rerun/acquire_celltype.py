import numpy as np
import pandas as pd
import scanpy as sc

single_cell = sc.read_h5ad("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/GSM5543482_10x_Visium_processed.h5ad")

cell_type = single_cell.obs['cell_type']
cell_type = np.array(cell_type)

df = pd.DataFrame(cell_type, columns=['cell_type'])
df.to_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/cell_type.csv", index=False)