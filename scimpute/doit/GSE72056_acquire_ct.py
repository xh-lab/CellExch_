import pandas as pd

data = pd.read_csv("/home/jby2/XH/scImpute/doit/GSE72056.csv", header=None, index_col=0)

first_row = data.iloc[0]

first_row_df = pd.DataFrame([first_row])

first_row_df.to_csv("/home/jby2/XH/scImpute/doit/GSE72056_ct.csv", index=False, header=False)


data = pd.read_csv("/home/jby2/XH/scImpute/doit/GSE72056_ct.csv", header=None)

cell_types = ['Melanoma cancer cells', 'T cells', 'B cells', 'Macrophages', 'Endothelial cells', 'CAFs ', 'NK cells']

first_row = data.iloc[0]

first_row_replaced = [cell_types[int(value)] for value in first_row]

first_row_df = pd.DataFrame([first_row_replaced])

first_row_df.to_csv("/home/jby2/XH/scImpute/doit/GSE72056_ct_.csv", index=False, header=False)