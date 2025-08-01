import pandas as pd


file1 = "/home/jby2/XH/scImpute/doit/GSE72056_ct_.csv"
df = pd.read_csv(file1, header=None)  


data = df.iloc[0].tolist()
new_df = pd.DataFrame(data, columns=['cell_type'])


new_df.to_csv("/home/jby2/XH/scImpute/doit/new_GSE72056_ct_.csv", index=False)