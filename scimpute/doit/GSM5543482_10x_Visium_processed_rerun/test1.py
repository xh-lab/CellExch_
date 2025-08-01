import pandas as pd

df = pd.read_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/CCC_Analysis/Medullary_Mesangial.csv", header=0)

df.iloc[:, 0] = df.iloc[:, 0].str.extract(r'(Mesangial)', expand=False).fillna(df.iloc[:, 0])

output_path = "/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/CCC_Analysis/Medullary_Mesangial_modified.csv"
df.to_csv(output_path, index=False)


df = pd.read_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/CCC_Analysis/Medullary_Collecting.csv", header=0)

df.iloc[:, 0] = df.iloc[:, 0].str.extract(r'(Collecting)', expand=False).fillna(df.iloc[:, 0])

output_path = "/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/CCC_Analysis/Medullary_Collecting_modified.csv"
df.to_csv(output_path, index=False)


df = pd.read_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/CCC_Analysis/Medullary_Tubule.csv", header=0)

df.iloc[:, 0] = df.iloc[:, 0].str.extract(r'(Tubule)', expand=False).fillna(df.iloc[:, 0])

output_path = "/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/CCC_Analysis/Medullary_Tubule_modified.csv"
df.to_csv(output_path, index=False)


df = pd.read_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/CCC_Analysis/Medullary_Progenitor.csv", header=0)

df.iloc[:, 0] = df.iloc[:, 0].str.extract(r'(Progenitor)', expand=False).fillna(df.iloc[:, 0])

output_path = "/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/CCC_Analysis/Medullary_Progenitor_modified.csv"
df.to_csv(output_path, index=False)