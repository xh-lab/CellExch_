import pandas as pd


df = pd.read_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/CCC_Analysis/Medullary_Progenitor.csv", index_col=0)
stacked_df = df.stack().reset_index()
stacked_df.columns = ["Source", "Target", "Value"]
stacked_df['Source_Target'] = stacked_df['Source'] + "_" + stacked_df['Target']
stacked_df = stacked_df[['Source_Target', 'Value']]
stacked_df.to_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/CCC_Analysis/Medullary_Progenitor.csv", index=False)



df = pd.read_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/CCC_Analysis/Medullary_Tubule.csv", index_col=0)
stacked_df = df.stack().reset_index()
stacked_df.columns = ["Source", "Target", "Value"]
stacked_df['Source_Target'] = stacked_df['Source'] + "_" + stacked_df['Target']
stacked_df = stacked_df[['Source_Target', 'Value']]
stacked_df.to_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/CCC_Analysis/Medullary_Tubule.csv", index=False)


df = pd.read_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/CCC_Analysis/Medullary_Mesangial.csv", index_col=0)
stacked_df = df.stack().reset_index()
stacked_df.columns = ["Source", "Target", "Value"]
stacked_df['Source_Target'] = stacked_df['Source'] + "_" + stacked_df['Target']
stacked_df = stacked_df[['Source_Target', 'Value']]
stacked_df.to_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/CCC_Analysis/Medullary_Mesangial.csv", index=False)


df = pd.read_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/CCC_Analysis/Medullary_Collecting.csv", index_col=0)
stacked_df = df.stack().reset_index()
stacked_df.columns = ["Source", "Target", "Value"]
stacked_df['Source_Target'] = stacked_df['Source'] + "_" + stacked_df['Target']
stacked_df = stacked_df[['Source_Target', 'Value']]
stacked_df.to_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/CCC_Analysis/Medullary_Collecting.csv", index=False)