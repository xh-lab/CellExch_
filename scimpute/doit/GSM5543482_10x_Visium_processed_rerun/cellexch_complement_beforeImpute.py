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



'''
# save gene expression matrix
tmp = pd.DataFrame(gene_exp)
tmp.to_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/after_impute.csv", header=False, index=False)
'''


#--------------------------------------------end of prepreocess-------------------------------------------------

#cell_type = np.delete(single_cell[0], 0)
cell_type = single_cell.obs['cell_type']
cell_type = np.array(cell_type)
gene_exp = gene_exp[1:]
Medullary = [i for i, x in enumerate(cell_type) if x == 'Medullary Cell']
Mesangial = [i for i, x in enumerate(cell_type) if x == 'Mesangial Cell']
Collecting = [i for i, x in enumerate(cell_type) if x == 'Collecting Duct Principal Cell']
Tubule = [i for i, x in enumerate(cell_type) if x == 'Tubule Cell']
Progenitor = [i for i, x in enumerate(cell_type) if x == 'Progenitor Cell']

'''
c4 = [i for i, x in enumerate(cell_type) if x == '4']
c5 = [i for i, x in enumerate(cell_type) if x == '5']
c6 = [i for i, x in enumerate(cell_type) if x == '6']
c7 = [i for i, x in enumerate(cell_type) if x == '7']
c8 = [i for i, x in enumerate(cell_type) if x == '8']
c9 = [i for i, x in enumerate(cell_type) if x == '9']
c10 = [i for i, x in enumerate(cell_type) if x == '10']
c11 = [i for i, x in enumerate(cell_type) if x == '11']
c12 = [i for i, x in enumerate(cell_type) if x == '12']
c13 = [i for i, x in enumerate(cell_type) if x == '13']
c14 = [i for i, x in enumerate(cell_type) if x == '14']
c15 = [i for i, x in enumerate(cell_type) if x == '15']
c16 = [i for i, x in enumerate(cell_type) if x == '16']
c17 = [i for i, x in enumerate(cell_type) if x == '17']
c18 = [i for i, x in enumerate(cell_type) if x == '18']
'''

gene = gene_exp[:, 0]
gene_data = np.delete(gene_exp, 0, axis=1)

row_sums = np.sum(gene_data, axis=1)/gene_data.shape[1]
gene_data_float = gene_data.astype(float)
std = np.std(gene_data_float, axis=1)
total = row_sums + std

Medullary_means = np.mean(gene_data[:, Medullary], axis=1)
new_arr_Medullary = [1 if x > y else 0 for x, y in zip(Medullary_means, total)]
Mesangial_means = np.mean(gene_data[:, Mesangial], axis=1)
new_arr_Mesangial = [1 if x > y else 0 for x, y in zip(Mesangial_means, total)]
Collecting_means = np.mean(gene_data[:, Collecting], axis=1)
new_arr_Collecting = [1 if x > y else 0 for x, y in zip(Collecting_means, total)]
Tubule_means = np.mean(gene_data[:, Tubule], axis=1)
new_arr_Tubule = [1 if x > y else 0 for x, y in zip(Tubule_means, total)]
Progenitor_means = np.mean(gene_data[:, Progenitor], axis=1)
new_arr_Progenitor = [1 if x > y else 0 for x, y in zip(Progenitor_means, total)]

'''
c4_means = np.mean(gene_data[:, c4], axis=1)
new_arr_4 = [1 if x > y else 0 for x, y in zip(c4_means, total)]
c5_means = np.mean(gene_data[:, c5], axis=1)
new_arr_5 = [1 if x > y else 0 for x, y in zip(c5_means, total)]
c6_means = np.mean(gene_data[:, c6], axis=1)
new_arr_6 = [1 if x > y else 0 for x, y in zip(c6_means, total)]
c7_means = np.mean(gene_data[:, c7], axis=1)
new_arr_7 = [1 if x > y else 0 for x, y in zip(c7_means, total)]
c8_means = np.mean(gene_data[:, c8], axis=1)
new_arr_8 = [1 if x > y else 0 for x, y in zip(c8_means, total)]
c9_means = np.mean(gene_data[:, c9], axis=1)
new_arr_9 = [1 if x > y else 0 for x, y in zip(c9_means, total)]
c10_means = np.mean(gene_data[:, c10], axis=1)
new_arr_10 = [1 if x > y else 0 for x, y in zip(c10_means, total)]
c11_means = np.mean(gene_data[:, c11], axis=1)
new_arr_11 = [1 if x > y else 0 for x, y in zip(c11_means, total)]
c12_means = np.mean(gene_data[:, c12], axis=1)
new_arr_12 = [1 if x > y else 0 for x, y in zip(c12_means, total)]
c13_means = np.mean(gene_data[:, c13], axis=1)
new_arr_13 = [1 if x > y else 0 for x, y in zip(c13_means, total)]
c14_means = np.mean(gene_data[:, c14], axis=1)
new_arr_14 = [1 if x > y else 0 for x, y in zip(c14_means, total)]
c15_means = np.mean(gene_data[:, c15], axis=1)
new_arr_15 = [1 if x > y else 0 for x, y in zip(c15_means, total)]
c16_means = np.mean(gene_data[:, c16], axis=1)
new_arr_16 = [1 if x > y else 0 for x, y in zip(c16_means, total)]
c17_means = np.mean(gene_data[:, c17], axis=1)
new_arr_17 = [1 if x > y else 0 for x, y in zip(c17_means, total)]
c18_means = np.mean(gene_data[:, c18], axis=1)
new_arr_18 = [1 if x > y else 0 for x, y in zip(c18_means, total)]
'''

value_thre = np.append(np.expand_dims(gene, axis=1), np.expand_dims(new_arr_Medullary, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_Mesangial, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_Collecting, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_Tubule, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_Progenitor, axis=0).T, axis=1)

'''
value_thre = np.append(value_thre, np.expand_dims(new_arr_4, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_5, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_6, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_7, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_8, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_9, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_10, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_11, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_12, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_13, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_14, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_15, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_16, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_17, axis=0).T, axis=1)
value_thre = np.append(value_thre, np.expand_dims(new_arr_18, axis=0).T, axis=1)
'''

value_thre = pd.DataFrame(value_thre)
value_thre.to_csv('/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/beforeImpute/The_expression_thresholding_data.csv', index=False, header=False)
print("-----Threshold data processing completed----")

value_pro = np.append(np.expand_dims(gene, axis=1), np.expand_dims(Medullary_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(Mesangial_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(Collecting_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(Tubule_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(Progenitor_means, axis=1), axis=1)

'''
value_pro = np.append(value_pro, np.expand_dims(c4_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(c5_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(c6_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(c7_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(c8_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(c9_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(c10_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(c11_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(c12_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(c13_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(c14_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(c15_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(c16_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(c17_means, axis=1), axis=1)
value_pro = np.append(value_pro, np.expand_dims(c18_means, axis=1), axis=1)
'''

value_pro = pd.DataFrame(value_pro)
value_pro.to_csv('/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/beforeImpute/The_expression_product_data.csv', index=False, header=False)
print("-----Product data processing completed----")


Medullary_than_zero = (gene_data[:, Medullary] > 0)
cell_Medullary = np.sum(Medullary_than_zero, axis=1)/gene_data[:, Medullary].shape[1]
Mesangial_than_zero = (gene_data[:, Mesangial] > 0)
cell_Mesangial = np.sum(Mesangial_than_zero, axis=1)/gene_data[:, Mesangial].shape[1]
Collecting_than_zero = (gene_data[:, Collecting] > 0)
cell_Collecting = np.sum(Collecting_than_zero, axis=1)/gene_data[:, Collecting].shape[1]
Tubule_than_zero = (gene_data[:, Tubule] > 0)
cell_Tubule = np.sum(Tubule_than_zero, axis=1)/gene_data[:, Tubule].shape[1]
Progenitor_than_zero = (gene_data[:, Progenitor] > 0)
cell_Progenitor = np.sum(Progenitor_than_zero, axis=1)/gene_data[:, Progenitor].shape[1]

'''
c4_than_zero = (gene_data[:, c4] > 0)
cell_c4 = np.sum(c4_than_zero, axis=1)/gene_data[:, c4].shape[1]
c5_than_zero = (gene_data[:, c5] > 0)
cell_c5 = np.sum(c5_than_zero, axis=1)/gene_data[:, c5].shape[1]
c6_than_zero = (gene_data[:, c6] > 0)
cell_c6 = np.sum(c6_than_zero, axis=1)/gene_data[:, c6].shape[1]
c7_than_zero = (gene_data[:, c7] > 0)
cell_c7 = np.sum(c7_than_zero, axis=1)/gene_data[:, c7].shape[1]
c8_than_zero = (gene_data[:, c8] > 0)
cell_c8 = np.sum(c8_than_zero, axis=1)/gene_data[:, c8].shape[1]
c9_than_zero = (gene_data[:, c9] > 0)
cell_c9 = np.sum(c9_than_zero, axis=1)/gene_data[:, c9].shape[1]
c10_than_zero = (gene_data[:, c10] > 0)
cell_c10 = np.sum(c10_than_zero, axis=1)/gene_data[:, c10].shape[1]
c11_than_zero = (gene_data[:, c11] > 0)
cell_c11 = np.sum(c11_than_zero, axis=1)/gene_data[:, c11].shape[1]
c12_than_zero = (gene_data[:, c12] > 0)
cell_c12 = np.sum(c12_than_zero, axis=1)/gene_data[:, c12].shape[1]
c13_than_zero = (gene_data[:, c13] > 0)
cell_c13 = np.sum(c13_than_zero, axis=1)/gene_data[:, c13].shape[1]
c14_than_zero = (gene_data[:, c14] > 0)
cell_c14 = np.sum(c14_than_zero, axis=1)/gene_data[:, c14].shape[1]
c15_than_zero = (gene_data[:, c15] > 0)
cell_c15 = np.sum(c15_than_zero, axis=1)/gene_data[:, c15].shape[1]
c16_than_zero = (gene_data[:, c16] > 0)
cell_c16 = np.sum(c16_than_zero, axis=1)/gene_data[:, c16].shape[1]
c17_than_zero = (gene_data[:, c17] > 0)
cell_c17 = np.sum(c17_than_zero, axis=1)/gene_data[:, c17].shape[1]
c18_than_zero = (gene_data[:, c18] > 0)
cell_c18 = np.sum(c18_than_zero, axis=1)/gene_data[:, c18].shape[1]
'''

value_cell = np.append(np.expand_dims(gene, axis=1), np.expand_dims(cell_Medullary, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_Mesangial, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_Collecting, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_Tubule, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_Progenitor, axis=1), axis=1)

'''
value_cell = np.append(value_cell, np.expand_dims(cell_c4, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_c5, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_c6, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_c7, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_c8, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_c9, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_c10, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_c11, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_c12, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_c13, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_c14, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_c15, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_c16, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_c17, axis=1), axis=1)
value_cell = np.append(value_cell, np.expand_dims(cell_c18, axis=1), axis=1)
'''

value_cell = pd.DataFrame(value_cell)
value_cell.to_csv('/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/beforeImpute/The_cell_expression_data.csv', index=False, header=False)
print("-----Cell data processing completed----")

