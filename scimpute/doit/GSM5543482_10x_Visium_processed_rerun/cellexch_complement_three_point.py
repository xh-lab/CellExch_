import pandas as pd
import numpy as np
import seaborn as sns
from time import time
import networkx as nx
import matplotlib.pyplot as plt
import os
import psutil
itime = time()

# Modify the following parameters as needed.
folders_to_create = [
    r"/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/expression_thresholding",
    r"/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/expression_product",
    r"/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/cell_expression",
    r"/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/Three",
    r"/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/Three/TOP",
    r"/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun"]
for folder in folders_to_create:
    os.makedirs(folder, exist_ok=True)
cancer = r'/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun'  # File directory for cancer species
x_ytick = ['Medullary', 'Mesangial', 'Collecting', 'Tubule', 'Progenitor']
cell_type = 4
thr = pd.read_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/The_expression_thresholding_data.csv", header=None, index_col=None).to_numpy()  # melanoma
pro = pd.read_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/The_expression_product_data.csv", header=None, index_col=None).to_numpy()  # melanoma
cell = pd.read_csv("/home/jby2/XH/scImpute/doit/GSM5543482_10x_Visium_processed_rerun/The_cell_expression_data.csv", header=None, index_col=None).to_numpy()  # melanoma
LRI = pd.read_csv("/home/jby2/XH/Cell_Dialog/mouse_brain/LRI.csv", header=None, index_col=None).to_numpy()  # LRI obtained after filtering


# --------------------------------------------------------------------------------------------------------------
#  The expression thresholding calculation method code
for w in range(0, LRI.shape[0]):
    zhibiao1 = 0
    zhibiao2 = 0
    flag = False
    for x in range(0, thr.shape[0]):
        if LRI[w][0] == thr[x][0]:
            b11 = thr[x]
            zhibiao1 = 1
        if LRI[w][1] == thr[x][0]:
            b12 = thr[x]
            zhibiao2 = 1
        if zhibiao1 == 1 and zhibiao2 == 1:
            flag = True
            break
    if not flag:
        continue
    b11_l = b11[0]
    b11_r = b12[0]
    LRI_com = "{}_{}".format(b11_l, b11_r)
    b11 = np.delete(b11, 0)
    b12 = np.delete(b12, 0)
    for i in range(0, cell_type + 1):
        for j in range(0, cell_type + 1):
            value0 = 0
            value1 = 1
            if b11[i] == b12[j] and b11[i] == 1 and b12[j] == 1:
                with open(cancer + "/expression_thresholding/" + str(i) + str(j) + ".csv", mode="a") as f:
                    f.write("{},{}\n".format(LRI_com, value1))
                f.close()
            else:
                with open(cancer + "/expression_thresholding/" + str(i) + str(j) + ".csv", mode="a") as f:
                    f.write("{},{}\n".format(LRI_com, value0))
                f.close()

for i in range(0, cell_type + 1):
    for j in range(0, cell_type + 1):
        a1 = pd.read_csv(cancer + "/expression_thresholding/" + str(i) + str(j) + ".csv", header=None,
                         index_col=None).to_numpy()  # Obtain the result
        SUM = sum(a1[:, 1])
        with open(cancer + "/thresholding_result.csv", mode="a") as f:
            f.write(str(i))
            f.write(str(j))
            f.write('____')
            f.write(str(SUM))
            f.write('\n')
        f.close()

#  The expression product calculation method code

for w in range(0, LRI.shape[0]):
    zhibiao1 = 0
    zhibiao2 = 0
    flag = False
    for x in range(0, pro.shape[0]):
        if LRI[w][0] == pro[x][0]:
            b11 = pro[x]
            zhibiao1 = 1
        if LRI[w][1] == pro[x][0]:
            b12 = pro[x]
            zhibiao2 = 1
        if zhibiao1 == 1 and zhibiao2 == 1:
            flag = True
            break
    if not flag:
        continue
    b11_l = b11[0]
    b11_r = b12[0]
    LRI_com = "{}_{}".format(b11_l, b11_r)
    b11 = np.delete(b11, 0)
    b12 = np.delete(b12, 0)
    for i in range(0, cell_type + 1):
        for j in range(0, cell_type + 1):
            Cheng = b11[i] * b12[j]
            with open(cancer + "/expression_product/" + str(i) + str(j) + ".csv", mode="a") as f:
                f.write("{},{}\n".format(LRI_com, Cheng))
            f.close()

for i in range(0, cell_type + 1):
    for j in range(0, cell_type + 1):
        a1 = pd.read_csv(cancer + "/expression_product/" + str(i) + str(j) + ".csv", header=None, index_col=None).to_numpy()  # Obtain the result
        SUM = sum(a1[:, 1])
        with open(cancer + "/product_result.csv", mode="a") as f:
            f.write(str(i))
            f.write(str(j))
            f.write('____')
            f.write(str(SUM))
            f.write('\n')
        f.close()

# The cell expression calculation method code

for w in range(0, LRI.shape[0]):
    zhibiao1 = 0
    zhibiao2 = 0
    flag = False
    for x in range(0, cell.shape[0]):
        if LRI[w][0] == cell[x][0]:
            b11 = cell[x]
            zhibiao1 = 1
        if LRI[w][1] == cell[x][0]:
            b12 = cell[x]
            zhibiao2 = 1
        if zhibiao1 == 1 and zhibiao2 == 1:
            flag = True
            break
    if not flag:
        continue
    b11_l = b11[0]
    b11_r = b12[0]
    LRI_com = "{}_{}".format(b11_l, b11_r)
    b11 = np.delete(b11, 0)
    b12 = np.delete(b12, 0)
    for i in range(0, cell_type + 1):
        for j in range(0, cell_type + 1):
            cell_ = b11[i] * b12[j]
            with open(cancer + "/cell_expression/" + str(i) + str(j) + ".csv", mode="a") as f:
                f.write("{},{}\n".format(LRI_com, cell_))
            f.close()

for i in range(0, cell_type + 1):
    for j in range(0, cell_type + 1):
        a1 = pd.read_csv(cancer + "/cell_expression/" + str(i) + str(j) + ".csv", header=None, index_col=None).to_numpy()  # Obtain the result
        SUM = sum(a1[:, 1])
        with open(cancer + "/cell_result.csv", mode="a") as f:
            f.write(str(i))
            f.write(str(j))
            f.write('____')
            f.write(str(SUM))
            f.write('\n')
        f.close()

explain = "The xxx_result.csv file indicates that 0 represents the Melanoma cancer cells,\n1 represents the T cells," \
          "\n2 represents the B cells,\n3 represents the Macrophages,\n4 represents the Endothelial cells," \
          "\n5 represents the CAFs,\n6 represents the NK cells.\nFor example: 12_xxx represents the communication " \
          "between T cells and B cells, and xxx is the calculated communication strength. "
print('--------------------------------------------------------------')
print(explain)


# Processing data
sum_data = 0
data = pd.read_csv(cancer + "/thresholding_result.csv", header=None,index_col=None).to_numpy()
data2 = []
for j in range(len(data)):
    data11 = data[j][0]
    data1 = float(data11[6:])
    a1 = data11[:2]
    if a1[0] == a1[1]:
        data1 = 0.0
    data2.append(data1)
# normalization
_range = np.max(data2) - np.min(data2)
aaa = (data2 - np.min(data2)) / _range
sum_data = pd.DataFrame(aaa)
totol = sum_data.values
result = totol.reshape((cell_type + 1, cell_type + 1))
result1 = pd.DataFrame(result)


sum_data = 0
data = pd.read_csv(cancer + "/product_result.csv", header=None,index_col=None).to_numpy()
data2 = []
for j in range(len(data)):
    data11 = data[j][0]
    data1 = float(data11[6:])
    a1 = data11[:2]
    if a1[0] == a1[1]:
        data1 = 0.0
    data2.append(data1)
# normalization
_range = np.max(data2) - np.min(data2)
aaa = (data2 - np.min(data2)) / _range
sum_data = pd.DataFrame(aaa)
totol = sum_data.values
result = totol.reshape((cell_type + 1, cell_type + 1))
result2 = pd.DataFrame(result)


sum_data = 0
data = pd.read_csv(cancer + "/cell_result.csv", header=None,index_col=None).to_numpy()
data2 = []
for j in range(len(data)):
    data11 = data[j][0]
    data1 = float(data11[6:])
    a1 = data11[:2]
    if a1[0] == a1[1]:
        data1 = 0.0
    data2.append(data1)
# normalization
_range = np.max(data2) - np.min(data2)
aaa = (data2 - np.min(data2)) / _range
sum_data = pd.DataFrame(aaa)
totol = sum_data.values
result = totol.reshape((cell_type + 1, cell_type + 1))
result3 = pd.DataFrame(result)


# The three-point estimation method
result_max = np.maximum(result1, result2)
result_med = np.median([result1, result2, result3], axis=0)
result_min = np.minimum(result1, result2)
result_matrix = np.maximum(result_max, result3) + result_med * 4 + np.minimum(result_min, result3)
result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)


# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='BuPu',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.margins(0.05, 0.05)                       
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "/case_study_blue.csv", index=True)
plt.savefig(cancer + '/case_study_heatmap_blue.pdf', dpi=1080,bbox_inches = 'tight')
plt.close(fig)
print("-----CellExch Run Completed----")
# plt.show()

#  network view
fig = plt.figure()
df = pd.read_csv(cancer + "/case_study.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]
node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
G = nx.DiGraph()
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])
pos = nx.circular_layout(G)
edge_widths = [6 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=14, node_size=300, connectionstyle='arc3,rad=0.1')
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=250)
nx.draw_networkx_labels(G, pos, font_size=12, font_family='sans-serif')
plt.axis('off')
plt.savefig(cancer + '/case_study_networkx_blue.pdf', dpi=1080,bbox_inches = 'tight')
# plt.show()
# print(' Time: {}seconds, mem: {}MB'.format((time() - itime), psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024))