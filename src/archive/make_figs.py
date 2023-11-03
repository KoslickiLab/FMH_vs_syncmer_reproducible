import argparse
import os
import re
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
from matplotlib.gridspec import GridSpec


# fig2 compression
file=os.path.abspath("./check_shift/compression_change_by_shift_k12s8t.csv")
df = pd.read_csv(file)
df.columns = ['name', 't=1', 't=2', 't=3']
# the FMH compression is in another file, will rerun to unify all parameters. Generate a random col here for coding
df['FMH'] = np.random.normal(0.2, 0.0002, 20)

# make multiple boxplot
fig, axs = plt.subplots(1,1, figsize=(6, 5))
sns.boxplot(data=df.iloc[:,1:], ax=axs, color="orange")
axs.set(ylabel="Ratio: 1 / compression factor", xlabel="Parameters")
axs.set_title('Compare compression factor')
axs.set_ylim([0.1, 0.3])
axs.axhline(y=0.2, color='red', linestyle='--', label='Expected value')
fig.savefig("compression.png", dpi=300)
plt.close(fig)



# fig3 compare CIJI
file=os.path.abspath("./result_random20_Brucella_K20s11ratio10_shift5/out_2_compare_jici.csv")
df = pd.read_csv(file)
# relative error in percent
df['fmh_rela_ji'] = (df['fmh_ji'] / df['gt_ji'] - 1)*100
df['os_rela_ji'] = (df['sync_ji'] / df['gt_ji'] - 1)*100
df['fmh_rela_ci'] = (df['fmh_ci'] / df['gt_ci'] - 1)*100
df['os_rela_ci'] = (df['sync_ci'] / df['gt_ci'] - 1)*100

fig, axs = plt.subplots(1,3, figsize=(15, 5))
# sub1: (ji, ji) dots
sns.scatterplot(x='fmh_ji', y='sync_ji', data=df, color='tan', marker='o', ax=axs[0])
# Set axis labels and title
axs[0].set_xlabel('FMH JI')
axs[0].set_ylabel('OS JI')
axs[0].set_title('Scatter plot for JI pairs')
# Plot the y=x line
x_line = np.linspace(0, 1, 100)
y_line = x_line
axs[0].plot(x_line, y_line, label='Y=X Line', color='grey', linestyle='--', linewidth=1)
axs[0].legend()

# sub2: (ci, ci) dots
sns.scatterplot(x='fmh_ci', y='sync_ci', data=df, color='gold', marker='o', ax=axs[1])
# Set axis labels and title
axs[1].set_xlabel('FMH CI')
axs[1].set_ylabel('OS CI')
axs[1].set_title('Scatter plot CI pairs')
# Plot the y=x line
x_line = np.linspace(0, 1, 100)
y_line = x_line
axs[1].plot(x_line, y_line, label='Y=X Line', color='grey', linestyle='--', linewidth=1)
axs[1].legend()

# sub3: boxplot of relative error in percent
sns.boxplot(data=df.loc[:,['fmh_rela_ji', 'os_rela_ji', 'fmh_rela_ci', 'os_rela_ci']], ax=axs[2], color="orange")
axs[2].set(ylabel="Relative error in percentage", xlabel="Category")
axs[2].set_xticklabels(['FMH JI', 'OS JI', 'FMH CI', 'OS CI'])
axs[2].set_title('Compare relative error')
axs[2].set_ylim([-5, 5])

# add abc title
for i, ax in enumerate(axs):
    ax.set_title(f'({chr(97 + i)})')

fig.savefig("ji_est.png", dpi=300)
plt.close(fig)



# fig 4, distance distribution and shift by t
file1 = os.path.abspath("./check_shift/coverage_change_by_shift_k12s5t.csv")
df1 = pd.read_csv(file1, index_col=0)

file_list = [x for x in glob.glob("./check_shift/*") if "distance_distribution_k12s5t" in x]
# read each and calculate average
temp_series_list=[]
for file in file_list:
 temp_df = pd.read_csv(file, index_col=0)
 temp_series_list.append(temp_df.mean(axis=1))
df2 = pd.DataFrame(temp_series_list, index=[re.sub('./check_shift/distance_distribution_k12s5', '', x) for x in file_list]).transpose()
df2.fillna(0, inplace=True)
df2 = df2.sort_index(axis=1)
df2.columns = ['t=1', 't=2', 't=3', 't=4']
df2['distance'] = df2.index
melted_df = pd.melt(df2, id_vars='distance', var_name='t value', value_name='freq_ratio')


# make 2 plot: 1 for coverage shift; 2 for distance dist. multiple line plot
fig, axs = plt.subplots(1,2, figsize=(15, 5), gridspec_kw={'width_ratios': [1, 2]})
# sub1: boxplot of genome coverage
sns.boxplot(data=df1, ax=axs[0], color="orange")
axs[0].set(ylabel="Ratio", xlabel="Category")
axs[0].set_xticklabels(["t=1", "t=2", "t=3", "t=4"])
axs[0].set_title('Base coverage by open syncmers')
axs[0].set_ylim([0.8,1])
# sub2: distance distribution, multiple lines
sns.lineplot(x='distance', y='freq_ratio', data=melted_df, hue='t value', ax=axs[1], marker='o', linewidth=0.2)
axs[1].set(ylabel="Frequency ratio", xlabel="Distance")
axs[1].set_title('Distance distribution')
axs[1].set_ylim([-0.01,0.2])
axs[1].set_xlim([1,36])
axs[1].legend(title='Category')
axs[1].axvline(x=8, color='red', linestyle='--', label='Average distance')
axs[1].legend()
fig.savefig("dist_coverage_compare.png", dpi=300)
plt.close(fig)


# fig 5: multi-resolution estimation
data_2_merge = {}
k_list = [30, 26, 22, 18, 14]
data_2_merge['p1_reg'] = [74.81006413,77.44670302,80.61172902,84.04426559,87.58865248]
data_2_merge['p5_reg'] = [22.0868278,27.4318294,33.3063701,40.9004024,50.9625127]
data_2_merge['p8_reg'] = [8.83078441,11.9633118,16.8301314,23.7323944,32.993921]
data_2_merge['p1_ext'] = [74.53658537,77.65365854,80.78536585,83.97560976,87.67691557]
data_2_merge['p5_ext'] = [22.1170732,27.3121951,33.6097561,41.3853659,51.3030747]
data_2_merge['p8_ext'] = [8.56585366,12.0731707,17.1073171,23.9219512,33.5334309]

df = pd.DataFrame(data_2_merge)
# make 3 sub fig for reg vs ext dot plot
fig, axs = plt.subplots(1,3, figsize=(15, 5))
fig.suptitle("K-mer conservation for k in [14, 18, 22, 26, 30]")
for idx, p_name in enumerate(['p1','p5','p8']):
 sns.scatterplot(x=p_name+"_reg", y=p_name+"_ext", data=df, color='tan', marker='o', ax=axs[idx])
 # Set axis labels and title
 axs[idx].set_xlabel('Conserved k-mers by Open Syncmer (%)')
 axs[idx].set_ylabel('Multi-resolution OS')
 axs[idx].set_title('Scatter plot for k-mer matches')
 # Plot the y=x line
 x_line = np.linspace(0, 100, 100)
 y_line = x_line
 axs[idx].plot(x_line, y_line, label='Y=X Line', color='grey', linestyle='--', linewidth=1)
 axs[idx].legend()
 axs[idx].set_title(f'({chr(97 + idx)})'+" mutation rate = " + re.sub('p','',p_name) + "%")

fig.savefig("mos.png", dpi=300)
plt.close(fig)


