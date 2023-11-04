import argparse
import os
import re
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
from matplotlib.gridspec import GridSpec

def local_variable():
	"""
	For test purposes
	:return:
	"""
	input_dir = os.path.abspath('./demo/test_fig')

def make_fig2_comparession_kmer_conservation(df1, df2):
	"""
	Make 2 subplots for compression and kmer conservation from df1 and df2, both have index=file, cols = FMH, t1,t2..
	"""
	fig, axs = plt.subplots(1, 2, figsize=(15, 6))
	sns.boxplot(data=df1, ax=axs[0], color="orange")
	#axs[0].set(xlabel="Condition")
	#axs[0].set_title('Compare compression factor')
	axs[0].set_ylim([9.5,10.5])
	axs[0].set_title(f'({chr(97 + 0)})' + ' Compression factor (expected 10)', pad=20)
	#axs[0].axhline(y=10, color='red', linestyle='--', label='Expected value')

	sns.boxplot(data=df2, ax=axs[1], color="orange")
	axs[1].set_ylim([0.3, 0.5])
	#axs[1].set(ylabel="K-mer conservation", xlabel="Condition")
	axs[1].set_title(f'({chr(97 + 1)})' + ' K-mer conservation ratio', pad=20)
	#axs[1].set_title('Compare k-mer conservation')
	#fig.suptitle('Compare compression and k-mer conservation', fontsize=24)

	# add footer
	# add more space at bottom
	plt.subplots_adjust(top=0.9, bottom=0.15)
	# footer
	fig.text(0.5, 0.05, 'Sketching conditions', horizontalalignment="center", fontsize=20)

	fig.savefig("compression.png", dpi=200)
	plt.close(fig)

def make_fig3_ciji_comparison(df):
	"""
	Input df: collected ciji results, with columns:
	['gt_ji', 'fmh_ji', 'sync_ji', 'gt_ci', 'fmh_ci', 'sync_ci', 'fmh_rela_ji', 'os_rela_ji', 'fmh_rela_ci', 'os_rela_ci']
	"""
	fig, axs = plt.subplots(1, 3, figsize=(16, 5))
	# sub1: (ji, ji) dots
	sns.scatterplot(x='fmh_ji', y='sync_ji', data=df, color='orange', marker='o', ax=axs[0], s=50, alpha=0.5)
	# Set axis labels and title
	axs[0].set_xlabel('FMH JI')
	axs[0].set_ylabel('OS JI')
	axs[0].set_title('(a) Scatter plot for JI pairs', pad=20)
	# Plot the y=x line
	x_line = np.linspace(0, 1, 100)
	y_line = x_line
	axs[0].plot(x_line, y_line, label='Y=X Line', color='grey', linestyle='--', linewidth=1)
	#axs[0].legend()

	# sub2: (ci, ci) dots
	sns.scatterplot(x='fmh_ci', y='sync_ci', data=df, color='orange', marker='o', ax=axs[1], s=50, alpha=0.5)
	# Set axis labels and title
	axs[1].set_xlabel('FMH CI')
	axs[1].set_ylabel('OS CI')
	axs[1].set_title('(b) Scatter plot CI pairs', pad=20)
	# Plot the y=x line
	x_line = np.linspace(0, 1, 100)
	y_line = x_line
	axs[1].plot(x_line, y_line, label='Y=X Line', color='grey', linestyle='--', linewidth=1)
	#axs[1].legend()

	# sub3: boxplot of relative error in percent
	sns.boxplot(data=df.loc[:, ['fmh_rela_ji', 'os_rela_ji', 'fmh_rela_ci', 'os_rela_ci']], ax=axs[2], color="orange")
	axs[2].set(ylabel="Relative error %", xlabel="Category")
	axs[2].set_xticklabels(['FMH JI', 'OS JI', 'FMH CI', 'OS CI'])
	axs[2].set_title('(c) Compare relative error', pad=20)
	axs[2].set_ylim([-5, 5])

	# # add abc title
	# for i, ax in enumerate(axs):
	# 	ax.set_title(f'({chr(97 + i)})', pad=20)

	fig.savefig("ji_est.png", dpi=200)
	plt.close(fig)

def make_fig4_dist_dist(df1, df2):
	"""
	Df1 is the genome coverage by condition, with index of file, col as condition
	df2 is a melted (by condition) freq table, with 3 cols [distance, condition, freq_ratio]
	"""
	# make 2 plot: 1 for coverage shift; 2 for distance dist. multiple line plot
	fig, axs = plt.subplots(1, 2, figsize=(16, 5), gridspec_kw={'width_ratios': [1, 2]})
	# sub1: boxplot of genome coverage
	sns.boxplot(data=df1, ax=axs[0], color="orange")
	axs[0].set(ylabel="Ratio", xlabel="Condition")
	# rotate xticks to look better
	axs[0].tick_params(axis='x', rotation=30)
	axs[0].set_title('(a) Genome coverage', pad=20)
	axs[0].set_ylim([0.8, 1])

	# sub2: distance distribution, multiple lines
	sns.lineplot(x='distance', y='freq_ratio', data=df2, hue='condition', ax=axs[1], marker='o', linewidth=0.2)
	axs[1].set(ylabel="Frequency ratio", xlabel="Distance")
	axs[1].set_title('(b) Distance distribution', pad=20)
	axs[1].set_ylim([-0.01, 0.2])
	axs[1].set_xlim([1, 36])
	axs[1].legend(title='Category')
	axs[1].axvline(x=8, color='red', linestyle='--', label='Average distance')
	axs[1].legend()
	fig.savefig("dist_coverage_compare.png", dpi=200)
	plt.close(fig)

def make_fig5_mos(df1, df2):
	"""
	df1: a file by k table for MOS kmer conservation
	df2: a file by k table for resular OS kmer conservation
	both should have the same colnames
	"""
	col_to_use = list(df1.columns)[0:min(4, len(df1.columns))]
	num_subplot = len(col_to_use)

	fig, axs = plt.subplots(1, num_subplot, figsize=(5*num_subplot, 5))
	#fig.suptitle("K-mer conservation for multiple k sizes")

	for idx, ksize in enumerate(col_to_use):
		sns.scatterplot(x=df1[ksize], y=df2[ksize], color='orange', marker='o', ax=axs[idx], s=50, alpha=0.5)
		# Set axis labels and title
		#axs[idx].set_xlabel('Conserved k-mers by MOS')
		axs[idx].set_ylabel('')
		axs[idx].set_xlabel('')
		# find a proper ylim region to present the dots
		ymin = min(min(df1[ksize]), min(df2[ksize]))
		ymax = max(max(df1[ksize]), max(df2[ksize]))
		min_lim = round(ymin, 2)-0.01
		max_lim = round(ymax, 2)+0.01
		# force to show stick by step of 0.01
		reset_tick = []
		current_value = min_lim
		while current_value <= max_lim:
			reset_tick.append(current_value)
			current_value += 0.01
		# some subfigure need 1 more tick
		if len(reset_tick) < 5:
			reset_tick.append(reset_tick[-1]+0.01)
		axs[idx].set_ylim(min_lim, max_lim)
		axs[idx].set_xlim(min_lim, max_lim)
		axs[idx].set_xticks(reset_tick)
		axs[idx].set_yticks(reset_tick[1:])
		x_line = np.linspace(0,1, 100)
		y_line = x_line
		axs[idx].plot(x_line, y_line, label='Y=X Line', color='grey', linestyle='--', linewidth=1)
		axs[idx].set_title(f'({chr(97 + idx)})' + " " + ksize, pad=20)

	# only use 1 ylable for all
	axs[0].set_ylabel('Ratio of conserved k-mers by OS')
	# add more space at bottom
	plt.subplots_adjust(top=0.9, bottom=0.2)
	# footer
	fig.text(0.5, 0.05, 'Ratio of conserved k-mers by MOS at different lengths', horizontalalignment="center")

	fig.savefig("mos.png", dpi=200)
	plt.close(fig)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Syncmer pilot",
	                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-d', '--input_dir', type=str, help="Path to dir containing all output files", default=".")

	### read parameters
	args = parser.parse_args()
	input_dir = args.inut_dir
	os.chdir(input_dir)

	# Set a global font size for text elements in Matplotlib
	mpl.rcParams['font.size'] = 16


	### fig2: comparison and k-mer conservation with regard to shift
	df2_compression = pd.read_csv("compression_change_by_shift_k20s11.csv", index_col=0)
	df2_kmer_conser = pd.read_csv("conservation_change_by_shift_k20s11_mut0.05.csv", index_col=0)
	# already fixed in code: forgot to change colnames for conservation table, here is a patch
	old_names = list(df2_kmer_conser.columns)
	new_names = []
	if 'OS_' not in old_names[-1]:
		for i,item in enumerate(old_names):
			if item == 'FMH':
				new_names.append(item)
			else:
				new_names.append('OS_t'+str(int(item)+1))
		df2_kmer_conser.columns = new_names
	# make fig
	make_fig2_comparession_kmer_conservation(df1=df2_compression, df2=df2_kmer_conser)



	### fig3: JI CI estimation
	df3 = pd.read_csv("compare_jici.csv")
	df3['fmh_rela_ji'] = (df3['fmh_ji'] / df3['gt_ji'] - 1) * 100
	df3['os_rela_ji'] = (df3['sync_ji'] / df3['gt_ji'] - 1) * 100
	df3['fmh_rela_ci'] = (df3['fmh_ci'] / df3['gt_ci'] - 1) * 100
	df3['os_rela_ci'] = (df3['sync_ci'] / df3['gt_ci'] - 1) * 100
	make_fig3_ciji_comparison(df=df3)



	### fig4: distance dist.
	df41 = pd.read_csv("coverage_change_by_shift_k20s11.csv", index_col=0)
	# distance: get a mean line for all conditions
	temp_series_list = []
	name_list = []
	# start from FMH file
	temp_df = pd.read_csv("distance_distribution_FMH_k20c10.csv", index_col=0)
	temp_series_list.append(temp_df.mean(axis=1))
	name_list.append('FMH')
	# then all OS file
	file_list = [x for x in glob.glob("distance_distribution_*.csv") if "S_k20s11t" in x]
	file_list.sort() # make sure start from t0
	for i,file in enumerate(file_list):
		temp_df = pd.read_csv(file, index_col=0)
		temp_series_list.append(temp_df.mean(axis=1))
		name_list.append('OS_t'+str(i+1))


	df42 = pd.DataFrame(temp_series_list, index=name_list)
	df42 = df42.transpose()
	df42.fillna(0, inplace=True)
	df42['distance'] = df42.index
	# melt df for plot
	melted_df42 = pd.melt(df42, id_vars='distance', var_name='condition', value_name='freq_ratio')
	# make fig
	make_fig4_dist_dist(df1=df41, df2=melted_df42)



	### fig5: multi-resolution
	df51 = pd.read_csv("compare_conservation_MOS_s11_mut0.05.csv", index_col=0)
	df52 = pd.read_csv("compare_conservation_regular_OS_s11_mut0.05.csv", index_col=0)
	make_fig5_mos(df1=df51, df2=df52)









