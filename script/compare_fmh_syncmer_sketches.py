import argparse
import os
from Bio import SeqIO
import gzip
import re
import numpy as np
import random
import khmer
import pandas as pd
import copy
import math
import subprocess
# import screed
import bisect
import inspect
import gc
import statistics
import seaborn as sns
import matplotlib.pyplot as plt


# ----------------------------------------------------------------
# define a syncmer
def is_kmer_open_syncmer_by_hash_value(kmer, s, shift=0, max_prime=9999999999971.):
	"""
	Determines whether a kmer is an open syncmer with position shift
	:param kmer: input kmer
	:param s: submer length
	:param shift: location of the min submer (default 0, which is the start)
	:return: boolean T/F
	"""
	target_submer = khmer.hash_no_rc_murmur3(kmer[shift:shift + s]) % max_prime
	for start_pos in list(range(len(kmer) - s + 1)):
		if target_submer > khmer.hash_no_rc_murmur3(kmer[start_pos:start_pos + s]) % max_prime:
			return False
	
	return True

# read input file containing paths of fasta files
def check_files(input_file):
	"""
	Read all file paths from an input file and check if they exist
	:param input_file:
	:return: a sorted list of file paths in the input file
	"""
	print("Reading paths from %s to a list." %input_file)
	out_list = list()
	with open(input_file, 'r') as temp:
		for line in temp.readlines():
			line = line.strip()
			if not os.path.exists(line):
				raise Exception("Input file %s does not exist." % line)
			out_list.append(os.path.abspath(line))
		out_list = sorted(out_list, key=os.path.basename)
		return (out_list)

# transfer a fasta to a list of seqs
def generate_string_list_from_genome(input_genome):
	"""
	Transfer input genome to a string
	"""
	# store fasta to a list (in case multiple records line)
	fasta_sequences = SeqIO.parse(input_genome, 'fasta')
	seq_list = []
	for fasta in fasta_sequences:
		seq_list.append(str(fasta.seq))
	# merge every record to a string
	if len(seq_list) == 0:
		raise Exception("couldn't read fasta from input file, please double check!")
	else:
		#print("There are " + str(len(seq_list)) + " record(s).")
		return seq_list

# mutate seq list by independent point mutation method
def mutate_strings_from_a_string_list(input_list: list, p: float):
	"""
	Mutate a string dict with mutation rate p

	:param input_list: a string list generate from previous helper function
	:param p: mutation rate p in [0, 1]. IID in every loci
	:return: a mutated string dict
	"""
	out_list = []

	for input_string in input_list:

		# simple mutation event by Unif(0,1)
		string_length = len(input_string)
		random_unif1 = list(np.random.uniform(size=string_length))
		mutation_id = [x < p for x in random_unif1]  # use collections.Counter(mutation_id) to check number of elements

		# generate mutate string
		out_string = ""
		counter = 0
		while counter < string_length:
			if not mutation_id[counter]:
				# no mutation
				out_string += input_string[counter]
			else:
				# choose from insertion / deletion / substitution of 1/3 prob
				mut_type = random.choice(['ins', 'del', 'sub'])
				if mut_type == 'ins':
					# inserted letter
					out_string += random.choice('ACGT')
					# original letter
					out_string += input_string[counter]
				elif mut_type == 'del':
					# just don't add anything
					continue
				else:
					# substitution
					dest_space = "ACGT".replace(input_string[counter], "")
					out_string += random.choice(dest_space)
			counter += 1

		# replace the string
		out_list.append(out_string)

	return out_list

# clean a merged dist distribution dataframe
def clean_dist_distribution_dataframe(input_dataframe):
	"""
	Clean the merged distance distribution dataframe for each individual files. The input df looks like:
	index: distance 0~xxx; column: filename.
	Sort, fill NA with 0, transfer number to ratio
	"""
	input_dataframe.index = input_dataframe.index.astype(int)
	input_dataframe = input_dataframe.sort_index()
	input_dataframe.fillna(0, inplace=True)
	
	# hard trim: remove 0, and remove all > 50
	index_list = list(range(1,51))
	input_dataframe = input_dataframe.loc[input_dataframe.index.isin(index_list)]
	
	# transfer to ratio
	input_dataframe = input_dataframe / input_dataframe.sum(axis=0)
	
	return input_dataframe


# ----------------------------------------------------------------
# alphabet-based syncmer is not appropriate (seq not random), below is to identify such syncmer and gives a counter example for window guarantees
def is_kmer_open_syncmer(kmer, s, shift=0):
	"""
	Determines whether a kmer is an open syncmer with position shift
	:param kmer: input kmer
	:param s: submer length
	:param shift: location of the min submer (default 0, which is the start)
	:return: boolean T/F
	"""
	target_submer = kmer[shift:shift+s]
	for start_pos in list(range(len(kmer)-s+1)):
		if target_submer > kmer[start_pos:start_pos+s]:
			return False
	
	return True

def get_syncmer_from_a_string(input_string, k, s):
	"""
	Get an input string, and find out if there exist a syncmer by alphabet order
	:param input_string:
	:param k: k-mer length
	:param s: sub-mer length
	:return:
	"""
	
	for i in range(len(input_string)-k+1):
		current_kmer = input_string[i:i+k]
		if is_kmer_open_syncmer(current_kmer, s=s, shift=0):
			print(i)
			print(current_kmer)

def example_long_string_without_syncmer():
	"""
	A manual case for weak window guarantee of open syncmer
	:return: a string that contains no open syncmer
	"""
	demo_string = "TTTTTGTTTTGTTTGTTGTGTTGGTTTGGTTGGTGTGGGTTGGGTGGGGG"

	# repeat the demo_string for TG - GC - CA, each of 2 will share 5 letters
	str2 = ["G" if x == "T" else "C" for x in demo_string]
	str3 = ["C" if x == "T" else "A" for x in demo_string]

	ultra_fake_str = demo_string + "".join(str2[5:]) + "".join(str3[5:])

	len(demo_string)
	get_syncmer_from_a_string(demo_string, k=10, s=5)
	get_syncmer_from_a_string(ultra_fake_str, 10, 5)
	len(ultra_fake_str)

	return ultra_fake_str


# ----------------------------------------------------------------
# build a submer object for downstream analysis
class kmer_obj_for_genome():
	def __init__(self, k, s, fraction_factor, filename, shift=0, max_prime=9999999999971., rev_comp=False):
		self.k = k
		self.s = s
		self.fraction_factor = fraction_factor
		self.shift = shift
		self.filename = filename
		self.max_prime = max_prime
		self.seq_list = generate_string_list_from_genome(filename)
		self.rev_comp = rev_comp
		self.cardinality = 0
		
		# kmer set
		self.all_kmers = set()
		self.kmers = set()
		self.syncmers = set()
		self.weight_kmer = dict()
		self.weight_syncmer = dict()
		
		# dist between 2 consecutive k-mers
		self.dist_kmer = dict()
		self.dist_syncmer = dict()
		
		# build sketch
		self.get_infor()
		self.build_kmer_sketch_and_count_dist()
		
	def get_infor(self):
		"""
		Print information of this object
		:return:
		"""
		print("Information of this object: \n\n")
		print("Filename: %s" %self.filename)
		print("K-mer length: %s" % self.k)
		print("submer length: %s" % self.s)
		print("open syncmer shift: %s" % self.shift)
		print("FMH fraction factor: %s" % self.fraction_factor)
		print("Reverse complement: %s" % str(self.rev_comp))
		print("Genome record and length: \n%s \n%s" %(len(self.seq_list), [len(x) for x in self.seq_list]))
		print("\n\n")
		
	def build_kmer_sketch_and_count_dist(self):
		"""
		Build the fracminhash kmer sketch and count the dist between 2 consecutive kmers
		:return:
		"""
		print("Building kmer sketch and count dist distribution")
		
		rev_comp = self.rev_comp
		all_kmers = self.all_kmers
		dict_fmh = self.weight_kmer
		dict_syncmer = self.weight_syncmer
		
		# for FMH kmers
		hash_threshold = self.max_prime / self.fraction_factor
		k_value = self.k
		kmer_set = self.kmers
		dict_kmer_dist = self.dist_kmer # use this to record distance between 2 consecutive kmers
		max_prime = self.max_prime
		
		# for syncmers
		s_value = self.s
		shift = self.shift
		syncmer_set = self.syncmers
		dict_syncmer_dist = self.dist_syncmer # use this to record distance between 2 consecutive kmers
		
		
		for seq in self.seq_list:
			# use this pointer to track distance between 2 consecutive kmers
			last_kmer_pos = 0
			last_syncmer_pos = 0
			
			# loop through each seq in the list
			for i in range(len(seq) - k_value + 1):
				kmer = seq[i:i + k_value]
				
				# rev_comp
				if rev_comp:
					kmer = min(kmer, khmer.reverse_complement(kmer))
					
				# for cardinality
				all_kmers.add(kmer)

				# build FMH kmer sketch
				if khmer.hash_no_rc_murmur3(kmer) % max_prime < hash_threshold:
					kmer_set.add(kmer)
					# track distance
					current_distance = i - last_kmer_pos
					last_kmer_pos = i
					if str(current_distance) in dict_kmer_dist:
						dict_kmer_dist[str(current_distance)] += 1
					else:
						dict_kmer_dist[str(current_distance)] = 1
					# add to weight count
					if kmer in dict_fmh:
						dict_fmh[kmer] += 1
					else:
						dict_fmh[kmer] = 1
						
				# build open syncmer sketch
				if is_kmer_open_syncmer_by_hash_value(kmer=kmer, s=s_value, shift=shift):
					syncmer_set.add(kmer)
					# track distance
					cur_sync_dist = str(i - last_syncmer_pos)
					last_syncmer_pos = i
					if cur_sync_dist in dict_syncmer_dist:
						dict_syncmer_dist[cur_sync_dist] += 1
					else:
						dict_syncmer_dist[cur_sync_dist] = 1
					# add to weight count
					if kmer in dict_syncmer:
						dict_syncmer[kmer] += 1
					else:
						dict_syncmer[kmer] = 1
		
		# update cardinality
		self.cardinality = len(all_kmers)
		
	def get_ji_ci_with_another_obj(self, other):
		"""
		Calculate ground truth JI, FMH JI and syncmer JI with another obj
		:return: a list of 3 JI values and 3 CI values
		"""
		gt_intersection = len(self.all_kmers.intersection(other.all_kmers))
		gt_ji = 1.0 * gt_intersection / len(self.all_kmers.union(other.all_kmers))
		gt_ci = 1.0 * gt_intersection / len(self.all_kmers)
		
		kmer_intersection = len(self.kmers.intersection(other.kmers))
		kmer_ji = 1.0 * kmer_intersection / len(self.kmers.union(other.kmers))
		kmer_ci = 1.0 * kmer_intersection / len(self.kmers)
		
		syncmer_intersection = len(self.syncmers.intersection(other.syncmers))
		syncmer_ji = 1.0 * syncmer_intersection / len(self.syncmers.union(other.syncmers))
		syncmer_ci = 1.0 * syncmer_intersection / len(self.syncmers)
		
		return [round(gt_ji, 3), round(kmer_ji,3), round(syncmer_ji, 3), round(gt_ci, 3), round(kmer_ci,3), round(syncmer_ci, 3)]
	
	def generate_weight_comparison(self):
		"""
		Get a weight profile for syncmer/fmh sketches of the given object
		"""
		dict_fmh = self.weight_kmer
		dict_sync = self.weight_syncmer
		
		freq_fmh = pd.DataFrame(pd.Series(dict_fmh.values()).value_counts(), columns=['fmh'])
		freq_sync = pd.DataFrame(pd.Series(dict_sync.values()).value_counts(), columns=['sync'])
		
		out_df = pd.concat([freq_fmh, freq_sync], axis=1)
		
		return out_df
		
		
	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Syncmer pilot",
	                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-f', '--input_file', type=str, help="Path to file containing genome files to check")
	parser.add_argument('-k', '--kvalue', type=int, help="k value to use", default=12)
	parser.add_argument('-s', '--svalue', type=int, help="s value to use", default=8)
	parser.add_argument('-a', '--fraction_factor', type=int, help="FMH fraction factor, ratio is 1/r", default=5)
	parser.add_argument('-r', '--rev_comp', type=str, help="Use canonical kmer", default='False')
	
	### read parameters
	args = parser.parse_args()
	genome_list = check_files(args.input_file)
	
	kvalue = args.kvalue
	svalue = args.svalue
	fraction_factor = args.fraction_factor
	
	rev_comp = args.rev_comp
	rev_comp = rev_comp == 'True'
	if rev_comp:
		print("Using canonical kmers!")
		
		
	
	###### generate object list
	object_list = []
	for file in genome_list:
		object_list.append(kmer_obj_for_genome(k=kvalue, s=svalue, fraction_factor=fraction_factor, filename=file))
	
	
	
	###### step1: compare compression ratio
	filename = []
	cardinality = []
	kmer_sketch = [] # sketch size
	syncmer_sketch = [] # sketch size
	
	for obj in object_list:
		filename.append(os.path.basename(obj.filename))
		cardinality.append(obj.cardinality)
		kmer_sketch.append(len(obj.kmers))
		syncmer_sketch.append(len(obj.syncmers))
		
	df1 = pd.DataFrame({'filename':filename, 'cardinality':cardinality, 'kmer_sketch': kmer_sketch, 'syncmer_sketch':syncmer_sketch})
	df1['compression_kmer'] = df1['kmer_sketch'] / df1['cardinality']
	df1['compression_syncmer'] = df1['syncmer_sketch'] / df1['cardinality']
	df1.to_csv("out_1_compare_compression_ratio.csv", sep=",", header=True, index=False)
	
	
	
	###### step2: compare pairwise JI/CIs
	out_list = []
	for i in range(len(object_list)):
		for j in range(len(object_list)):
			if i < j:
				out_list.append(object_list[i].get_ji_ci_with_another_obj(object_list[j]))
				
	df2 = pd.DataFrame(out_list)
	df2.columns = ['gt_ji', 'fmh_ji', 'sync_ji', 'gt_ci', 'fmh_ci', 'sync_ci']
	df2.to_csv("out_2_compare_jici.csv", sep=",", header=True, index=False)
	
	
	
	###### step3: compare dist distribution by mean frequency. May have some 0 records (beginning of contig, just ignore it here)
	
	# temp new
	temp_list_of_df_dist_freq_kmer = []
	temp_list_of_df_dist_freq_syncmer = []
	
	# get all values for box plot
	for obj in object_list:
		temp_list_of_df_dist_freq_kmer.append(pd.DataFrame.from_dict(obj.dist_kmer, orient='index', columns=[os.path.basename(obj.filename)]))
		temp_list_of_df_dist_freq_syncmer.append(pd.DataFrame.from_dict(obj.dist_syncmer, orient='index', columns=[os.path.basename(obj.filename)]))
		
	# merge df
	out_df_dist_kmer = pd.concat(temp_list_of_df_dist_freq_kmer, axis=1)
	out_df_dist_kmer = clean_dist_distribution_dataframe(out_df_dist_kmer)
	out_df_dist_kmer.to_csv("out_3_dist_distribution_kmer.csv", header=True, index=True)
	
	out_df_dist_syncmer = pd.concat(temp_list_of_df_dist_freq_syncmer, axis=1)
	out_df_dist_syncmer = clean_dist_distribution_dataframe(out_df_dist_syncmer)
	out_df_dist_syncmer.to_csv("out_3_dist_distribution_syncmer.csv", header=True, index=True)
	
	
	
	# # old code
	# dist_xaxis = list(range(1, 52))  # 1~51
	#
	# kmer_dist_freq = [0] * 51
	# syncmer_dist_freq = [0] * 51
	#
	# for obj in object_list:
	# 	# denominator for frequency
	# 	total_dist_count = sum(obj.dist_kmer.values())
	# 	total_syncmer_dist_count = sum(obj.dist_syncmer.values())
	#
	# 	for i in range(1,51): # note we record all >50 into last column
	# 		if str(i) in obj.dist_kmer:
	# 			kmer_dist_freq[i-1] += obj.dist_kmer[str(i)]/total_dist_count
	# 		if str(i) in obj.dist_syncmer:
	# 			syncmer_dist_freq[i-1] += obj.dist_syncmer[str(i)]/total_syncmer_dist_count
	#
	# # normalize (as we just sum in the loop)
	# kmer_dist_freq = [x/len(object_list) for x in kmer_dist_freq]
	# kmer_dist_freq[-1] = 1-sum(kmer_dist_freq)
	#
	# syncmer_dist_freq = [x/len(object_list) for x in syncmer_dist_freq]
	# syncmer_dist_freq[-1] = 1-sum(syncmer_dist_freq)
	#
	# # make plot
	# fig, axs = plt.subplots(1, 2, figsize=(16, 8))
	# sns.lineplot(x=dist_xaxis, y=kmer_dist_freq, ax=axs[0], color="orange")
	# axs[0].set_title("FMH kmer dist")
	# axs[0].set_ylim(0,0.22)
	# sns.lineplot(x=dist_xaxis, y=syncmer_dist_freq, ax=axs[1], color="black")
	# axs[1].set_title("Syncmer dist")
	# axs[1].set_ylim(0, 0.22)
	# fig.savefig("compare_dist_distribution.png", dpi=300)
	# plt.close(fig)



	### step4: check k-mer conservation between k-mer and syncmer
	test_file = genome_list[0]

	original_obj = kmer_obj_for_genome(k=kvalue, s=svalue, fraction_factor=fraction_factor, filename=test_file)
	seq_list = generate_string_list_from_genome(test_file)

	mut_rate = []
	kmer_conservation = []
	syncmer_conservation = []

	for mutation_rate in [x/100 for x in range(1, 20, 2)]:
		mut_rate.append(mutation_rate)
		mutated_seq_list = mutate_strings_from_a_string_list(seq_list, p=mutation_rate)
		temp_obj = copy.deepcopy(original_obj)
		# rebuild kmers with mutated sequences
		temp_obj.kmers = set()
		temp_obj.syncmers = set()
		temp_obj.all_kmers = set()
		temp_obj.seq_list = []

		temp_obj.seq_list = mutated_seq_list
		temp_obj.build_kmer_sketch_and_count_dist()

		kmer_conservation.append(len(original_obj.kmers.intersection(temp_obj.kmers)) / len(original_obj.kmers))
		syncmer_conservation.append(len(original_obj.syncmers.intersection(temp_obj.syncmers)) / len(original_obj.syncmers))

	df4 = pd.DataFrame({'mut_rate': mut_rate,  'kmer_conservation': kmer_conservation, 'syncmer_conservation': syncmer_conservation})
	df4.to_csv("out_4_compare_conservation_level.csv", sep=",", header=True, index=False)


	### step 5: generate weight comparison
	for obj in object_list:
		temp_df = obj.generate_weight_comparison()
		temp_name = re.sub('_genomic.fna.gz', '', os.path.basename(obj.filename))
		temp_df.to_csv("out_5_"+ temp_name + "_compare_kmer_weight.csv", sep=",", header=True, index=True)


		
		


	








