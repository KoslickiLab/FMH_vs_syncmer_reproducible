import argparse
import os
from Bio import SeqIO
import re
import numpy as np
import random
import khmer
import pandas as pd
import copy
import math
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
	print("Reading paths from %s to a list." % input_file)
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
		temp_seq = str(fasta.seq)
		seq_list.append(temp_seq.replace('N', ''))  # remove N letters
	# merge every record to a string
	if len(seq_list) == 0:
		raise Exception("couldn't read fasta from input file, please double check!")
	else:
		# print("There are " + str(len(seq_list)) + " record(s).")
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
	index_list = list(range(1, 51))
	input_dataframe = input_dataframe.loc[input_dataframe.index.isin(index_list)]
	
	# transfer to ratio
	input_dataframe = input_dataframe / input_dataframe.sum(axis=0)
	
	return input_dataframe

# obj construction
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
		print("Filename: %s" % self.filename)
		print("K-mer length: %s" % self.k)
		print("submer length: %s" % self.s)
		print("open syncmer shift: %s" % self.shift)
		print("FMH fraction factor: %s" % self.fraction_factor)
		print("Reverse complement: %s" % str(self.rev_comp))
		print("Genome record and length: \n%s \n%s" % (len(self.seq_list), [len(x) for x in self.seq_list]))
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
		dict_kmer_dist = self.dist_kmer  # use this to record distance between 2 consecutive kmers
		max_prime = self.max_prime
		
		# for syncmers
		s_value = self.s
		shift = self.shift
		syncmer_set = self.syncmers
		dict_syncmer_dist = self.dist_syncmer  # use this to record distance between 2 consecutive kmers
		
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
		
		# adj algorithm: JI should be S(AUB) \cap S(A) \cap S(B), not just A\capB
		# ci alg is correct, S(A) \cap B / S(A), is equivalent to S(A) \cap S(B) in FMH
		kmer_intersection = len(self.kmers.intersection(other.kmers))
		kmer_ji = 1.0 * kmer_intersection / len(self.kmers.union(other.kmers))
		kmer_ci = 1.0 * kmer_intersection / len(self.kmers)
		
		syncmer_intersection = len(self.syncmers.intersection(other.syncmers))
		syncmer_ji = 1.0 * syncmer_intersection / len(self.syncmers.union(other.syncmers))
		syncmer_ci = 1.0 * syncmer_intersection / len(self.syncmers)
		
		return [round(gt_ji, 3), round(kmer_ji, 3), round(syncmer_ji, 3), round(gt_ci, 3), round(kmer_ci, 3),
		        round(syncmer_ci, 3)]
	
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

def local_variable():
	"""
	For test purposes
	:return:
	"""
	genome_list = os.path.abspath('./demo/file_paths.txt')
	kvalue = 20
	svalue = 11
	fraction_factor = 10
	shift_high = int((kvalue - svalue) / 2)
	genome_list = check_files(genome_list)
	genome_list = genome_list[0:2]

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Syncmer pilot",
	                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-f', '--input_file', type=str, help="Path to file containing genome files to check")
	parser.add_argument('-k', '--kvalue', type=int, help="k value to use", default=20)
	parser.add_argument('-s', '--svalue', type=int, help="s value to use", default=11)
	parser.add_argument('-a', '--fraction_factor', type=int, help="FMH fraction factor, ratio is 1/r", default=10)

	### read parameters
	args = parser.parse_args()
	genome_list = check_files(args.input_file)
	kvalue = args.kvalue
	svalue = args.svalue
	fraction_factor = args.fraction_factor
	if fraction_factor != kvalue-svalue+1:
		raise Exception("Invalid fraction factor")
	shift_high = int((kvalue - svalue) / 2)

		
	################################################################
	### step1: get compression and genome coverages and distance distributions
	out_df_compression = [] # store a list of DFs to merge
	out_df_coverage = []
	
	print("Step1, getting compression and seq coverage......")
	# for FMH
	max_prime = 9999999999971.  #use this to shuffle hash values one more time and fix the hash range
	hash_threshold = max_prime / fraction_factor
	# temp list for results
	list_filename = []
	list_compression = []
	list_coverage = []
	out_df_distance = [] # this if for each single FMH or t value, so put inside the loop
	for fasta in genome_list:
		# get a merged string
		list_filename.append(os.path.basename(fasta))
		temp_seq_list = generate_string_list_from_genome(fasta)
		temp_seq = "".join(temp_seq_list)
		
		# loop through seq
		temp_location = []
		for i in range(len(temp_seq) - kvalue + 1):
			current_kmer = temp_seq[i:i + kvalue]
			if khmer.hash_no_rc_murmur3(current_kmer) % max_prime < hash_threshold:
				temp_location.append(i)
		
		# check distance
		temp_location = np.array(temp_location)
		temp_distance = temp_location[1:] - temp_location[:-1]
		temp_name = os.path.basename(fasta)
		# change to freq dist.
		out_df_distance.append(pd.DataFrame(pd.Series(temp_distance).value_counts(), columns=[temp_name]))
		
		# compression ratio
		temp_compression = (len(temp_seq) - kvalue + 1) / len(temp_location) * 1.0
		list_compression.append(temp_compression)
		
		# get seq coverage
		temp_seq_cov = np.zeros(len(temp_seq), dtype=int)
		for loc in temp_location:
			temp_seq_cov[loc:loc + kvalue] = 1
		temp_coverage = np.sum(temp_seq_cov) * 1.0 / len(temp_seq)
		list_coverage.append(temp_coverage)
	
	# merge and output distance: this will not be merged with other conditions
	merged_df_dist = pd.concat(out_df_distance, axis=1)
	merged_df_dist = clean_dist_distribution_dataframe(merged_df_dist)
	merged_df_dist.to_csv("distance_distribution_FMH_k%dc%d.csv" % (kvalue, fraction_factor), header=True,
	                           index=True)
	# merge to the dict
	out_df_compression.append(pd.DataFrame(list_compression, index=list_filename, columns=['FMH']))
	out_df_coverage.append(pd.DataFrame(list_coverage, index=list_filename, columns=['FMH']))
	
	# for OS with different shift values
	for shift_value in range(shift_high + 1):
		print("shift value is " + str(shift_value))
		list_filename = []
		list_compression = []
		list_coverage = []
		out_df_distance = []  # this if for each single FMH or t value, so put inside the loop
		for fasta in genome_list:
			# get a merged string
			list_filename.append(os.path.basename(fasta))
			temp_seq_list = generate_string_list_from_genome(fasta)
			temp_seq = "".join(temp_seq_list)
			
			# loop through seq
			temp_location = []
			for i in range(len(temp_seq) - kvalue + 1):
				current_kmer = temp_seq[i:i + kvalue]
				if is_kmer_open_syncmer_by_hash_value(current_kmer, s=svalue, shift=shift_value):
					temp_location.append(i)
			
			# check distance
			temp_location = np.array(temp_location)
			temp_distance = temp_location[1:] - temp_location[:-1]
			temp_name = os.path.basename(fasta)
			# change to freq dist.
			temp_df_of_freq = pd.DataFrame(pd.Series(temp_distance).value_counts())  #, columns=[temp_name]
			# not sure what happened: if I assign colname within the pd.DataFrame command, the df becomes empty
			temp_df_of_freq.columns = [temp_name]
			out_df_distance.append(temp_df_of_freq)
			
			# compression ratio
			temp_compression = (len(temp_seq) - kvalue + 1) / len(temp_location) * 1.0
			list_compression.append(temp_compression)
			
			# get seq coverage
			temp_seq_cov = np.zeros(len(temp_seq), dtype=int)
			for loc in temp_location:
				temp_seq_cov[loc:loc + kvalue] = 1
			temp_coverage = np.sum(temp_seq_cov) * 1.0 / len(temp_seq)
			list_coverage.append(temp_coverage)
		
		# merge and output distance: this will not be merged with other conditions
		merged_df_dist = pd.concat(out_df_distance, axis=1)
		merged_df_dist = clean_dist_distribution_dataframe(merged_df_dist)
		merged_df_dist.to_csv("distance_distribution_OS_k%ds%dt%d.csv" % (kvalue, svalue, shift_value),  header=True, index=True)
		# merge to the dict
		out_df_compression.append(pd.DataFrame(list_compression, index=list_filename, columns=['OS_t'+str(shift_value+1)]))
		out_df_coverage.append(pd.DataFrame(list_coverage, index=list_filename, columns=['OS_t'+str(shift_value+1)]))
	
	# merge df
	merged_df_compression = pd.concat(out_df_compression, axis=1)
	merged_df_compression.to_csv("compression_change_by_shift_k%ds%d.csv" % (kvalue, svalue), header=True, index=True)
	
	merged_df_coverage = pd.concat(out_df_coverage, axis=1)
	merged_df_coverage.to_csv("coverage_change_by_shift_k%ds%d.csv" % (kvalue, svalue), header=True, index=True)



	################################################################
	# step2: k-mer conservation by FMH and OS with different k values
	# setup: mutation each string and check conserved k-mers
	print("Step2, get kmer conservation......")
	# for FMH
	max_prime = 9999999999971.  # use this to shuffle hash values one more time and fix the hash range
	hash_threshold = max_prime / fraction_factor
	
	for mutation_rate in [0.01, 0.05, 0.08]:
		out_df_kmer_conservation = []
		list_filename = []
		list_conservation = []
	
		# FMH
		for fasta in genome_list:
			kmer_set = set()
			mut_set = set()
			# get a merged string
			list_filename.append(os.path.basename(fasta))
			temp_seq_list = generate_string_list_from_genome(fasta)
			mutated_seq_list = mutate_strings_from_a_string_list(temp_seq_list, p=mutation_rate)
			temp_seq = "".join(temp_seq_list)
			mut_seq = "".join(mutated_seq_list)
			
			# get kmer conservation
			for i in range(len(temp_seq) - kvalue + 1):
				current_kmer = temp_seq[i:i + kvalue]
				if khmer.hash_no_rc_murmur3(current_kmer) % max_prime < hash_threshold:
					kmer_set.add(current_kmer)
			
			for i in range(len(mut_seq) - kvalue + 1):
				current_kmer = mut_seq[i:i + kvalue]
				if khmer.hash_no_rc_murmur3(current_kmer) % max_prime < hash_threshold:
					mut_set.add(current_kmer)
			
			# conservation
			temp_conservation = len(kmer_set.intersection(mut_set)) / len(kmer_set)
			list_conservation.append(temp_conservation)
		# add to out df
		out_df_kmer_conservation.append(pd.DataFrame(list_conservation, index=list_filename, columns=['FMH']))
		
		
		# OS with different t value
		for shift_value in range(shift_high + 1):
			print("shift value is " + str(shift_value))
			list_filename = []
			list_conservation = []

			
			for fasta in genome_list:
				kmer_set = set()
				mut_set = set()
				# get a merged string
				list_filename.append(os.path.basename(fasta))
				temp_seq_list = generate_string_list_from_genome(fasta)
				mutated_seq_list = mutate_strings_from_a_string_list(temp_seq_list, p=mutation_rate)
				temp_seq = "".join(temp_seq_list)
				mut_seq = "".join(mutated_seq_list)
				
				# get kmer conservation
				for i in range(len(temp_seq) - kvalue + 1):
					current_kmer = temp_seq[i:i + kvalue]
					if is_kmer_open_syncmer_by_hash_value(current_kmer, s=svalue, shift=shift_value):
						kmer_set.add(current_kmer)
				
				for i in range(len(mut_seq) - kvalue + 1):
					current_kmer = mut_seq[i:i + kvalue]
					if is_kmer_open_syncmer_by_hash_value(current_kmer, s=svalue, shift=shift_value):
						mut_set.add(current_kmer)
				
				# conservation
				temp_conservation = len(kmer_set.intersection(mut_set)) / len(kmer_set)
				list_conservation.append(temp_conservation)
			# add to out df
			out_df_kmer_conservation.append(pd.DataFrame(list_conservation, index=list_filename, columns=['OS_t'+str(shift_value+1)]))
	
		# merge df
		merged_df_conservation = pd.concat(out_df_kmer_conservation, axis=1)
		merged_df_conservation.to_csv("conservation_change_by_shift_k%ds%d_mut%s.csv" % (kvalue, svalue, mutation_rate), header=True, index=True)
		
		

	################################################################
	# step3: JI CI estimation
	print("Step3, get JI CI estimation")
	# from now on, we use the optimal t value shift_high
	object_list = []
	for file in genome_list:
		object_list.append(
			kmer_obj_for_genome(k=kvalue, s=svalue, fraction_factor=fraction_factor, filename=file, shift=shift_high))

	# get CI, JI compare
	out_list = []
	for i in range(len(object_list)):
		for j in range(len(object_list)):
			if i < j:
				out_list.append(object_list[i].get_ji_ci_with_another_obj(object_list[j]))
	
	df2 = pd.DataFrame(out_list)
	df2.columns = ['gt_ji', 'fmh_ji', 'sync_ji', 'gt_ci', 'fmh_ci', 'sync_ci']
	df2.to_csv("compare_jici.csv", sep=",", header=True, index=False)
	
	
	################################################################
	# step4: multi-resolution on k-mer conservation only
	print("Step4, MOS comparison")
	# only compare k-mer conservation
	kmax = kvalue + 10 # change to 10 later
	kvalue_list = list(range(kvalue, kmax+1, 2))
	
	# mutation variance may affect results a little
	
	# get regular OS k-mer conservation
	for mutation_rate in [0.05]:
		# OS with different k size
		out_df_kmer_conservation = []
		for temp_k in kvalue_list:
			print(temp_k)
			list_filename = []
			list_conservation = []
			
			for fasta in genome_list:
				kmer_set = set()
				mut_set = set()
				# get a merged string
				list_filename.append(os.path.basename(fasta))
				temp_seq_list = generate_string_list_from_genome(fasta)
				mutated_seq_list = mutate_strings_from_a_string_list(temp_seq_list, p=mutation_rate)
				temp_seq = "".join(temp_seq_list)
				mut_seq = "".join(mutated_seq_list)
				
				# get kmer conservation
				for i in range(len(temp_seq) - temp_k + 1):
					current_kmer = temp_seq[i:i + temp_k]
					if is_kmer_open_syncmer_by_hash_value(current_kmer, s=svalue, shift=int((kvalue - svalue) / 2)):
						kmer_set.add(current_kmer)
				
				for i in range(len(mut_seq) - temp_k + 1):
					current_kmer = mut_seq[i:i + temp_k]
					if is_kmer_open_syncmer_by_hash_value(current_kmer, s=svalue, shift=int((kvalue - svalue) / 2)):
						mut_set.add(current_kmer)
				
				# conservation
				temp_conservation = len(kmer_set.intersection(mut_set)) / len(kmer_set)
				list_conservation.append(temp_conservation)
			# add to out df
			out_df_kmer_conservation.append(
				pd.DataFrame(list_conservation, index=list_filename, columns=['k='+str(temp_k)]))
		
		# merge df for all k values
		merged_df_conservation = pd.concat(out_df_kmer_conservation, axis=1)
		merged_df_conservation.to_csv("compare_conservation_regular_OS_s%d_mut%s.csv" % (svalue, mutation_rate),
		                              header=True, index=True)

	# get MOS k-mer conservation in one run
	for mutation_rate in [0.05]:
		# build MOS on kmax
		out_df_kmer_conservation = []
		list_filename = []
		list_mos_set_original_seq = []
		list_mos_set_mutated_seq = []
		# loop through genomes
		for fasta in genome_list:
			kmer_set = set()
			mut_set = set()
			# get a merged string
			list_filename.append(os.path.basename(fasta))
			temp_seq_list = generate_string_list_from_genome(fasta)
			mutated_seq_list = mutate_strings_from_a_string_list(temp_seq_list, p=mutation_rate)
			temp_seq = "".join(temp_seq_list)
			mut_seq = "".join(mutated_seq_list)
			# generate MOS and add store in the list
			# we find kvalue-based OS and extend them to kmax
			for i in range(len(temp_seq) - kmax + 1):
				current_kmer = temp_seq[i:i + kmax]
				# note the different index for middle substring
				mid_substring = current_kmer[int((kmax-kvalue)/2):int((kmax+kvalue)/2)]
				if is_kmer_open_syncmer_by_hash_value(mid_substring, s=svalue, shift=int((kvalue - svalue) / 2)):
					kmer_set.add(current_kmer)
			
			for i in range(len(mut_seq) - kmax + 1):
				current_kmer = mut_seq[i:i + kmax]
				# note the different index for middle substring
				mid_substring = current_kmer[int((kmax - kvalue) / 2):int((kmax + kvalue) / 2)]
				if is_kmer_open_syncmer_by_hash_value(mid_substring, s=svalue, shift=int((kvalue - svalue) / 2)):
					mut_set.add(current_kmer)
			
			# store this set into the list for future comparison
			list_mos_set_original_seq.append(kmer_set)
			list_mos_set_mutated_seq.append(mut_set)
		
		# compare conservation at multiple k values
		out_df_kmer_conservation = []
		for temp_k in kvalue_list:
			print(temp_k)
			list_conservation = []
			for i in range(len(list_mos_set_original_seq)):
				temp_kmer_set = set([x[int((kmax-temp_k)/2):int((kmax+temp_k)/2)] for x in list_mos_set_original_seq[i]])
				temp_mut_set = set([x[int((kmax-temp_k)/2):int((kmax+temp_k)/2)] for x in list_mos_set_mutated_seq[i]])
				list_conservation.append(len(temp_kmer_set.intersection(temp_mut_set)) / len(temp_kmer_set))
			# add to out df
			out_df_kmer_conservation.append(
				pd.DataFrame(list_conservation, index=list_filename, columns=['k=' + str(temp_k)]))
			
		# merge df for all k values
		merged_df_conservation = pd.concat(out_df_kmer_conservation, axis=1)
		merged_df_conservation.to_csv(
			"compare_conservation_MOS_s%d_mut%s.csv" % (svalue, mutation_rate),
			header=True, index=True)
		
		
		
		
		
	
	
	
