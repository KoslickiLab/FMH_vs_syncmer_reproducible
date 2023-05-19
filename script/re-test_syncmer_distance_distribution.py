import os
from Bio import SeqIO
import numpy as np
import pandas as pd
import argparse
import khmer


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


def local_variable():
	"""
	For test purposes
	:return:
	"""
	genome_list = os.path.abspath('input_genomes.txt')
	#genome_list = os.path.abspath('mac_studio_input_file.txt')
	kvalue = 12
	svalue = 8
	fraction_factor = 5
	rev_comp = False
	shift_high = int((kvalue - svalue) / 2)
	genome_list = check_files(genome_list)
	
	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Manual validation of syncmer distance",
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
	shift_high = int((kvalue - svalue) / 2)
	
	rev_comp = args.rev_comp
	rev_comp = rev_comp == 'True'
	if rev_comp:
		print("Using canonical kmers!")
		
	
	### process file in a brute force manner
	print("Getting distance distribution......")
	for shift_value in range(shift_high + 1):
		print("shift value is "+str(shift_value))
		df_to_merge = []
		
		for fasta in genome_list:
			temp_seq_list = generate_string_list_from_genome(fasta)
			# merge seqs for convenience, will introduce few extra kmers, but won't affect results as there are few contigs
			temp_seq = "".join(temp_seq_list)
			# record syncmer location
			temp_location = []
			
			# loop through seq
			for i in range(len(temp_seq) - kvalue + 1):
				current_kmer = temp_seq[i:i + kvalue]
				if is_kmer_open_syncmer_by_hash_value(current_kmer, s=svalue, shift=shift_value):
					temp_location.append(i)
			
			# check distance
			temp_location = np.array(temp_location)
			temp_distance = temp_location[1:] - temp_location[:-1]
			
			# change to freq dist.
			temp_name = os.path.basename(fasta)
			df_to_merge.append(pd.DataFrame(pd.Series(temp_distance).value_counts(), columns=[temp_name]))
		
		out_df_dist_syncmer = pd.concat(df_to_merge, axis=1)
		out_df_dist_syncmer = clean_dist_distribution_dataframe(out_df_dist_syncmer)
		out_df_dist_syncmer.to_csv("distance_distribution_k%ds%dt%d.csv" % (kvalue, svalue, shift_value), header=True,
		                           index=True)
	

	### added: check compression and sequence coverage change for different shift values
	out_df_compression = []
	out_df_coverage = []
	
	print("Getting compression and seq coverage......")
	for shift_value in range(shift_high+1):
		print("shift value is " + str(shift_value))
		list_filename = []
		list_compression = []
		list_coverage = []
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
					
			# compression ratio
			temp_compression = len(temp_location) * 1.0 / (len(temp_seq)-kvalue+1)
			list_compression.append(temp_compression)
					
			# get seq coverage
			temp_seq_cov = np.zeros(len(temp_seq), dtype=int)
			for loc in temp_location:
				temp_seq_cov[loc:loc+kvalue] = 1
			temp_coverage = np.sum(temp_seq_cov)*1.0 / len(temp_seq)
			list_coverage.append(temp_coverage)
		
		# merge to the dict
		out_df_compression.append(pd.DataFrame(list_compression, index=list_filename, columns=[str(shift_value)]))
		out_df_coverage.append(pd.DataFrame(list_coverage, index=list_filename, columns=[str(shift_value)]))
		
	# merge df
	merged_df_compression = pd.concat(out_df_compression, axis=1)
	merged_df_compression.to_csv("compression_change_by_shift_k%ds%dt.csv" %(kvalue, svalue), header=True, index=True)
	
	merged_df_coverage = pd.concat(out_df_coverage, axis=1)
	merged_df_coverage.to_csv("coverage_change_by_shift_k%ds%dt.csv" % (kvalue, svalue), header=True, index=True)
			
			
