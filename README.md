## FracMinHash vs Open syncmer in random seeding

### Environment setup

Please use conda environment with Python 3.7 to avoid dependency conflicts. Run the following code:

```
conda create -n fmh_vs_syncmer python=3.7
conda activate fmh_vs_syncmer
conda install -c conda-forge -c bioconda -c anaconda --file requirements.txt
```



### Download random genomes for comparison

```
# Here are 20 random Brucella rep genomes from GTDB r214 here: https://data.gtdb.ecogenomic.org/releases/release214/214.0/
# grep Brucella bac120_taxonomy_r214.tsv | shuf --random-source=<(yes 42) | head -20 | cut -f 1 | sed 's/.*_GC/GC/g'

mkdir -p demo
cd demo
echo "GCF_000413755.1
GCF_001866335.1
GCF_002191835.1
GCF_000366825.1
GCF_009733375.1
GCF_900236485.1
GCF_002150355.1
GCA_009733355.1
GCF_023651275.1
GCF_000292205.1
GCF_009664925.1
GCF_023651255.1
GCF_000413595.1
GCF_900236395.1
GCF_003994435.1
GCF_006507335.1
GCA_002291225.1
GCA_014884585.1
GCF_002191555.1
GCF_000480275.1" > rand_20_brucella_genome.txt

datasets download genome accession $(cat rand_20_brucella_genome.txt | tr "\n" "\t" )
mkdir -p brucella_genomes
unzip ncbi_dataset.zip
find . -name "*.fna" | xargs -I{} mv {} ./brucella_genomes
rm -r ncbi_dataset*
rm README.md
readlink -f brucella_genomes/*.fna > file_paths.txt
```



### Run the main script to generate FMH vs syncmer comparison

You can also use arbitrary fasta files as input. Store file paths in one file and pass it to the script by "-f" option.

```
# assume we are in the demo folder now
mkdir -p result_random20_Brucella_k20s16ratio5
cd result_random20_Brucella_k20s16ratio5
python ../../script/compare_fmh_syncmer_sketches.py -f ../file_paths.txt -k 20 -s 16 -a 5

cd ..
mkdir -p result_random20_Brucella_k20s11ratio10
cd result_random20_Brucella_k20s11ratio10
python ../../script/compare_fmh_syncmer_sketches.py -f ../file_paths.txt -k 20 -s 11 -a 10
```

More parameters:

| Parameter | Note                                                         |
| --------- | ------------------------------------------------------------ |
| -f        | Input file with paths of all input genomes, one per line     |
| -k        | k-mer length, default 12                                     |
| -s        | submer length, default 8                                     |
| -a        | fraction factor for FracMinHash, default 5. Please make sure k-s+1=a |
| -r        | use canonical k-mers, default False.                         |



### Check how submer shift affect sketch performance

```
# assume we are in the demo folder now
mkdir -p test_submer_shift
cd test_submer_shift
python ../../re-test_syncmer_distance_distribution.py -f ../file_paths.txt -k 12 -s 8
```

Only requires input k and s value besides genome files. This script will compare distance distribution, compression, and sequence coverage of given k and s regarding gradiant shift values from 0 to (k-s)/2





### Results interpretation

List of output files and contents from the main script

| Filename                             | Content                                                      |
| ------------------------------------ | ------------------------------------------------------------ |
| out_1_compare_compression_ratio.csv  | Compare sketch compression ratio                             |
| out_2_compare_jici.csv               | Compare pairwise JI/CI estimated by FMH, syncmer and ground truth |
| out_3_dist_distribution              | Frequency ratio of distance between 2 consecutive seeds      |
| out_4_compare_conservation_level.csv | Compare mean seed conservation ratio under gradient mutation rates. This step is time consuming, so only applied to the 1st input file. |
| out_5_weight_dist                    | Compare weight of 2 sketches                                 |



List of output files and contents from the added script for affection due to submer shift

| Filename                              | Content                              |
| ------------------------------------- | ------------------------------------ |
| compression_change_by_shift_k_s_t.csv | Compression ratio by shift value     |
| coverage_change_by_shift_k_s_t.csv    | Sequence coverate by shift value     |
| distance_distribution_k_s_t_.csv      | Distance distribution by shift value |





