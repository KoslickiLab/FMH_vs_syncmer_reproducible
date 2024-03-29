## FracMinHash vs Open syncmers in genomic comparisons

Reproducible scripts for the manuscript `Connecting Syncmers to FracMinHash: similarities and advantages`. The associated preprint will be added soon.
</br>

<!-- TOC start -->

### Contents:

- [Step1](#environment-setup): environment setup
- [Step2](#download-random-genomes-for-comparison): download random genomes for comparison
- [Step3](#generate-k-mer-sketches-and-compare-the-performance-of-fmh-and-open-syncmers): generate k-mer sketches and compare the performance of FMH and open syncmers
- [Step4](#generate-figures-based-on-comparison-results): generate figures based on comparison results
- [Parameter explanation and output files](#parameter-explanation-and-output-files)
  * [Parameters](#parameters)
  * [Output files](#output-files)
  * [Output figures](#output-figures)

<!-- TOC end -->
<!-- TOC --><a name="environment-setup"></a>

</br>

### Environment setup

```
git clone https://github.com/KoslickiLab/FMH_vs_syncmer_reproducible.git
cd FMH_vs_syncmer_reproducible
conda create -n fmh_vs_syncmer
conda activate fmh_vs_syncmer
conda install -c conda-forge -c bioconda -c anaconda --file ./src/requirements.txt
```

</br>

<!-- TOC --><a name="download-random-genomes-for-comparison"></a>

### Download random genomes for comparison

```
mkdir -p demo
cd demo

# download GTDB taxonomy file
wget https://data.gtdb.ecogenomic.org/releases/release214/214.0/bac120_taxonomy_r214.tsv

# get 20 brucella genomes
grep Brucella bac120_taxonomy_r214.tsv | shuf --random-source=<(yes 42) | head -20 | cut -f 1 | sed 's/.*_GC/GC/g' > rand_20_brucella_genome.txt

datasets download genome accession $(cat rand_20_brucella_genome.txt | tr "\n" "\t" )
mkdir -p brucella_genomes
unzip ncbi_dataset.zip
find . -name "*.fna" | xargs -I{} mv {} ./brucella_genomes
rm -r ncbi_dataset*
rm README.md
readlink -f brucella_genomes/*.fna > filepath_20_brucella.txt


# get more genomes from randomly selected families
grep -oP '(?<=;f__)[^;]*(?=;)' bac120_taxonomy_r214.tsv | sort | uniq | sed '/[0-9]/d' > gtdb_family_list.txt
> rand_genome_list.txt
for random_family in $(cat gtdb_family_list.txt | shuf --random-source=<(yes 42) | head -50 ); do
 grep ${random_family} bac120_taxonomy_r214.tsv | shuf --random-source=<(yes 42) | head -10 | cut -f 1 | sed 's/.*_GC/GC/g' >> rand_genome_list.txt 
done

mkdir -p more_genomes_from_50_families
cd more_genomes_from_50_families
datasets download genome accession $(cat ../rand_genome_list.txt | tr "\n" "\t")
unzip ncbi_dataset.zip
find . -name "*.fna" | xargs -I{} mv {} .
rm -r ncbi_dataset*
rm README.md
readlink -f *.fna > ../filepath_more_genomes_from_rand_50_families.txt
cd ..


# merge file paths together
cat filepath_20_brucella.txt filepath_more_genomes_from_rand_50_families.txt > file_paths.txt
```

</br>

<!-- TOC --><a name="generate-k-mer-sketches-and-compare-the-performance-of-fmh-and-open-syncmers"></a>

### Generate k-mer sketches and compare the performance of FMH and open syncmers

The default parameters are `k=20, s=11, c=10` and we will use `shift value` from 1 to 5.  
This step may run ~2h on a single-thread laptop, reducing input genomes (e.g. keep only the first 3 records) can make it much faster for testing purpose.

```
# assume we are in the demo folder now
# if want a faster test run, try:
# head -2 file_paths.txt > test_file.txt  
python ../src/reproducible.py -f file_paths.txt -k 20 -s 11 -a 10
```

</br>

<!-- TOC --><a name="generate-figures-based-on-comparison-results"></a>

### Generate figures based on comparison results  

We will obtain all results in the previous step, this is the last step and is for visulization.

```
# assume we are in the demo folder now
python ../src/make_figs_for_manuscript.py
```

</br>
</br>

---

<!-- TOC --><a name="parameter-explanation-and-output-files"></a>

### Parameter explanation and output files

After running the steps above, you will have a list of output files in the `demo` folder.
<!-- TOC --><a name="parameters"></a>

#### Parameters:

| Parameter | Note                                                         |
| --------- | ------------------------------------------------------------ |
| -f        | Input file with paths of all input genomes, one per line     |
| -k        | k-mer length, default 20                                     |
| -s        | submer length, default 11                                    |
| -a        | fraction factor for FracMinHash, default 10. Please make sure k-s+1=a |

</br>

<!-- TOC --><a name="output-files"></a>

#### Output files:

List of output files and contents from the main script.

| Filename                                        | Content                                                      |
| ----------------------------------------------- | ------------------------------------------------------------ |
| compare_conservation_MOS_s11_mut0.05.csv        | Multi-resolution estimation by MOS                           |
| compare_conservation_regular_OS_s11_mut0.05.csv | Multi-resolution estimation by regular open syncmers         |
| compare_jici.csv                                | Estimation of Jaccard and containment index by OS and FMH    |
| conservation_change_by_shift_k20s11_mut0.01.csv | K-mer conservation status by FMH and OS                      |
| coverage_change_by_shift_k20s11.csv             | Genome coverage status by FMH and OS                         |
| coverage_change_by_shift_k20s11.csv             | Genome coverage status by FMH and OS                         |
| distance_distribution_FMH_k20c10.csv            | Frequency distribution of distances between adjacent k-mers in the genome for FMH and OS |

</br>

<!-- TOC --><a name="output-figures"></a>

#### Output figures:

List of output figures.

| Filename                       | Content                                                      |
| ------------------------------ | ------------------------------------------------------------ |
| Fig2_compression.png           | Figure2: Open syncmer vs FracMinHash: compression ratio and k-mer conservation are consistent. |
| Fig3_ji_est.png                | Figure 3: The open syncmer sketch fits the MinHash algorithm. |
| Fig4_dist_coverage_compare.png | Figure 4: Genome coverage and distance distribution          |
| Fig5_mos.png                   | Figure 5: Multi-resolution open syncmers are comparable to regular open syncmers |


