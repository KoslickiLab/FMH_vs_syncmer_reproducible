## FracMinHash vs Open syncmers in genomic comparisons
Reproducible scripts for the manuscript `Connecting Syncmers to FracMinHash: similarities and advantages`. The associated preprint will be added soon.
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

### Download random genomes for comparison
```
# Here are 20 random Brucella representative genomes from GTDB r214: https://data.gtdb.ecogenomic.org/releases/release214/214.0/
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

datasets download genome accession $(cat rand_20_brucella_genome.txt | tr "\n" "\t")
mkdir -p brucella_genomes
unzip ncbi_dataset.zip
find . -name "*.fna" | xargs -I{} mv {} ./brucella_genomes
rm -r ncbi_dataset*
rm README.md
readlink -f brucella_genomes/*.fna > file_paths.txt
```

</br>

### Generate k-mer sketches and compare the performance of FMH and open syncmers
The default parameters are `k=20, s=11, c=10` and we will use `shift value` from 1 to 5.  
This step may run ~2h on a single-thread laptop, reducing input genomes (e.g. keep only first 3 records) can make it much faster for testing purpose.
```
# assume we are in the demo folder now
python ../src/reproducible.py -f file_paths.txt -k 20 -s 11 -a 10
```

</br>

### Generate figures based on comparison results  
We will obtain all results in the previous step, this is the last step and is for visulization.
```
# assume we are in the demo folder now
python ../src/make_figs_for_manuscript.py
```
</br>
</br>

---

### Parameter explanation and output files
After running the steps above, you will have a list of output files in the `demo` folder.
#### Parameters:
| Parameter | Note                                                         |
| --------- | ------------------------------------------------------------ |
| -f        | Input file with paths of all input genomes, one per line     |
| -k        | k-mer length, default 20                                     |
| -s        | submer length, default 11                                     |
| -a        | fraction factor for FracMinHash, default 10. Please make sure k-s+1=a |

</br>

#### Output files:
List of output files and contents from the main script.
| Filename                             | Content                                                      |
| ------------------------------------ | ------------------------------------------------------------ |
| compare_conservation_MOS_s11_mut0.05.csv  | Multi-resolution estimation by MOS                             |
| compare_conservation_regular_OS_s11_mut0.05.csv               | Multi-resolution estimation by regular open syncmers |
| compare_jici.csv            | Estimation of Jaccard and containment index by OS and FMH |
| conservation_change_by_shift_k20s11_mut0.01.csv | K-mer conservation status by FMH and OS |
| coverage_change_by_shift_k20s11.csv                | Genome coverage status by FMH and OS     |
| coverage_change_by_shift_k20s11.csv                | Genome coverage status by FMH and OS     |
| distance_distribution_FMH_k20c10.csv | Frequency distribution of distances between adjacent k-mers in the genome for FMH and OS |

</br>

#### Output figures:
List of output figures.

| Filename                              | Content                              |
| ------------------------------------- | ------------------------------------ |
| Fig2_compression.png | Figure2: Open syncmer vs FracMinHash: compression ratio and k-mer conservation are consistent. |
| Fig3_ji_est.png | Figure 3: The open syncmer sketch fits the MinHash algorithm. |
| Fig4_dist_coverage_compare.png    | Figure 4: Genome coverage and distance distribution |
| Fig5_mos.png | Figure 5: Multi-resolution open syncmers are comparable to regular open syncmers |






