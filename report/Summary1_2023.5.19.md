#### Terminology

| Item      | Explanation                                                  |
| --------- | ------------------------------------------------------------ |
| Syncmer   | We only discuss open syncmer here. Open syncmer is a k-mer of which the minimal s-mer (defined by hash value) appear at a specific position t. A syncemer is defined by 3 parameters: k, s and t. |
| k         | k-mer length                                                 |
| s         | submer length, s<k.                                          |
| t / shift | The location of minimal s-mer for a k-mer to be considered as an open syncmer. Usually t=0 (at the beginning), but we can get a different syncmer set by using different shift values. In ERS-mer, we use t = floor((k-s+1)/2) |
| Distance  | Distance refers to the distance between 2 consecutive k-mers in a FMH k-mer set or an open syncmer set |





#### Observed in experiment:

1. Open syncmer set and FracMinhash sketches are similar in the following aspects:
   1. same compression ratio (sketch size / number of total k-mers), no matter what shift value to use.
   2. same weight profile for all selected k-mers
   3. same k-mer-based JI/CI similarity estimates
   4. same seed conservations at gradiant mutation rates
2. But they are different in:
   1. **Open syncmer has different distance distribution. And the distance distribution depends on shift value**. 
   2. Specifically, using a large shift value makes a much metter distance distribution: **it erases out short distances** (i.e. no adjacent k-mers) and gives a larger sequence coverage.





#### Theoretical analysis based on observations:

- 2.1 already indicates that syncmer algorithm is NOT identical to FMH algorithm, because overlapped syncmers are NOT independent. **This contracdicts my previous thought, and I don't know where I got wrong for the proof of independence**.
- Combined 1.1 with 2.1 and 2.2. It's hard to explain why we still have the same compression ratio for syncmers with large shift values given that they have distinct distance distribution.





#### Benchmark results from 20 random Brucella genomes

1. Distance distribution (please see attached files in the email)

   1. for FMH k-mer set, the distance follows binomial distribution, and its decreasing
   2. **for syncmer with shift=0, there is always a peak at distance = k-s+1 (no idea why, but may worth digging)**
   3. for syncmer with shift=t, there are very few distances in [1, t-1], and there is no significant decreasing in several locations. This makes a more evenly distributed seed

2. Compression NOT affected by shift: k=12, s=5, shift=0,1,2,3

   | Shift                                                        | 0          | 1          | 2          | 3          |
   | ------------------------------------------------------------ | ---------- | ---------- | ---------- | ---------- |
   | GCA_002291225.1_ASM229122v1_genomic.fna                      | 0.12851773 | 0.12833959 | 0.12813234 | 0.12892379 |
   | GCA_009733355.1_Brucella_abortus_C-577_genomic.fna           | 0.12570208 | 0.12548439 | 0.12525999 | 0.12620636 |
   | GCA_014884585.1_ASM1488458v1_genomic.fna                     | 0.12570583 | 0.12541376 | 0.12512199 | 0.12616977 |
   | GCF_000292205.1_V02_for_genome_sequence_of_BCB013_genomic.fna | 0.12585807 | 0.12557885 | 0.12531766 | 0.12621793 |
   | GCF_000366825.1_Bruc_cani_UK10_02_V1_genomic.fna             | 0.12572277 | 0.12540987 | 0.12513113 | 0.12617807 |
   | GCF_000413595.1_Bruc_abor_biovar_1_85-1058_V1_genomic.fna    | 0.12578939 | 0.12547445 | 0.12524807 | 0.12621082 |
   | GCF_000413755.1_Bruc_abor_biovar_1_01-0648_V1_genomic.fna    | 0.12578082 | 0.12546776 | 0.1252481  | 0.12621404 |
   | GCF_000480275.1_Bruc_cani_96-7258_V1_genomic.fna             | 0.12571599 | 0.12540569 | 0.12512557 | 0.12617148 |
   | GCF_001866335.1_ASM186633v1_genomic.fna                      | 0.13523206 | 0.13498998 | 0.134758   | 0.13558983 |
   | GCF_002150355.1_ASM215035v1_genomic.fna                      | 0.12584941 | 0.12553654 | 0.1252902  | 0.12616926 |
   | GCF_002191555.1_ASM219155v1_genomic.fna                      | 0.12586034 | 0.12555459 | 0.1252772  | 0.12631642 |
   | GCF_002191835.1_ASM219183v1_genomic.fna                      | 0.12573816 | 0.12542391 | 0.1251486  | 0.12618645 |
   | GCF_003994435.1_ASM399443v1_genomic.fna                      | 0.13470175 | 0.13439164 | 0.13416163 | 0.13511383 |
   | GCF_006507335.1_ASM650733v1_genomic.fna                      | 0.12582431 | 0.12558147 | 0.12537849 | 0.12617641 |
   | GCF_009664925.1_ASM966492v1_genomic.fna                      | 0.12575745 | 0.12544799 | 0.12534263 | 0.12621429 |
   | GCF_009733375.1_Brucella_abortus_C-551_genomic.fna           | 0.12574437 | 0.12552043 | 0.12531587 | 0.12623316 |
   | GCF_023651255.1_ASM2365125v1_genomic.fna                     | 0.13285523 | 0.13260354 | 0.13237234 | 0.13312258 |
   | GCF_023651275.1_ASM2365127v1_genomic.fna                     | 0.13363801 | 0.13339516 | 0.13316875 | 0.13391458 |
   | GCF_900236395.1_Bru_genome41_genomic.fna                     | 0.12574411 | 0.12555626 | 0.12542395 | 0.12616958 |
   | GCF_900236485.1_Bru_genome24_genomic.fna                     | 0.12575287 | 0.12556531 | 0.1254512  | 0.1261823  |

3. Sequence coverage increases when shift increase: k=12, s=5, shift=0,1,2,3

4. | Shift                                                        | 0          | 1          | 2          | 3          |
   | ------------------------------------------------------------ | ---------- | ---------- | ---------- | ---------- |
   | GCA_002291225.1_ASM229122v1_genomic.fna                      | 0.89512081 | 0.93079543 | 0.95585899 | 0.9694484  |
   | GCA_009733355.1_Brucella_abortus_C-577_genomic.fna           | 0.89453661 | 0.93022036 | 0.95542368 | 0.96932338 |
   | GCA_014884585.1_ASM1488458v1_genomic.fna                     | 0.8944967  | 0.9301876  | 0.95531065 | 0.96934032 |
   | GCF_000292205.1_V02_for_genome_sequence_of_BCB013_genomic.fna | 0.89490204 | 0.93039027 | 0.95566041 | 0.96957973 |
   | GCF_000366825.1_Bruc_cani_UK10_02_V1_genomic.fna             | 0.89461313 | 0.93030249 | 0.95541763 | 0.9694159  |
   | GCF_000413595.1_Bruc_abor_biovar_1_85-1058_V1_genomic.fna    | 0.89422975 | 0.92997558 | 0.95532379 | 0.96944796 |
   | GCF_000413755.1_Bruc_abor_biovar_1_01-0648_V1_genomic.fna    | 0.89422503 | 0.92998544 | 0.95532558 | 0.96946083 |
   | GCF_000480275.1_Bruc_cani_96-7258_V1_genomic.fna             | 0.8945965  | 0.93028781 | 0.95540322 | 0.96940965 |
   | GCF_001866335.1_ASM186633v1_genomic.fna                      | 0.89590935 | 0.93140317 | 0.95622675 | 0.96977351 |
   | GCF_002150355.1_ASM215035v1_genomic.fna                      | 0.8946166  | 0.93070112 | 0.95579498 | 0.96969345 |
   | GCF_002191555.1_ASM219155v1_genomic.fna                      | 0.89454681 | 0.9302771  | 0.9554166  | 0.96941539 |
   | GCF_002191835.1_ASM219183v1_genomic.fna                      | 0.89457057 | 0.93027752 | 0.95541402 | 0.96941074 |
   | GCF_003994435.1_ASM399443v1_genomic.fna                      | 0.89528987 | 0.93068351 | 0.95577867 | 0.96979338 |
   | GCF_006507335.1_ASM650733v1_genomic.fna                      | 0.89451799 | 0.93030169 | 0.95535467 | 0.96958927 |
   | GCF_009664925.1_ASM966492v1_genomic.fna                      | 0.89445738 | 0.93012991 | 0.95567963 | 0.96955351 |
   | GCF_009733375.1_Brucella_abortus_C-551_genomic.fna           | 0.89428949 | 0.93009632 | 0.95531813 | 0.96951505 |
   | GCF_023651255.1_ASM2365125v1_genomic.fna                     | 0.89545488 | 0.93100043 | 0.95595128 | 0.96951106 |
   | GCF_023651275.1_ASM2365127v1_genomic.fna                     | 0.89554314 | 0.93110102 | 0.95601794 | 0.96956339 |
   | GCF_900236395.1_Bru_genome41_genomic.fna                     | 0.89446407 | 0.93038369 | 0.95559233 | 0.96953192 |
   | GCF_900236485.1_Bru_genome24_genomic.fna                     | 0.89445705 | 0.93034353 | 0.95564879 | 0.96957044 |

   Both share similar weight profile, pairwise CI/JI, and conservations (previously showed): k=20, s=11

   3.1 Partial CI/JI estimates

   | gt_ji | fmh_ji | sync_ji | gt_ci | fmh_ci | sync_ci |
   | ----- | ------ | ------- | ----- | ------ | ------- |
   | 0.263 | 0.263  | 0.263   | 0.421 | 0.421  | 0.421   |
   | 0.105 | 0.105  | 0.105   | 0.193 | 0.193  | 0.193   |
   | 0.316 | 0.316  | 0.316   | 0.485 | 0.485  | 0.484   |
   | 0.105 | 0.106  | 0.105   | 0.193 | 0.194  | 0.193   |
   | 0.218 | 0.218  | 0.219   | 0.361 | 0.361  | 0.362   |
   | 0.218 | 0.218  | 0.219   | 0.361 | 0.36   | 0.362   |
   | 0.105 | 0.106  | 0.105   | 0.193 | 0.193  | 0.193   |
   | 0.366 | 0.367  | 0.366   | 0.538 | 0.539  | 0.538   |
   | 0.249 | 0.249  | 0.249   | 0.403 | 0.404  | 0.403   |
   | 0.099 | 0.099  | 0.099   | 0.182 | 0.183  | 0.182   |
   | 0.099 | 0.099  | 0.099   | 0.182 | 0.183  | 0.182   |
   | 0.218 | 0.217  | 0.218   | 0.36  | 0.359  | 0.36    |
   | 0.307 | 0.306  | 0.307   | 0.475 | 0.473  | 0.474   |
   | 0.176 | 0.175  | 0.175   | 0.304 | 0.304  | 0.303   |
   | 0.321 | 0.321  | 0.322   | 0.489 | 0.489  | 0.49    |

​		

​		3.2 Conservation of 1 file at different mutation rate (repeat: 10 times)

| mut_rate | kmer_conservation | syncmer_conservation |
| -------- | ----------------- | -------------------- |
| 0.01     | 0.83451016        | 0.8354163            |
| 0.03     | 0.57639097        | 0.5772787            |
| 0.05     | 0.39070772        | 0.3896685            |
| 0.07     | 0.26109717        | 0.26093048           |
| 0.09     | 0.17508241        | 0.17494391           |
| 0.11     | 0.11376548        | 0.11352225           |
| 0.13     | 0.07444623        | 0.07456851           |
| 0.15     | 0.0486621         | 0.04784105           |
| 0.17     | 0.0302948         | 0.03005395           |
| 0.19     | 0.01957534        | 0.01931304           |



​		3.3 Weight profile of 1 file

| Weight | fmh    | sync   |
| ------ | ------ | ------ |
| 1      | 327120 | 327384 |
| 2      | 501    | 474    |
| 3      | 91     | 115    |
| 4      | 41     | 43     |
| 5      | 18     | 18     |
| 6      | 9      | 13     |
| 7      | 5      | 10     |
| 8      | 5      | 3      |
| 9      | 1      | 3      |
| 10     |        | 1      |
| 11     | 2      | 2      |
| 12     | 4      | 3      |
| 13     | 1      | 2      |
| 14     |        | 1      |
| 15     | 3      | 3      |
| 17     |        | 1      |
| 19     | 1      | 1      |
| 20     | 1      | 1      |
| 21     |        | 1      |
| 25     | 1      |        |



