# ResynPy

ResynPy is a time-efficient pipeline written in python that interrogates the genotyping data of a segregating population of plants, aiming to the detection of individuals with complementary genotyping profiles. This tool has been created with the aim to automatize, in an time-efficient way, the selection process of individuals as part of the Resynthesis method, as explained by [Eduardo _et al._ (2020)](https://www.frontiersin.org/articles/10.3389/fpls.2020.01205/full). The aim of Resynthesis is to obtain a variety that is a slightly modified variant of an elite commercial line but with new traits of agronomic importance. This method can greatly reduce the required time of variety production in plant species with high intergenerational periods, such as fruit trees. 

## Prerequisites

Apart from python, all the other packages are optional for the user to install. These packages are required for graphical representation of data generated by the pipeline.

- python (version 3.8.10 or higher)
- pandas (version 0.25.3 or higher)
- matplotlib (version 3.1.2 or higher)
- seaborn (0.11.2 or higher)

**_Note for Windows users_**: instead of using the Command Prompt window, a user can also use the "Terminal" tab in the lower left area of the [Rstudio](https://www.rstudio.com/) environment.

1. **Python installation**

   1. **In Linux**
   
      Follow this [tutorial](https://docs.python-guide.org/starting/install3/linux/) for installing Python 3 in Linux.
      
   2. **In Windows**
   
      Follow this [tutorial](https://phoenixnap.com/kb/how-to-install-python-3-windows) (Steps 1 to 5), for installing Python 3 in Windows 10. 
      
2. **Pandas, matplotlib and seaborn installation (Linux/Windows)**

   install [pandas](https://pandas.pydata.org/): `pip install pandas`
   
   install [matplotlib](https://matplotlib.org/): `pip install matplotlib`
   
   install [seaborn](https://seaborn.pydata.org/): `pip install seaborn`

## Running the pipeline

ResynPy performs pairwise individual comparison of the genotyping data of a segragating population for detecting complementary pairs of individuals.

**obselete**: ResynPy offers to the user two tasks of analysis: 1) Pairwise individual comparison of the genotyping data of a segragating population for detecting complementary pairs of individuals, and 2) Detection of heterozygous regions (Runs Of Heterozygosity; ROHet), using [detectRUNS](https://cran.r-project.org/web/packages/detectRUNS/vignettes/detectRUNS.vignette.html), in the parental line of the segregating population and annotation of the corresponding markers based on their presence inside or outside of ROHet.

To see all the avaiable options for running the pipeline, type `python (or python3) ResynPy.py -h`:

```
usage: resynthesis_master.py [-h] [--genos TAB FILE] [--markers TAB FILE] [--results_dir STR] [--not_phased] [--scores_file FILE] [--invariable FLOAT] [--hetero FLOAT] [--out_prefix STR]

Analysis of F2 genotyping data for detection of highly complementary individuals.

optional arguments:
  -h, --help          show this help message and exit
  --genos TAB FILE    Tab-delimited file (.tab) containing the genotyping data of the F2 population. Markers should be in columns and individuals in rows. Genotypes should be in the format "A,B,H,-", where "-" represents missing
                      data. This argument has to be used together with --markers
  --markers TAB FILE  A 2-column tab-delimited file with the markers used for the F2 genotyping, in the format of "chromosome<tab>marker name". Marker order per chromosome should match that of markers in the genotyping file.
  --results_dir STR   Name of the results directory [Default: ResynPy_results]
  --not_phased        Use this argument if your genotyping data are not phased. [Default: FALSE]
  --scores_file FILE  A tab delimited file containing user-defined scores for the different combinations of genotypes, assigned during the comparison of the individuals. Nucleotides ingenotypes should be separated with a "/".
                      e.g.: A/H<tab>0.75. [Default: scores_default.tab]
  --invariable FLOAT  Keep individual pairs that have a percentage of AA or BB combinations that is lower than the argument's value. [Default: 0.1]
  --hetero FLOAT      Keep individuals that have a heterozygosity lower than the argument's value. [Default: 0.5]
  --out_prefix STR    Prefix for the output files.[Default: resynpyOut]
```
- Example data provided in the repository:

   - `SDmarkers.tab`: markers file
   - `SD_F2data.tab`: genotyping file originating from [Eduardo et al., 2020](https://www.frontiersin.org/articles/10.3389/fpls.2020.01205/full)

### Pairwise comparison

During this process the genotypic profile of each individual is compared with the corresponding profile of the rest of the individuals in the F2 population. Each pair of genotypes in a specific marker has a score assigned to it. All scores assigned at each genotype pair are summed up and the total score is tagged to the respective pair of individuals. Scores for each combination of genotypes is provided in the file `scores_default.tab`. The user can also provide a tab-delimited file with different scores (see below). For the pairwise comparison analysis, the user should provide a file with the genotyping data (`--genos`) and a file with markers information (`--markers`) (see below). The pipeline considers that genotyping data are phased. If data are not phased then the user should declare it by using the argument `--not_phased`.

- Example command:

   ```
   python ResynPy.py --genos SD_F2data.tab --markers SDmarkers.tab [--results_dir ResynPy_results] [--filter_invariable 0.1] [--hetero 0.5] [--not_phased] [--scores_file scores_default.tab]
   ```

#### Input data

##### Mandatory arguments

`--genos`: A tab-delimited file with individuals in rows and markers in columns. Markers should ordered first by chromosome and then by their physical or genetic position. Individual column should have a title. Data should be in A, B, H format; i.e. "A" represents the homogyzous genotype for one allele, "B" the homozygous genotype for the other allele and "H" the heterozygous genotype. Missing data should be represented as "-".

   - example of genotyping data.

         indv	SNP1	SNP2  SNP3  SNP4  SNP5  SNP6  SNP7  SNP8  SNP9  SNP10 ...
         indv_1	A	-	A	-	B	A	A	H	B	B	...
         indv_2	A	-	B	-	B	H	B	H	H	H	...
         indv_3	A	-	H	-	B	A	B	A	H	B	...
         indv_4	B	-	H	-	A	H	H	A	H	B	...
         ...      ....
         indv_n	A	-	H	-	H	H	A	H	A	A	...


`--markers`: A 2-column tab-delimited file. Column 1: marker name used for genotyping the population. Column 2: chromosome the marker belongs to. Marker order should match that of the genotyping file.

##### Optional arguments

`--not_phased`: Should be used if genotyping data are not phased. The use of this argument will cause the exclusion of recombination information form the output file (see _Output data_ Section).

`--scores_file`: A 2-column tab-delimited file containing scores assigned to each possible genotype combination, when two individuals are compared. Column 1: genotype pair, Column 2: score. By default, the pipeline will use the score information present inside file `scores_default.tab`, provided in the repository's main directory. However if the user wants to apply a different scoring, a tab-delimited file can be provided by using this argument. 

   - contents of `scores_default.tab`:

         A/A	0
         B/B	0
         -/-	0.5
         H/H	0.5
         H/-	0.625
         A/-	0.75
         B/-	0.75
         A/H	0.75
         B/H	0.75
         A/B	1

`--invariable`: We also offer the possibility to filter out pairs of individuals that have a number of genotype identities (A/A or B/B) higher than the threshold float value set by this argument. The use of this argument can speed up the whole comparison process, since the pipeline does not proceed to score assignment if the pair is discarded. In theory, there shouldn't be any genotype identity between the two individuals, since the presence of genotype identity implies that the corresponding region of the genome will not contribute to obtaining the near-identical variety when the two individuals are crossed. Therefore, the value of this argument should be as low as possible. Default value: _0.01_

`--hetero`: A value defining the maximum level of heterozygosity for any indvidual in the pair. Low heterozygosity of an individual increases the probability to find a matching individual with a complementary genotype in the following generation. Default value: _0.5_

`--out_prefix`: Prefix to be used for all the output files. Default value: _resynpyOut_


#### Output data

All the output files are saved into the `--results_dir` directory (Default results directory: *ResynPy_results*).

- `selected_pairs.tab`

   File containing all the accepted pairs of individuals. This file can contain up to seven columns.

   column 1: selected pair of individuals

   column 2: total score of genotyping

   column 3: total number of invariable pairs (AA or BB)

   column 4: average pair heterozygosity
   
   _#### The following three columns are included in the output if data are phased._

   column 5: number of recombinations present in the first individual

   column 6: number of recombinations present in the second individual

   column 7: total number of recombinations in the pair

   Example of output file for _unphased_ data:

   ```
   #indiv_pair	score	invariable_pairs	average_heterozygosity
   SDF2-131|SDF2-180	52.875	0	0.29
   SDF2-180|SDF2-190	51.375	0	0.31
   SDF2-131|SDF2-146	50.625	2	0.23
   SDF2-156|SDF2-222	50.375	0	0.35
   SDF2-188|SDF2-300	50.375	0	0.36
   ```

   Example of output file for _phased_ data:

   ```
   #indiv_pair	score	invariable_pairs	average_heterozygosity	recomb_no1	recomb_no2	sum_recomb
   SDF2-131|SDF2-180	52.875	0	0.29	1	0	1
   SDF2-180|SDF2-190	51.375	0	0.31	0	0	0
   SDF2-131|SDF2-146	50.625	2	0.23	1	5	6
   SDF2-156|SDF2-222	50.375	0	0.35	0	0	0
   SDF2-188|SDF2-300	50.375	0	0.36	0	0	0
   ```

- `discarded_individuals_and_pairs.tab`

   File containing all individuals and/or pairs that have not met the filters set by the user, supplemented with the parameter that did not meet the criteria and its value.
   
   Example:

   ```
   PR20CG15-1840-3988	heterozygosity: 0.43
   PR20CG15-1840-3995	heterozygosity: 0.46
   PR20CG15-1840-002|PR20CG15-1840-003	invariable genotype ratio: 0.27
   PR20CG15-1840-002|PR20CG15-1840-004	invariable genotype ratio: 0.23
   PR20CG15-1840-002|PR20CG15-1840-204	pair score: 19.75
   ```

- `<prefix>.png`, `<prefix>.pdf`

   A figure showing the haplotypes of the individuals in the first 10 selected pairs, in a color-coded mode. `<prefix>` corresponds to argument `--out_prefix`.
   

- `PairComparison.log`

   Log file with the arguments used by the user, the different steps implemented by the pipeline and counts of accepted and discarded individuals and pairs. When the user runs again the pipeline the file is appended with the new information. This can help the user to evaluate different parameters of the pipeline.
