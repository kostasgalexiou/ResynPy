# ResynPy

ResynPy is a time-efficient pipeline written in python that interrogates the genotyping data of a segregating population of plants, aiming to the detection of individuals with complementary genotyping profiles. This tool has been created with the aim to automatize, in an efficient way, the selection process of individuals as part of the Resynthesis method, as explained by [Eduardo _et al._ (2020)](https://www.frontiersin.org/articles/10.3389/fpls.2020.01205/full). Briefly, the aim of Resynthesis is to obtain a variety that is a slightly modified variant of an elite commercial line but with new traits of agronomic importance. This method can greatly reduce the required time of variety production in plant species with high intergenerational periods, such as fruit trees. 

## Prerequisites

Apart from python, all the other packages are optional for the user to install. These packages are required for graphical representation of data generated by the pipeline.

- python (version 3.8.10 or higher)
- bcftools (version 1.10.2 or higher) # **need more info**
- pandas (version 0.25.3 or higher)
- matplotlib (version 3.1.2 or higher)
- seaborn (0.11.2)

**_Note for Windows users_**: instead of using the Command Prompt window, a user can also use the "Terminal" tab in the lower left area of the [Rstudio](https://www.rstudio.com/) environment.

1. **Python installation**

   1. **In Linux**
   
      Follow this [tutorial](https://docs.python-guide.org/starting/install3/linux/) for installing Python 3 in Linux.
      
   2. **In Windows**
   
      Follow this [tutorial](https://phoenixnap.com/kb/how-to-install-python-3-windows) for installing Python 3 in Windows 10. 
      
2. **Pandas, matplotlib and seaborn installation**

   Type `pip install pandas` in the Terminal (Linux and Windows)
   
   **to be continued**...probably all these can be installed from the terminal window in Rstudio, for Windows users
   
## Running the pipeline

XXXXX offers to the user two tasks of analysis: 1) Pairwise individual comparison in the genotyping data of a segragating population for detecting complementary pairs of individuals, and 2) Detection of homozygous regions (Runs Of Homozygosity; ROH) in the parental line of the segregating population and annotation of the corresponding markers based on their presence inside or outside of ROH.

To see all the avaiable options for running the pipeline, type `python3 XXXXXX.py -h`:

```
usage: resynthesis_master_noImpute.py [-h] [--vcf VCF FILE] [--results_dir STR] [--genos TAB FILE] [--markers TAB FILE] [--not_phased] [--scores_file FILE] [--filter_invariable FLOAT] [--hetero FLOAT]

Detection of ROH regions and analysis of F2 genotyping data for detection of highly complementary individuals.

optional arguments:
  -h, --help            show this help message and exit
  --genos TAB FILE      Tab-delimited file (.tab) containing the genotyping data of the F2 population. Markers should be in columns and individuals in rows. Genotypes should be in the format "A,B,H,-", where
                        "-" represents missing data. Incomapatible with --vcf. This arguments is used together with --markers
  --markers TAB FILE    A 2-column tab-delimited file with the markers used for the F2 genotyping, in the format of "chromosome<tab>marker name". Incombatible with "--vcf"
  --results_dir STR     Name of the results directory [Default: ROH_results]
  --not_phased          Use this argument if your genotyping data are not phased. [Default: FALSE]
  --scores_file FILE    A tab delimited file containing user-defined scores for the different combinations of genotypes, observed during the comparison of the individuals. Nucleotides ingenotypes should be
                        separated with a "/". e.g.: A/H<tab>0.75. [Default: scores_default.tab]
  --filter_invariable FLOAT
                        Keep individual pairs that have a percentage of AA or BB combinations that is lower than the argument value. [Default: 0.1]
  --hetero FLOAT        Keep individuals that have a heterozygosity lower than the argument value. [Default: 0.5]

  --vcf VCF FILE        Name of the VCF file for the line, for which the ROH regions will be detected.

```

### Pairwise comparison

During this process the genotypic profile of each individual is compared with the corresponding profile of the rest of the individuals in the segregating population and the process is completed when all the possible one-way pairwise comparisons are made. For this type of analysis the user must provide a file with the genotyping data (`--genos`), a file with markers information (`--markers`) (see below). The pipeline considers that genotyping data are phased. If data are not phased then the user should declare it by using the argument `--not_phased`.

#### Input data

##### Mandatory arguments

`--genos`: a tab-delimited file with individuals in rows and markers in columns. Individual column should have a title. Data should be in A, B, H format; i.e. "A" represents the homogyzous genotype for one allele, "B" the homozygous genotype for the other allele and "H" the heterozygous genotype. Missing data should be represented as "-".

   - example of genotyping data.

         indv	SNP1	SNP2  SNP3  SNP4  SNP5  SNP6  SNP7  SNP8  SNP9  SNP10 ...
         indv_1	A	-	A	-	B	A	A	H	B	B	...
         indv_2	A	-	B	-	B	H	B	H	H	H	...
         indv_3	A	-	H	-	B	A	B	A	H	B	...
         indv_4	B	-	H	-	A	H	H	A	H	B	...
         ...      ....
         indv_n	A	-	H	-	H	H	A	H	A	A	...


`--markers`: a 2-column tab-delimited file. Column 1: marker name used for genotyping the population. Column 2: chromosome the marker belongs to. File should be sorted by chromosome.

##### Optional arguments

`--not_phased`: this argument should be used we used if genotyping data are not phased. The use of this argument will cause the omitance of recombination information in the output file (see _Output data_ Section).

`--scores_file`: a 2-column tab-delimited file containing scores assigned to each possible genotype combination, when two individuals are compared. Column 1: genotype pair, Column 2: score. If no file is provided, the pipeline will use the score information present inside file `scores_default.tab`, provided in the repository's main directory.

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

`--filter_invariable`: we also offer the possibility to filter out pairs of individuals that have a number of genotype identities (A/A or B/B) higher than the threshold float value set by this argument. The use of this argument can speed up the whole comparison process, since the pipeline does not proceed to score assignment if the pair is to be discarded. In theory, there shouldn't be any genotype identity between the two individuals, since the presence of genotype identity implies that the corresponding region of the genome will not contribute to obtaining the near-identical variety when the two individuals are crossed. Therefore, the value of this argument should be as low as possible.

`--hetero`: a value defining the maximum level of heterozygosity for any indvidual in the pair. Low heterozygosity of an individual increases the probability to find a matching individual with a complementary genotype in the following generation, within a small pool....XXXXXX 

#### Output data

- `selected_pairs.tab`

   File containing all the accepted pairs of individuals. This file can contain up to seven columns.

   column 1: selected pair of individuals

   column 2: total score of genotyping

   column 3: number of AA genotype pairs

   column 4: number of BB genotype pairs

   column 5: number of recombinations present in the first individual

   column 6: number of recombinations present in the second individual

   column 7: total number of recombinations in the pair

   The last four columns are included in the output if data are phased.

   Example of output file for _unphased_ data:

   ```
   #indiv_pair	score	AA pairs	BB pairs
   PR20CG15-1840-1675|PR20CG15-1840-3523	22.5	0	0
   PR20CG15-1840-059|PR20CG15-1840-1911	22.25	0	0
   PR20CG15-1840-2382|PR20CG15-1840-3853	22.25	0	0
   PR20CG15-1840-1717|PR20CG15-1840-3848	22.0	0	0
   PR20CG15-1840-1872|PR20CG15-1840-3288	22.0	0	0
   PR20CG15-1840-1872|PR20CG15-1840-3433	22.0	0	0
   PR20CG15-1840-1455|PR20CG15-1840-3588	21.875	0	0
   ```

   Example of output file for _phased_ data:

   ```
   #indiv_pair	score	AA pairs	BB pairs	recomb_no1	recomb_no2	sum_recomb
   PR20CG15-1840-1675|PR20CG15-1840-3523	22.5	0	0	0	0	0
   PR20CG15-1840-059|PR20CG15-1840-1911	22.25	0	0	0	0	0
   PR20CG15-1840-2382|PR20CG15-1840-3853	22.25	0	0	0	0	0
   PR20CG15-1840-1717|PR20CG15-1840-3848	22.0	0	0	0	0	0
   PR20CG15-1840-1872|PR20CG15-1840-3288	22.0	0	0	0	0	0
   PR20CG15-1840-1872|PR20CG15-1840-3433	22.0	0	0	0	0	0
   PR20CG15-1840-1455|PR20CG15-1840-3588	21.875	0	0	0	0	0
   ```

- `discarded_individuals_and_pairs.tab`

   File containing all individuals and/or pairs that have not met the filters set by the user, supplemented with the parameter that did not meet the criteria and its value.

   ```
   PR20CG15-1840-3988	heterozygosity: 0.43
   PR20CG15-1840-3995	heterozygosity: 0.46
   PR20CG15-1840-002|PR20CG15-1840-003	invariable genotype ratio: 0.27
   PR20CG15-1840-002|PR20CG15-1840-004	invariable genotype ratio: 0.23
   PR20CG15-1840-002|PR20CG15-1840-204	pair score: 19.75
   ```

- `PairComparison.log`

   Log file with the arguments used by the user and the different steps run by the pipeline. When the user runs again the pipeline the file is appended with the new information, so that the user can have a calendar of the different runs that can serve for evaluating the different parameters of the pipeline.


### ROH detection

In the case where a user has not yet decided on the set of markers that will be used for screening the segregating population, ResynPy pipeline is providing an additional functionality that is to detect regions of homozygosity (Runs Of Homozygosity; ROH) from a VCF file resulted from the analysis of resequencing data from the elite cultivar line. The output of this process is a tab-delimited file containing the annotation of each marker based on its presence inside a ROH region or not.

#### Input data

##### Mandatory arguments

--vcf: VCF file containing SNPs of the elite line for which a segregating population will be generated.

#### Output data

- `markers_region-anno.tab`



