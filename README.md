![Release](https://img.shields.io/github/v/release/SansamLab/Process_RepTiming_Snakemake?include_prereleases)
![ReleaseDate](https://img.shields.io/github/release-date-pre/SansamLab/Process_RepTiming_Snakemake?include_prereleases)
![Size](https://img.shields.io/github/repo-size/SansamLab/Process_RepTiming_Snakemake)
![License](https://img.shields.io/github/license/SansamLab/Process_RepTiming_Snakemake)
![LastCommit](https://img.shields.io/github/last-commit/SansamLab/Process_RepTiming_Snakemake)
![Downloads](https://img.shields.io/github/downloads/SansamLab/Process_RepTiming_Snakemake/total)
![OpenIssues](https://img.shields.io/github/issues-raw/SansamLab/Process_RepTiming_Snakemake)


# Process_RepTiming_Snakemake

## Project Description:

Process_RepTiming_Snakemake describes how to process short whole-genome sequening reads from G1 and S phase cells into replication timing values. Each of the individual data processing steps are described, which enables step-by-step processing to be done. Alternatively, a Snakemake pipeline with clearly defined dependencies and Anaconda environments is also provided so that the data processing pipeline can be automated.

## Table of contents:
* [Examples](https://github.com/SansamLab/Process_RepTiming_Snakemake/blob/main/README.md#examples)
  * [Run Snakemake pipeline with sample data included in repository](https://github.com/SansamLab/Process_RepTiming_Snakemake/blob/main/RunSnakemakeOnTestData.md)
  * [Run Snakemake pipeline with Siefert, 2017 data.]()
* [Requirements](https://github.com/SansamLab/Process_RepTiming_Snakemake/blob/main/README.md#requirements)
* [Description of individual steps in pipeline](https://github.com/SansamLab/Process_HiC_SnakeMake/edit/main/README.md#description-of-individual-steps-in-pipeline)
  * [1.  trim_reads_with_fastp](https://github.com/SansamLab/Process_RepTiming_Snakemake#1--trim_reads_with_fastp)
  * [2.  align_reads_with_bwamem](https://github.com/SansamLab/Process_RepTiming_Snakemake#2--run_bwa_mem)
  * [3.  mark_duplicates_with_picard](https://github.com/SansamLab/Process_RepTiming_Snakemake#3--mark_duplicates_with_picard)
  * [4.  quality_filter_with_bamtools](https://github.com/SansamLab/Process_RepTiming_Snakemake#4--quality_filter_with_bamtools)
  * [5.  blacklist_filter_with_bedtools](https://github.com/SansamLab/Process_RepTiming_Snakemake#5--blacklist_filter_with_bedtools)
  * [6.  count_reads_in_RT_windows](https://github.com/SansamLab/Process_RepTiming_Snakemake#6--count_reads_in_rt_windows)
  * [7.  Merge count tables](https://github.com/SansamLab/Process_RepTiming_Snakemake#7--merge-count-tables)
  * [8.  Process count tables](https://github.com/SansamLab/Process_RepTiming_Snakemake#8--process-count-tables)
  * [9.  Make bedgraphs](https://github.com/SansamLab/Process_RepTiming_Snakemake#9--make-bedgraphs)
* [Step-by-step instructions on running Snakemake pipeline:](https://github.com/SansamLab/Process_RepTiming_Snakemake#step-by-step-instructions-on-running-snakemake-pipeline)
  * [1.  Load slurm and miniconda](https://github.com/SansamLab/Process_RepTiming_Snakemake#1--load-slurm-and-miniconda)
  * [2.  Clone repository](https://github.com/SansamLab/Process_RepTiming_Snakemake#2--clone-repository)
  * [3.  Start the conda environment](https://github.com/SansamLab/Process_RepTiming_Snakemake#3--start-the-conda-environment)
    * [3A.  FIRST TIME ONLY:  Setup conda environment](https://github.com/SansamLab/Process_RepTiming_Snakemake#3a--first-time-only--setup-conda-environment)
    * [3B.  Activate conda environment](https://github.com/SansamLab/Process_RepTiming_Snakemake#3b--activate-conda-environment)
  * [4.  Modify the job-specific configuration files.](https://github.com/SansamLab/Process_RepTiming_Snakemake#4--modify-the-job-specific-configuration-files)
    * [4A.  Modify the config/config.yml file](https://github.com/SansamLab/Process_RepTiming_Snakemake#4a--modify-the-configconfigyml-file)
    * [4B.  Modify the config/samples.csv file](https://github.com/SansamLab/Process_RepTiming_Snakemake#4b--modify-the-configsamplescsv-file)
    * [4C.  IF SLURM RESOURCE CHANGES ARE NEEDED. Modify the config/cluster_config.yml file](https://github.com/SansamLab/Process_RepTiming_Snakemake#4c--if-slurm-resource-changes-are-needed-modify-the-configcluster_configyml-file)
  * [5.  Do a dry run](https://github.com/SansamLab/Process_RepTiming_Snakemake#5--do-a-dry-run)
  * [6.  Make a DAG diagram](https://github.com/SansamLab/Process_RepTiming_Snakemake#6--make-a-dag-diagram)
  * [7.  Run on cluster with slurm](https://github.com/SansamLab/Process_RepTiming_Snakemake#7--run-on-cluster-with-slurm)
    * [7A.  Use conda environments](https://github.com/SansamLab/Process_RepTiming_Snakemake/blob/main/README.md#7a-use-conda-environments)
    * [7B.  Use modules](https://github.com/SansamLab/Process_RepTiming_Snakemake/blob/main/README.md#7b-use-environment-modules)
* [Output structure](https://github.com/SansamLab/Process_RepTiming_Snakemake/blob/main/README.md#output-structure)
* [References](https://github.com/SansamLab/Process_RepTiming_Snakemake/blob/main/README.md#references)

## Examples
[Run Snakemake pipeline with sample data included in repository](https://github.com/SansamLab/Process_RepTiming_Snakemake/blob/main/RunSnakemakeOnTestData.md)

## Requirements

### General requirements:
**Requirement**|**Description**
:-----:|:-----:
Raw RT data|separate .fastq files from paired-end whole genome sequencing of S phase and matched G1 cells
Indexed genome|your genome of interest indexed for bwa
RT windows|.bed file with coordinates of replication timing windows. RT windows file for the zebrafish GRCz11 genome is provided.
Blacklisted regions|.bed file with blacklisted regions. Blacklisted file provided for zebrafish GRCz11 is provided.
Greylisted regions|.bed file with greylisted regions. Greylisted file for zebrafish GRCz11 that works well for the Tab5 strain is provided.
Software|fastp/0.23.2, bwa/0.7.15, samtools/1.14, picard/2.21.2, bamtools/2.5.1, bedtools/2.30.0, R/4.1.2-mkl
R Packages|GenomicRanges, SummarizedExperiment, pspline, bamsignals, matrixStats

### Snakemake pipeline requirements:
**Requirement**|**Description**
:-----:|:-----:
Software|miniconda/4.11.0 (for Conda environments), python/3.10.2, numpy/1.2.2, pandas/1.4.1, snakemake-minimal/7.3.1
Snakemake config file|.csv file with snakemake parameters in the config/ directory
Snakemake cluster config file|.csv file with cluster configuration parameters in the config/directory

## Description of individual steps in pipeline:

The data are processed according to Siefert, 2017 and Siefert, 2018 with some modifications. Paired-end sequencing reads are trimmed with fastp and then aligned to the genome using BWA-mem. Aligned reads are filtered using picard tools, bamtools, and bedtools. Reads are counted in user selected replication timing windows using the R bamsignals package with a custom R script. The S/G1 quotients are calculated for each window using a custom R script and then smoothed with a with a cubic smoothing spline through a custom R script and the pspline R package. The smoothed quotients are then transformed into log2 values or z-scores. If replicate data is provided, median values are calculated. Individual and median data for all samples are provided in an R RangedSummarizedExperiment (.rds) file and in individual .bedgraph files. 

![DAG of Test Pipeline](https://github.com/SansamLab/Process_RepTiming_Snakemake/blob/main/resources/dag.png)

### 1.  trim_reads_with_fastp
```bash
fastp \
  -i {input.fq1} \
  -I {input.fq2} \
  -o {output.trimmed1} \
  -O {output.trimmed2} \
  -h {output.fastp_report} \
  --json {output.fastp_json} \
  -R "{wildcards.sample}" \
  -w {params.threads}
```
### 2.  run_bwa_mem
```bash
# align pair of .fastq files to the genome and convert the .sam output to .bam
bwa mem \
  -M \
  -t {params.threads} \
  {params.genome} \
  {input.R1} {input.R2} | \
samtools sort \
  -@ {params.threads} > {output.bam}
    
samtools index \
  -@ {params.threads} \
    {output.bam} > {output.bai}
```
### 3.  mark_duplicates_with_picard
```bash
JAVA_MEM_OPTS={params.memory}

RAM=$(echo $JAVA_MEM_OPTS | rev | cut -c2- | rev)

picard MarkDuplicates \
  COMPRESSION_LEVEL=9 \
  VALIDATION_STRINGENCY=LENIENT \
  MAX_RECORDS_IN_RAM=$((200000*RAM)) \
  CREATE_INDEX=true \
  I={input.bam} \
  O={output.marked_bam} \
  M={output.marked_metrics}

samtools index {output.marked_bam}
```

### 4.  quality_filter_with_bamtools
```bash
bamtools filter \
  -in {input.marked_bam} \
  -forceCompression \
  -mapQuality '{params.mapQ}' \
  -isDuplicate false \
  -isFailedQC false \
  -isMapped true \
  -isMateMapped true \
  -isPaired true \
  -isPrimaryAlignment true \
  -isProperPair true \
  -out {output.filtered_bam}

samtools index {output.filtered_bam}
```

### 5.  blacklist_filter_with_bedtools
```bash
bedtools intersect \
  -v \
  -split \
  -abam {input.filtered_bam} \
  -b {params.HiCount_BED_FILENAME} > {input.filtered_bam}_FILTERED1.BAM

bedtools intersect \
  -v \
  -split \
  -f 0.3 \
  -abam {input.filtered_bam}_FILTERED1.BAM \
  -b {params.Sat_BED_FILENAME} > {output.doubleFiltered_bam}

rm {input.filtered_bam}_FILTERED1.BAM

samtools index {output.doubleFiltered_bam}
```
### 6.  count_reads_in_RT_windows
```bash
Rscript \
  workflow/scripts/CountReadsOverBed_ver02.R \
    {params.RT_windows} \
    {input.doubleFiltered_bam} \
    {output.counts}          
```
### 7.  Merge count tables
```bash
SAMPLES=( {params.sample_list} )
BEDGRAPHS=( {params.counts_bedgraphs_list} )
TEMP_COUNTS=( {params.counts_temp_list} )
mkdir results/temp_merge
for i in ${{!BEDGRAPHS[@]}}; do 
 echo -e "\t\t\t${{SAMPLES[i]}}" > ${{TEMP_COUNTS[i]}}.header
 cut -f4 ${{BEDGRAPHS[i]}} > ${{TEMP_COUNTS[i]}}.counts
 cat ${{TEMP_COUNTS[i]}}.header ${{TEMP_COUNTS[i]}}.counts > ${{TEMP_COUNTS[i]}}
 rm ${{TEMP_COUNTS[i]}}.header
 rm ${{TEMP_COUNTS[i]}}.counts
done 
paste {params.counts_temp_list} > {output.merged}
rm -rf results/temp_merge/
Rscript workflow/scripts/MakeCountsRSE_ver01.R \
  {output.merged} \
  {params.samples_table} \
  {params.RT_windows} \
  {output.rse_counts}
```
### 8.  Process count tables
```bash
Rscript workflow/scripts/CalculateQuotientsSmoothScale_ver01.R \
  {input.rse_counts} \
  {output.rse_processed}
```

### 9.  Make bedgraphs
```bash
mkdir Log2Ratios
mkdir ZScores
mkdir Smoothed
mkdir Quotients
Rscript workflow/scripts/Generate_RT_Bedgraphs_ver01.R \
  {input.rse_processed} \
  "Log2Ratios"
Rscript workflow/scripts/Generate_RT_Bedgraphs_ver01.R \
  {input.rse_processed} \
  "ZScores"
Rscript workflow/scripts/Generate_RT_Bedgraphs_ver01.R \
  {input.rse_processed} \
  "Smoothed"
Rscript workflow/scripts/Generate_RT_Bedgraphs_ver01.R \
  {input.rse_processed} \
  "Quotients"
tar -czvf {output.Log2Ratios_bedgraphs} Log2Ratios/
tar -czvf {output.ZScores_bedgraphs} ZScores/
tar -czvf {output.Smoothed_bedgraphs} Smoothed/
tar -czvf {output.Quotients_bedgraphs} Quotients/
rm -rf Log2Ratios
rm -rf ZScores
rm -rf Smoothed
rm -rf Quotients
```

## Step-by-step instructions on running Snakemake pipeline:

### 1.  Load slurm and miniconda
Note. The commands to do this will be different on your machine. These commands are specific to an HPC using slurm with these modules installed.

```bash
ml slurm/20.02
ml miniconda/4.11.0
```
### 2.  Clone repository
```bash
git clone https://github.com/SansamLab/Process_RepTiming_Snakemake.git
# rename folder with project name
mv Process_RepTiming_Snakemake/ My_RT_Project_Folder/
# change directory into root of your project folder
cd My_RT_Project_Folder
```
### 3.  Start the conda environment
#### 3A.  FIRST TIME ONLY:  Setup conda environment with snakemake
```bash
# -f is the location of the environment .yml file. 
## The relative path assumes that you are in the root directory of this repository.
# -p is the path where you want to install this environment
conda env create -f workflow/envs/SnakemakeEnv.yml -p /s/sansam-lab/SnakemakeEnv 
```

#### 3B.  Activate conda environment with snakemake
```bash
conda activate /s/sansam-lab/SnakemakeEnv
```

### 4.  Modify the job-specific configuration files.
#### 4A.  Modify the config/config.yml file

You must enter the following:
* samples_table:
  * path and name of the .csv file that lists sample names and paths to all fastq files
  * example:  "config/testSamples.csv"
* fastp_cpus:
  * number of processors that fastp will use to trim your sequencing reads
  * example:  8
* bwa_cpus:
  * number of processors that bwa mem will use to align reads to genome
  * example:  12
* picard_memory:
  * amount of memory allocated to picard for marking duplicates
  * example:  "24G"
* bwa_genome:
  * location of bwa indexed genome for the alignment
  * example:  "/Volumes/hts_core/Shared/zebrafish/danRer11/noAlts/danRer11_noAlts.fa"
* bamtools_filter_mapQ:
  * mapQ score cutoff for quality filtering reads
  * example:  ">10"
* bedtools_intersect_blacklist:
  * .bed file with blacklisted genomic regions for filtering out reads
  * example:  "resources/Satellites.bed"
    * for zebrafish GRCz11
    * file has regions with certain repeat types that we find are associated with anomalously high read counts
* bedtools_intersect_greylist:
  * .bed file with greylisted genomic regions for filtering out reads
  * example:  "resources/HiCountWindows.bed"
    * for zebrafish GRCz11
    * file has regions that have anomalously high read counts in multiple Tab5 whole genome sequencing runs
* RT_windows:
  * .bed file with regions used for calculating replication timing values
  * example:  "resources/RTWindows_danRer11_noAlts_ver01.bed"
    * for zebrafish GRCz11

#### 4B.  Modify the config/samples.csv file

The samples.csv file in the config folder has paths to the test fastq files. You must replace those paths with those for your own fastq files. The first column of each row is the sample name. This name will be used for all output files.

#### 4C.  IF SLURM RESOURCE CHANGES ARE NEEDED. Modify the config/cluster_config.yml file

CPU and memory requests for each rule in the pipeline are detailed in this file. If you are using SLURM, you may need to alter this file to fit your needs/system.

### 5.  Do a dry run.
A dry run produces a text output showing exactly what commands will be executed. Look this over carefully before submitting the full job. It is normal to see warnings about changes made to the code, input, and params.
```bash
snakemake -npr
```

### 6.  Make a DAG diagram.
```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

### 7.  Run on cluster with slurm.
This snakemake pipeline could be executed without slurm, but if an hpc with slurm is used, the following will start the pipeline with the parameters defined in the config/cluster_config.yml file.

#### 7A. Use conda environments
If conda is to be used for rule-specific environments, you may find it useful to create the environments first. Running 'snakemake' with the '--conda-create-envs-only' option will create the environments without running the pipeline. The '--conda-prefix' option is used to set a directory in which the ‘conda’ and ‘conda-archive’ directories are created. This directory may be changed to a stable or shared location.
```bash
sbatch --mem 32G \
--wrap="\
snakemake \
--cores all \
--use-conda \
--conda-prefix ../condEnvs/ \
--conda-create-envs-only \
--conda-frontend conda"
```

Once the environments are setup, you may execute pipeline with conda environments using the following command:
```bash
sbatch --constraint=westmere \
--wrap="\
snakemake \
-R \
-j 999 \
--use-conda \
--conda-prefix ../condEnvs/ \
--conda-frontend conda \
--latency-wait 100 \
--cluster-config config/cluster_config.yml \
--cluster '\
sbatch \
-A {cluster.account} \
-p {cluster.partition} \
--cpus-per-task {cluster.cpus-per-task}  \
--mem {cluster.mem} \
--output {cluster.output}'"
```

#### 7B. Use environment modules.
Rather than using conda environments, you may prefer to use modules installed on your computing cluster. These modules are defined for each rule in 'workflow/Snakefile'. This must be customized for your environment, and you must modify the Snakefile yourself.

To execute the pipeline with environment modules, enter the following:
```bash
sbatch --constraint=westmere \
--wrap="\
snakemake \
-R \
-j 999 \
--use-envmodules \
--latency-wait 100 \
--cluster-config config/cluster_config.yml \
--cluster '\
sbatch \
-A {cluster.account} \
-p {cluster.partition} \
--cpus-per-task {cluster.cpus-per-task}  \
--mem {cluster.mem} \
--output {cluster.output}'"
```

### 8.  Check results, and when finished, exit environment.
The results will be saved to the "results" folder. Look over log files generated in either the logs/ or logs/snakelogs folders (depending on whether slurm was used).
```bash
conda deactivate
```

## Output structure:

**Folder**|**Description**
:-----:|:-----:
logs|Log file for each snakemake rule executed
qc|Quality control files generated by fastp
trimmed|.fastq files trimmed by fastq
aligned|.bam files aligned by bwamem
duplicatesMarkedBam|.bam files with duplicates marked by picardtools
filteredBams|.bam files filtered for read quality
doubleFilteredBams|.bam files also filtered using blacklisted and greylisted regions
counts|.bedgraph files with counts for each sample
merged|All counts in a .txt tab delimited table
rse|.rds file with RangedSummarized Experiment Object with Counts
processed\_rse|.rds file with RangedSummarizedExperiment Object with quotients, smoothed, Log2, and ZScores for all S phase samples. Medians values for replicates are also included
bedgraphs|Zipped bedgraphs with calculated values from processed\_rse

## References:

### Replication timing data processing

Siefert JC, Georgescu C, Wren JD, Koren A, Sansam CL. DNA replication timing during development anticipates transcriptional programs and parallels enhancer activation. Genome Res. 2017;27(8):1406-1416.

Siefert JC, Clowdus EA, Goins D, Koren A, Sansam CL. Profiling dna replication timing using zebrafish as an in vivo model system. JoVE. 2018;(134):57146.

### fastp

Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018;34(17):i884-i890. doi:10.1093/bioinformatics/bty560

### Samtools

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of samtools and bcftools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008

### BWA

Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. ArXiv:1303.3997 [q-Bio]. http://arxiv.org/abs/1303.3997

### Picard tools

http://broadinstitute.github.io/picard/

### bamtools

Barnett DW, Garrison EK, Quinlan AR, Stromberg MP, Marth GT. BamTools: a C++ API and toolkit for analyzing and managing BAM files. Bioinformatics. 2011;27(12):1691-1692.

### bedtools

Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010;26(6):841-842.

### R

R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.
URL https://www.R-project.org/.

### bamsignals

Alessandro Mammana and Johannes Helmuth (2021). bamsignals: Extract read count signals from bam files. R package version 1.26.0.
https://github.com/lamortenera/bamsignals

### pspline

S original by Jim Ramsey. R port by Brian Ripley <ripley@stats.ox.ac.uk>. (2022). pspline: Penalized Smoothing Splines. R package version 1.0-19. https://CRAN.R-project.org/package=pspline

### SummarizedExperiment

Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0. https://bioconductor.org/packages/SummarizedExperiment

### GenomicRanges

Lawrence M, Huber W, Pag\`es H, Aboyoun P, Carlson M, et al. (2013) Software for Computing and Annotating Genomic Ranges. PLoS Comput Biol 9(8): e1003118. doi:10.1371/journal.pcbi.1003118

### matrixstats

Henrik Bengtsson (2021). matrixStats: Functions that Apply to Rows and Columns of Matrices (and to Vectors). R package version 0.61.0. https://CRAN.R-project.org/package=matrixStats

### Snakemake
Köster, J., & Rahmann, S. (2012). Snakemake—A scalable bioinformatics workflow engine. Bioinformatics (Oxford, England), 28(19), 2520–2522. https://doi.org/10.1093/bioinformatics/bts480

### Python
Anaconda Software Distribution. (2020). Anaconda Documentation. Anaconda Inc. Retrieved from https://docs.anaconda.com/
