# Snakemake run on included test data

## Description

Here we describe the step-by-step process used to run the snakemake pipeline on the included test data.

## 

### 1.  Load slurm and miniconda
Note. The commands to do this will be different on your machine. These commands are specific to an HPC with slurm and miniconda modules installed.

```bash
ml slurm/20.02
ml miniconda/4.11.0
```
### 2.  Clone repository
```bash
git clone https://github.com/SansamLab/Process_RepTiming_Snakemake.git
# rename folder with project name
mv Process_RepTiming_Snakemake/ RunTest/
# change directory into root of your project folder
cd RunTest
```
### 3.  Start the snakemake conda environment
#### 3A.  Setup the snakemake conda environment
```bash
# sbatch submits a job script to slurm
# the --wrap option wraps the quoted command string in a sh script and submits
# -f is the location of the environment .yml file. 
## The relative path assumes that you are in the root directory of this repository.
# -p is the path where you want to install this environment
sbatch --mem 16G --wrap="conda env create -f workflow/envs/SnakemakeEnv.yml -p SnakemakeEnv" 
```

#### 3B.  Activate the snakemake conda environment
```bash
conda activate SnakemakeEnv/
```

### 4. Modify the job-specific configuration files.

#### 4A. Modify the config/config.yml file
The config.yml file should be preconfigured for the test data set. Check it to be sure:
![Config File Image](https://github.com/SansamLab/Process_RepTiming_Snakemake/blob/make-test-vignette/resources/configFileImage.png)

#### 4B. Modify the config/samples.csv file
The testSamples.csv file in the config folder has paths to the test fastq files.
![Sample Table Image](https://github.com/SansamLab/Process_RepTiming_Snakemake/blob/make-test-vignette/resources/sampleTableImage.png)

### 5A. Run pipeline with conda environments (Alternative 1)
#### Install necessary conda environments
```
sbatch --mem 32G \
--wrap="\
snakemake \
--cores all \
--use-conda \
--conda-prefix condEnvs/ \
--conda-create-envs-only \
--conda-frontend conda"
```
#### Run pipeline with conda environments
```bash
While within the root directory of the repository clone, enter the following command.
sbatch --constraint=westmere \
--wrap="\
snakemake \
-R \
-j 999 \
--use-conda \
--conda-prefix condEnvs/ \
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

### 5B. Run pipeline with installed modules (Alternative 2)
#### Modify Snakefile with modules installed on your hpc
Each rule in the workflow/Snakefile file has modules listed. These should be changed to match the names of the modules on your hpc. For example:
![rule change example](https://github.com/SansamLab/Process_RepTiming_Snakemake/blob/make-test-vignette/resources/ruleChangeExample.png)

#### Run pipeline with modules installed on hpc
While within the root directory of the repository clone, enter the following command.
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
### 6.  Check output in the results/ directory
Use the tree program to list all of the files with sizes.
```bash
tree -hF results
```

You should get an output like this:
```
results
├── [   0]  aligned/
│   ├── [ 21M]  testInput.bam
│   ├── [1.9M]  testInput.bam.bai
│   ├── [ 21M]  testInput2.bam
│   ├── [1.9M]  testInput2.bam.bai
│   ├── [ 21M]  testSample.bam
│   ├── [1.9M]  testSample.bam.bai
│   ├── [ 21M]  testSample2.bam
│   └── [1.9M]  testSample2.bam.bai
├── [   0]  bedgraphs/
│   ├── [5.6M]  test_Log2Ratios_bedgraphs.gz
│   ├── [2.2M]  test_Quotients_bedgraphs.gz
│   ├── [5.3M]  test_Smoothed_bedgraphs.gz
│   └── [6.3M]  test_ZScores_bedgraphs.gz
├── [   0]  counts/
│   ├── [2.9M]  testInput2_counts.bedgraph
│   ├── [2.9M]  testInput_counts.bedgraph
│   ├── [2.9M]  testSample2_counts.bedgraph
│   └── [2.9M]  testSample_counts.bedgraph
├── [   0]  doubleFilteredBams/
│   ├── [ 17M]  testInput.bam
│   ├── [1.7M]  testInput.bam.bai
│   ├── [ 17M]  testInput2.bam
│   ├── [1.7M]  testInput2.bam.bai
│   ├── [ 17M]  testSample.bam
│   ├── [1.7M]  testSample.bam.bai
│   ├── [ 17M]  testSample2.bam
│   └── [1.7M]  testSample2.bam.bai
├── [   0]  duplicatesMarkedBam/
│   ├── [2.0M]  testInput.bai
│   ├── [ 21M]  testInput.bam
│   ├── [1.9M]  testInput.bam.bai
│   ├── [2.0M]  testInput2.bai
│   ├── [ 21M]  testInput2.bam
│   ├── [1.9M]  testInput2.bam.bai
│   ├── [2.0M]  testSample.bai
│   ├── [ 21M]  testSample.bam
│   ├── [1.9M]  testSample.bam.bai
│   ├── [2.0M]  testSample2.bai
│   ├── [ 21M]  testSample2.bam
│   └── [1.9M]  testSample2.bam.bai
├── [   0]  filteredBams/
│   ├── [ 18M]  testInput.bam
│   ├── [1.7M]  testInput.bam.bai
│   ├── [ 18M]  testInput2.bam
│   ├── [1.7M]  testInput2.bam.bai
│   ├── [ 17M]  testSample.bam
│   ├── [1.7M]  testSample.bam.bai
│   ├── [ 17M]  testSample2.bam
│   └── [1.7M]  testSample2.bam.bai
├── [   0]  logs/
│   ├── [2.1K]  align_reads_with_bwamem.testInput.log
│   ├── [2.1K]  align_reads_with_bwamem.testInput2.log
│   ├── [2.1K]  align_reads_with_bwamem.testSample.log
│   ├── [2.1K]  align_reads_with_bwamem.testSample2.log
│   ├── [ 661]  blacklist_filter_with_bedtools.testInput.log
│   ├── [ 665]  blacklist_filter_with_bedtools.testInput2.log
│   ├── [ 665]  blacklist_filter_with_bedtools.testSample.log
│   ├── [ 669]  blacklist_filter_with_bedtools.testSample2.log
│   ├── [1.5K]  count_reads_in_RT_windows.testInput.log
│   ├── [1.5K]  count_reads_in_RT_windows.testInput2.log
│   ├── [1.6K]  count_reads_in_RT_windows.testSample.log
│   ├── [1.6K]  count_reads_in_RT_windows.testSample2.log
│   ├── [ 26K]  make_bedgraphs.test.log
│   ├── [4.8K]  mark_duplicates_with_picard.testInput.log
│   ├── [4.8K]  mark_duplicates_with_picard.testInput2.log
│   ├── [4.8K]  mark_duplicates_with_picard.testSample.log
│   ├── [4.8K]  mark_duplicates_with_picard.testSample2.log
│   ├── [3.4K]  merge_count_tables.test.log
│   ├── [   0]  picard_metrics/
│   │   ├── [3.5K]  testInput.marked.duplicates_metrics.txt
│   │   ├── [3.5K]  testInput2.marked.duplicates_metrics.txt
│   │   ├── [3.5K]  testSample.marked.duplicates_metrics.txt
│   │   └── [3.5K]  testSample2.marked.duplicates_metrics.txt
│   ├── [6.9K]  process_count_tables.test.log
│   ├── [ 658]  quality_filter_with_bamtools.testInput.log
│   ├── [ 662]  quality_filter_with_bamtools.testInput2.log
│   ├── [ 662]  quality_filter_with_bamtools.testSample.log
│   ├── [ 666]  quality_filter_with_bamtools.testSample2.log
│   ├── [   0]  snakelogs/
│   │   ├── [   0]  align_reads_with_bwamem.testInput.log
│   │   ├── [   0]  align_reads_with_bwamem.testInput2.log
│   │   ├── [   0]  align_reads_with_bwamem.testSample.log
│   │   ├── [   0]  align_reads_with_bwamem.testSample2.log
│   │   ├── [   0]  blacklist_filter_with_bedtools.testInput.log
│   │   ├── [   0]  blacklist_filter_with_bedtools.testInput2.log
│   │   ├── [   0]  blacklist_filter_with_bedtools.testSample.log
│   │   ├── [   0]  blacklist_filter_with_bedtools.testSample2.log
│   │   ├── [   0]  count_reads_in_RT_windows.testInput.log
│   │   ├── [   0]  count_reads_in_RT_windows.testInput2.log
│   │   ├── [   0]  count_reads_in_RT_windows.testSample.log
│   │   ├── [   0]  count_reads_in_RT_windows.testSample2.log
│   │   ├── [   0]  mark_duplicates_with_picard.testInput.log
│   │   ├── [   0]  mark_duplicates_with_picard.testInput2.log
│   │   ├── [   0]  mark_duplicates_with_picard.testSample.log
│   │   ├── [   0]  mark_duplicates_with_picard.testSample2.log
│   │   ├── [   0]  quality_filter_with_bamtools.testInput.log
│   │   ├── [   0]  quality_filter_with_bamtools.testInput2.log
│   │   ├── [   0]  quality_filter_with_bamtools.testSample.log
│   │   ├── [   0]  quality_filter_with_bamtools.testSample2.log
│   │   ├── [   0]  test_make_bedgraphs.log
│   │   ├── [   0]  test_merge_count_tables.log
│   │   ├── [   0]  test_process_count_tables.log
│   │   ├── [   0]  trim_reads_with_fastp.testInput.log
│   │   ├── [   0]  trim_reads_with_fastp.testInput2.log
│   │   ├── [   0]  trim_reads_with_fastp.testSample.log
│   │   └── [   0]  trim_reads_with_fastp.testSample2.log
│   ├── [1.9K]  trim_reads_with_fastp.testInput.log
│   ├── [2.0K]  trim_reads_with_fastp.testInput2.log
│   ├── [2.0K]  trim_reads_with_fastp.testSample.log
│   └── [2.0K]  trim_reads_with_fastp.testSample2.log
├── [   0]  merged/
│   └── [941K]  test_counts.txt
├── [   0]  processed_rse/
│   └── [8.9M]  test_processed_rse.rds
├── [   0]  qc/
│   └── [   0]  fastp_reports/
│       ├── [467K]  testInput.html
│       ├── [122K]  testInput.json
│       ├── [467K]  testInput2.html
│       ├── [122K]  testInput2.json
│       ├── [467K]  testSample.html
│       ├── [121K]  testSample.json
│       ├── [467K]  testSample2.html
│       └── [121K]  testSample2.json
├── [   0]  rse/
│   └── [356K]  test_rse.rds
└── [   0]  trimmed/
    ├── [7.8M]  testInput2_trimmed_R1.fastq.gz
    ├── [8.5M]  testInput2_trimmed_R2.fastq.gz
    ├── [7.8M]  testInput_trimmed_R1.fastq.gz
    ├── [8.5M]  testInput_trimmed_R2.fastq.gz
    ├── [7.8M]  testSample2_trimmed_R1.fastq.gz
    ├── [8.5M]  testSample2_trimmed_R2.fastq.gz
    ├── [7.8M]  testSample_trimmed_R1.fastq.gz
    └── [8.5M]  testSample_trimmed_R2.fastq.gz

15 directories, 121 files
```
