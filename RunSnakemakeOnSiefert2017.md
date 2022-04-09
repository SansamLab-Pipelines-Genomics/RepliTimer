# Snakemake run on Siefert 2017

## Description

Here we describe the step-by-step process used to run the snakemake pipeline on the Siefert 2017 data, which can be obtained from SRA.

### 1.  Transfer .fastq files from SRA.
This package lacks the functionality to do this. We recommend using the fasterq-dump to transfer the files from SRA. In this example, we create a directory in which the .fastq files are transferred to a subdirectory called "fastqs". The locations of the .fastq files defined in the "Siefert_Samples.csv" provided in this repo as "../fastqs/SRR4036047_1.fastq".

```bash
# load the module with fasterq-dump. note:  this will likely differ on your system
ml ncbi_sra
# create directory for the entire snakemake run
mkdir Siefert2017
cd Siefert2017/
# create a subdirectory to which the .fastq files will be transferred
mkdir fastqs
cd fastqs
# transfer each pair of files for the following SRR entries. 
# In this example we use "sbatch --wrap" to run the fasterq-dump command on our hpc.
SRAIDS=( "SRR4036047" \
"SRR4036048" \
"SRR4036049" \
"SRR4036050" \
"SRR4036051" \
"SRR4036052" \
"SRR4036059" \
"SRR4036060" \
"SRR4036061" \
"SRR4036062" )
for t in ${SRAIDS[@]}; do
  sbatch --cpus-per-task 12 --wrap="fasterq-dump --split-files --threads 12 $t"
  done
# change back to the Siefert2017/ directory
cd ..
```

### 2.  Load slurm and miniconda
Note. The commands to do this will be different on your machine. These commands are specific to an HPC with slurm and miniconda modules installed.

```bash
ml slurm/20.02
ml miniconda/4.11.0
```
### 3.  Clone repository
```bash
git clone https://github.com/SansamLab/RepliTimer.git
# rename folder with project name
mv RepliTimer/ Siefert2017_RT_Project_Folder/
# change directory into root of your project folder
cd Siefert2017_RT_Project_Folder
```
### 4.  Start the snakemake conda environment
#### 4A.  Setup the snakemake conda environment
```bash
# sbatch submits a job script to slurm
# the --wrap option wraps the quoted command string in a sh script and submits
# -f is the location of the environment .yml file. 
## The relative path assumes that you are in the root directory of this repository.
# -p is the path where you want to install this environment
sbatch --mem 16G --wrap="conda env create -f workflow/envs/SnakemakeEnv.yml -p SnakemakeEnv" 
```

#### 4B.  Activate the snakemake conda environment
```bash
conda activate SnakemakeEnv/
```

### 5. Modify the job-specific configuration files.

#### 5A. Modify the config/config.yml file
The config.yml file is preconfigured for the test data set, so the path to the .csv file must be changed.
![Config File Image](https://github.com/SansamLab/RepliTimer/blob/main/resources/SiefertSamplesTableImage.png)

#### 5B. Modify the config/Siefert_Samples.csv file
The Siefert_Samples.csv file in the config folder has relative paths to the fastq files transferred from SRA. If the path differs in your system update it in Siefert_Samples.csv.
![Sample Table Image](https://github.com/SansamLab/RepliTimer/blob/main/resources/SiefertSamplesImage.png)

### 6A. Run pipeline with conda environments (Alternative 1)
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

### 6B. Run pipeline with installed modules (Alternative 2)
#### Modify Snakefile with modules installed on your hpc
Each rule in the workflow/Snakefile file has modules listed. These should be changed to match the names of the modules on your hpc. For example:
![rule change example](https://github.com/SansamLab/RepliTimer/blob/main/resources/ruleChangeExample.png)

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
### 7.  Check output in the results/ directory
Use the tree program to list all of the files with sizes.
```bash
tree -hF results
```

You should get an output like this:
```
results
├── [   0]  aligned/
│   ├── [3.0G]  SRR4036047_28hpf_G1.bam
│   ├── [1.1M]  SRR4036047_28hpf_G1.bam.bai
│   ├── [3.2G]  SRR4036048_28hpf_S.bam
│   ├── [1.1M]  SRR4036048_28hpf_S.bam.bai
│   ├── [3.5G]  SRR4036049_28hpf_G1.bam
│   ├── [1.1M]  SRR4036049_28hpf_G1.bam.bai
│   ├── [4.6G]  SRR4036050_28hpf_S.bam
│   ├── [1.5M]  SRR4036050_28hpf_S.bam.bai
│   ├── [8.1G]  SRR4036051_28hpf_G1.bam
│   ├── [3.7M]  SRR4036051_28hpf_G1.bam.bai
│   ├── [9.2G]  SRR4036052_28hpf_S.bam
│   ├── [3.8M]  SRR4036052_28hpf_S.bam.bai
│   ├── [1.5G]  SRR4036059_Bud.bam
│   ├── [1.1M]  SRR4036059_Bud.bam.bai
│   ├── [1.8G]  SRR4036060_Bud.bam
│   ├── [1.1M]  SRR4036060_Bud.bam.bai
│   ├── [4.3G]  SRR4036061_ZTF_G1.bam
│   ├── [1.4M]  SRR4036061_ZTF_G1.bam.bai
│   ├── [4.3G]  SRR4036062_ZTF_S.bam
│   └── [1.3M]  SRR4036062_ZTF_S.bam.bai
├── [   0]  bedgraphs/
│   ├── [ 17M]  Siefert2017_Log2Ratios_bedgraphs.gz
│   ├── [ 16M]  Siefert2017_Quotients_bedgraphs.gz
│   ├── [ 16M]  Siefert2017_Smoothed_bedgraphs.gz
│   └── [ 18M]  Siefert2017_ZScores_bedgraphs.gz
├── [   0]  counts/
│   ├── [3.1M]  SRR4036047_28hpf_G1_counts.bedgraph
│   ├── [3.1M]  SRR4036048_28hpf_S_counts.bedgraph
│   ├── [3.1M]  SRR4036049_28hpf_G1_counts.bedgraph
│   ├── [3.1M]  SRR4036050_28hpf_S_counts.bedgraph
│   ├── [3.2M]  SRR4036051_28hpf_G1_counts.bedgraph
│   ├── [3.2M]  SRR4036052_28hpf_S_counts.bedgraph
│   ├── [3.1M]  SRR4036059_Bud_counts.bedgraph
│   ├── [3.1M]  SRR4036060_Bud_counts.bedgraph
│   ├── [3.1M]  SRR4036061_ZTF_G1_counts.bedgraph
│   └── [3.1M]  SRR4036062_ZTF_S_counts.bedgraph
├── [   0]  doubleFilteredBams/
│   ├── [2.3G]  SRR4036047_28hpf_G1.bam
│   ├── [1.1M]  SRR4036047_28hpf_G1.bam.bai
│   ├── [2.5G]  SRR4036048_28hpf_S.bam
│   ├── [1.1M]  SRR4036048_28hpf_S.bam.bai
│   ├── [2.7G]  SRR4036049_28hpf_G1.bam
│   ├── [1.1M]  SRR4036049_28hpf_G1.bam.bai
│   ├── [3.5G]  SRR4036050_28hpf_S.bam
│   ├── [1.3M]  SRR4036050_28hpf_S.bam.bai
│   ├── [5.9G]  SRR4036051_28hpf_G1.bam
│   ├── [3.3M]  SRR4036051_28hpf_G1.bam.bai
│   ├── [6.7G]  SRR4036052_28hpf_S.bam
│   ├── [3.4M]  SRR4036052_28hpf_S.bam.bai
│   ├── [1.1G]  SRR4036059_Bud.bam
│   ├── [1.1M]  SRR4036059_Bud.bam.bai
│   ├── [1.4G]  SRR4036060_Bud.bam
│   ├── [1.1M]  SRR4036060_Bud.bam.bai
│   ├── [3.4G]  SRR4036061_ZTF_G1.bam
│   ├── [1.2M]  SRR4036061_ZTF_G1.bam.bai
│   ├── [3.3G]  SRR4036062_ZTF_S.bam
│   └── [1.2M]  SRR4036062_ZTF_S.bam.bai
├── [   0]  duplicatesMarkedBam/
│   ├── [3.6M]  SRR4036047_28hpf_G1.bai
│   ├── [3.0G]  SRR4036047_28hpf_G1.bam
│   ├── [1.1M]  SRR4036047_28hpf_G1.bam.bai
│   ├── [3.6M]  SRR4036048_28hpf_S.bai
│   ├── [3.3G]  SRR4036048_28hpf_S.bam
│   ├── [1.1M]  SRR4036048_28hpf_S.bam.bai
│   ├── [3.6M]  SRR4036049_28hpf_G1.bai
│   ├── [3.5G]  SRR4036049_28hpf_G1.bam
│   ├── [1.2M]  SRR4036049_28hpf_G1.bam.bai
│   ├── [3.7M]  SRR4036050_28hpf_S.bai
│   ├── [4.6G]  SRR4036050_28hpf_S.bam
│   ├── [1.7M]  SRR4036050_28hpf_S.bam.bai
│   ├── [3.9M]  SRR4036051_28hpf_G1.bai
│   ├── [8.1G]  SRR4036051_28hpf_G1.bam
│   ├── [3.7M]  SRR4036051_28hpf_G1.bam.bai
│   ├── [3.9M]  SRR4036052_28hpf_S.bai
│   ├── [9.2G]  SRR4036052_28hpf_S.bam
│   ├── [3.8M]  SRR4036052_28hpf_S.bam.bai
│   ├── [3.2M]  SRR4036059_Bud.bai
│   ├── [1.5G]  SRR4036059_Bud.bam
│   ├── [1.1M]  SRR4036059_Bud.bam.bai
│   ├── [3.3M]  SRR4036060_Bud.bai
│   ├── [1.9G]  SRR4036060_Bud.bam
│   ├── [1.1M]  SRR4036060_Bud.bam.bai
│   ├── [3.8M]  SRR4036061_ZTF_G1.bai
│   ├── [4.3G]  SRR4036061_ZTF_G1.bam
│   ├── [1.6M]  SRR4036061_ZTF_G1.bam.bai
│   ├── [3.8M]  SRR4036062_ZTF_S.bai
│   ├── [4.3G]  SRR4036062_ZTF_S.bam
│   └── [1.5M]  SRR4036062_ZTF_S.bam.bai
├── [   0]  filteredBams/
│   ├── [2.4G]  SRR4036047_28hpf_G1.bam
│   ├── [1.1M]  SRR4036047_28hpf_G1.bam.bai
│   ├── [2.6G]  SRR4036048_28hpf_S.bam
│   ├── [1.1M]  SRR4036048_28hpf_S.bam.bai
│   ├── [2.8G]  SRR4036049_28hpf_G1.bam
│   ├── [1.1M]  SRR4036049_28hpf_G1.bam.bai
│   ├── [3.7G]  SRR4036050_28hpf_S.bam
│   ├── [1.3M]  SRR4036050_28hpf_S.bam.bai
│   ├── [6.2G]  SRR4036051_28hpf_G1.bam
│   ├── [3.4M]  SRR4036051_28hpf_G1.bam.bai
│   ├── [7.1G]  SRR4036052_28hpf_S.bam
│   ├── [3.5M]  SRR4036052_28hpf_S.bam.bai
│   ├── [1.2G]  SRR4036059_Bud.bam
│   ├── [1.1M]  SRR4036059_Bud.bam.bai
│   ├── [1.5G]  SRR4036060_Bud.bam
│   ├── [1.1M]  SRR4036060_Bud.bam.bai
│   ├── [3.5G]  SRR4036061_ZTF_G1.bam
│   ├── [1.2M]  SRR4036061_ZTF_G1.bam.bai
│   ├── [3.4G]  SRR4036062_ZTF_S.bam
│   └── [1.2M]  SRR4036062_ZTF_S.bam.bai
├── [   0]  logs/
│   ├── [ 42K]  align_reads_with_bwamem.SRR4036047_28hpf_G1.log
│   ├── [ 44K]  align_reads_with_bwamem.SRR4036048_28hpf_S.log
│   ├── [ 48K]  align_reads_with_bwamem.SRR4036049_28hpf_G1.log
│   ├── [ 62K]  align_reads_with_bwamem.SRR4036050_28hpf_S.log
│   ├── [127K]  align_reads_with_bwamem.SRR4036051_28hpf_G1.log
│   ├── [146K]  align_reads_with_bwamem.SRR4036052_28hpf_S.log
│   ├── [ 16K]  align_reads_with_bwamem.SRR4036059_Bud.log
│   ├── [ 23K]  align_reads_with_bwamem.SRR4036060_Bud.log
│   ├── [ 66K]  align_reads_with_bwamem.SRR4036061_ZTF_G1.log
│   ├── [ 65K]  align_reads_with_bwamem.SRR4036062_ZTF_S.log
│   ├── [1.1K]  blacklist_filter_with_bedtools.SRR4036047_28hpf_G1.log
│   ├── [1.1K]  blacklist_filter_with_bedtools.SRR4036048_28hpf_S.log
│   ├── [1.1K]  blacklist_filter_with_bedtools.SRR4036049_28hpf_G1.log
│   ├── [1.1K]  blacklist_filter_with_bedtools.SRR4036050_28hpf_S.log
│   ├── [1.1K]  blacklist_filter_with_bedtools.SRR4036051_28hpf_G1.log
│   ├── [1.1K]  blacklist_filter_with_bedtools.SRR4036052_28hpf_S.log
│   ├── [1.1K]  blacklist_filter_with_bedtools.SRR4036059_Bud.log
│   ├── [1.1K]  blacklist_filter_with_bedtools.SRR4036060_Bud.log
│   ├── [1.1K]  blacklist_filter_with_bedtools.SRR4036061_ZTF_G1.log
│   ├── [1.1K]  blacklist_filter_with_bedtools.SRR4036062_ZTF_S.log
│   ├── [2.1K]  count_reads_in_RT_windows.SRR4036047_28hpf_G1.log
│   ├── [1.9K]  count_reads_in_RT_windows.SRR4036048_28hpf_S.log
│   ├── [2.1K]  count_reads_in_RT_windows.SRR4036049_28hpf_G1.log
│   ├── [2.1K]  count_reads_in_RT_windows.SRR4036050_28hpf_S.log
│   ├── [2.1K]  count_reads_in_RT_windows.SRR4036051_28hpf_G1.log
│   ├── [2.1K]  count_reads_in_RT_windows.SRR4036052_28hpf_S.log
│   ├── [2.0K]  count_reads_in_RT_windows.SRR4036059_Bud.log
│   ├── [2.0K]  count_reads_in_RT_windows.SRR4036060_Bud.log
│   ├── [2.1K]  count_reads_in_RT_windows.SRR4036061_ZTF_G1.log
│   ├── [2.1K]  count_reads_in_RT_windows.SRR4036062_ZTF_S.log
│   ├── [ 27K]  make_bedgraphs.Siefert2017.log
│   ├── [ 14K]  mark_duplicates_with_picard.SRR4036047_28hpf_G1.log
│   ├── [ 14K]  mark_duplicates_with_picard.SRR4036048_28hpf_S.log
│   ├── [ 15K]  mark_duplicates_with_picard.SRR4036049_28hpf_G1.log
│   ├── [ 19K]  mark_duplicates_with_picard.SRR4036050_28hpf_S.log
│   ├── [ 33K]  mark_duplicates_with_picard.SRR4036051_28hpf_G1.log
│   ├── [ 36K]  mark_duplicates_with_picard.SRR4036052_28hpf_S.log
│   ├── [9.3K]  mark_duplicates_with_picard.SRR4036059_Bud.log
│   ├── [ 10K]  mark_duplicates_with_picard.SRR4036060_Bud.log
│   ├── [ 19K]  mark_duplicates_with_picard.SRR4036061_ZTF_G1.log
│   ├── [ 18K]  mark_duplicates_with_picard.SRR4036062_ZTF_S.log
│   ├── [4.0K]  merge_count_tables.Siefert2017.log
│   ├── [   0]  picard_metrics/
│   │   ├── [3.4K]  SRR4036047_28hpf_G1.marked.duplicates_metrics.txt
│   │   ├── [3.4K]  SRR4036048_28hpf_S.marked.duplicates_metrics.txt
│   │   ├── [3.4K]  SRR4036049_28hpf_G1.marked.duplicates_metrics.txt
│   │   ├── [3.4K]  SRR4036050_28hpf_S.marked.duplicates_metrics.txt
│   │   ├── [3.5K]  SRR4036051_28hpf_G1.marked.duplicates_metrics.txt
│   │   ├── [3.5K]  SRR4036052_28hpf_S.marked.duplicates_metrics.txt
│   │   ├── [3.4K]  SRR4036059_Bud.marked.duplicates_metrics.txt
│   │   ├── [3.4K]  SRR4036060_Bud.marked.duplicates_metrics.txt
│   │   ├── [3.4K]  SRR4036061_ZTF_G1.marked.duplicates_metrics.txt
│   │   └── [3.4K]  SRR4036062_ZTF_S.marked.duplicates_metrics.txt
│   ├── [7.2K]  process_count_tables.Siefert2017.log
│   ├── [1.1K]  quality_filter_with_bamtools.SRR4036047_28hpf_G1.log
│   ├── [1.1K]  quality_filter_with_bamtools.SRR4036048_28hpf_S.log
│   ├── [1.1K]  quality_filter_with_bamtools.SRR4036049_28hpf_G1.log
│   ├── [1.1K]  quality_filter_with_bamtools.SRR4036050_28hpf_S.log
│   ├── [1.1K]  quality_filter_with_bamtools.SRR4036051_28hpf_G1.log
│   ├── [1.1K]  quality_filter_with_bamtools.SRR4036052_28hpf_S.log
│   ├── [1.1K]  quality_filter_with_bamtools.SRR4036059_Bud.log
│   ├── [1.1K]  quality_filter_with_bamtools.SRR4036060_Bud.log
│   ├── [1.1K]  quality_filter_with_bamtools.SRR4036061_ZTF_G1.log
│   ├── [1.1K]  quality_filter_with_bamtools.SRR4036062_ZTF_S.log
│   ├── [   0]  snakelogs/
│   │   ├── [   0]  Siefert2017_make_bedgraphs.log
│   │   ├── [   0]  Siefert2017_merge_count_tables.log
│   │   ├── [   0]  Siefert2017_process_count_tables.log
│   │   ├── [   0]  align_reads_with_bwamem.SRR4036047_28hpf_G1.log
│   │   ├── [   0]  align_reads_with_bwamem.SRR4036048_28hpf_S.log
│   │   ├── [   0]  align_reads_with_bwamem.SRR4036049_28hpf_G1.log
│   │   ├── [   0]  align_reads_with_bwamem.SRR4036050_28hpf_S.log
│   │   ├── [   0]  align_reads_with_bwamem.SRR4036051_28hpf_G1.log
│   │   ├── [   0]  align_reads_with_bwamem.SRR4036052_28hpf_S.log
│   │   ├── [   0]  align_reads_with_bwamem.SRR4036059_Bud.log
│   │   ├── [   0]  align_reads_with_bwamem.SRR4036060_Bud.log
│   │   ├── [   0]  align_reads_with_bwamem.SRR4036061_ZTF_G1.log
│   │   ├── [   0]  align_reads_with_bwamem.SRR4036062_ZTF_S.log
│   │   ├── [   0]  blacklist_filter_with_bedtools.SRR4036047_28hpf_G1.log
│   │   ├── [   0]  blacklist_filter_with_bedtools.SRR4036048_28hpf_S.log
│   │   ├── [   0]  blacklist_filter_with_bedtools.SRR4036049_28hpf_G1.log
│   │   ├── [   0]  blacklist_filter_with_bedtools.SRR4036050_28hpf_S.log
│   │   ├── [   0]  blacklist_filter_with_bedtools.SRR4036051_28hpf_G1.log
│   │   ├── [   0]  blacklist_filter_with_bedtools.SRR4036052_28hpf_S.log
│   │   ├── [   0]  blacklist_filter_with_bedtools.SRR4036059_Bud.log
│   │   ├── [   0]  blacklist_filter_with_bedtools.SRR4036060_Bud.log
│   │   ├── [   0]  blacklist_filter_with_bedtools.SRR4036061_ZTF_G1.log
│   │   ├── [   0]  blacklist_filter_with_bedtools.SRR4036062_ZTF_S.log
│   │   ├── [   0]  count_reads_in_RT_windows.SRR4036047_28hpf_G1.log
│   │   ├── [   0]  count_reads_in_RT_windows.SRR4036048_28hpf_S.log
│   │   ├── [   0]  count_reads_in_RT_windows.SRR4036049_28hpf_G1.log
│   │   ├── [   0]  count_reads_in_RT_windows.SRR4036050_28hpf_S.log
│   │   ├── [   0]  count_reads_in_RT_windows.SRR4036051_28hpf_G1.log
│   │   ├── [   0]  count_reads_in_RT_windows.SRR4036052_28hpf_S.log
│   │   ├── [   0]  count_reads_in_RT_windows.SRR4036059_Bud.log
│   │   ├── [   0]  count_reads_in_RT_windows.SRR4036060_Bud.log
│   │   ├── [   0]  count_reads_in_RT_windows.SRR4036061_ZTF_G1.log
│   │   ├── [   0]  count_reads_in_RT_windows.SRR4036062_ZTF_S.log
│   │   ├── [   0]  mark_duplicates_with_picard.SRR4036047_28hpf_G1.log
│   │   ├── [   0]  mark_duplicates_with_picard.SRR4036048_28hpf_S.log
│   │   ├── [   0]  mark_duplicates_with_picard.SRR4036049_28hpf_G1.log
│   │   ├── [   0]  mark_duplicates_with_picard.SRR4036050_28hpf_S.log
│   │   ├── [   0]  mark_duplicates_with_picard.SRR4036051_28hpf_G1.log
│   │   ├── [   0]  mark_duplicates_with_picard.SRR4036052_28hpf_S.log
│   │   ├── [   0]  mark_duplicates_with_picard.SRR4036059_Bud.log
│   │   ├── [   0]  mark_duplicates_with_picard.SRR4036060_Bud.log
│   │   ├── [   0]  mark_duplicates_with_picard.SRR4036061_ZTF_G1.log
│   │   ├── [   0]  mark_duplicates_with_picard.SRR4036062_ZTF_S.log
│   │   ├── [   0]  quality_filter_with_bamtools.SRR4036047_28hpf_G1.log
│   │   ├── [   0]  quality_filter_with_bamtools.SRR4036048_28hpf_S.log
│   │   ├── [   0]  quality_filter_with_bamtools.SRR4036049_28hpf_G1.log
│   │   ├── [   0]  quality_filter_with_bamtools.SRR4036050_28hpf_S.log
│   │   ├── [   0]  quality_filter_with_bamtools.SRR4036051_28hpf_G1.log
│   │   ├── [   0]  quality_filter_with_bamtools.SRR4036052_28hpf_S.log
│   │   ├── [   0]  quality_filter_with_bamtools.SRR4036059_Bud.log
│   │   ├── [   0]  quality_filter_with_bamtools.SRR4036060_Bud.log
│   │   ├── [   0]  quality_filter_with_bamtools.SRR4036061_ZTF_G1.log
│   │   ├── [   0]  quality_filter_with_bamtools.SRR4036062_ZTF_S.log
│   │   ├── [   0]  trim_reads_with_fastp.SRR4036047_28hpf_G1.log
│   │   ├── [   0]  trim_reads_with_fastp.SRR4036048_28hpf_S.log
│   │   ├── [   0]  trim_reads_with_fastp.SRR4036049_28hpf_G1.log
│   │   ├── [   0]  trim_reads_with_fastp.SRR4036050_28hpf_S.log
│   │   ├── [   0]  trim_reads_with_fastp.SRR4036051_28hpf_G1.log
│   │   ├── [   0]  trim_reads_with_fastp.SRR4036052_28hpf_S.log
│   │   ├── [   0]  trim_reads_with_fastp.SRR4036059_Bud.log
│   │   ├── [   0]  trim_reads_with_fastp.SRR4036060_Bud.log
│   │   ├── [   0]  trim_reads_with_fastp.SRR4036061_ZTF_G1.log
│   │   └── [   0]  trim_reads_with_fastp.SRR4036062_ZTF_S.log
│   ├── [2.5K]  trim_reads_with_fastp.SRR4036047_28hpf_G1.log
│   ├── [2.5K]  trim_reads_with_fastp.SRR4036048_28hpf_S.log
│   ├── [2.5K]  trim_reads_with_fastp.SRR4036049_28hpf_G1.log
│   ├── [2.5K]  trim_reads_with_fastp.SRR4036050_28hpf_S.log
│   ├── [2.5K]  trim_reads_with_fastp.SRR4036051_28hpf_G1.log
│   ├── [2.5K]  trim_reads_with_fastp.SRR4036052_28hpf_S.log
│   ├── [2.4K]  trim_reads_with_fastp.SRR4036059_Bud.log
│   ├── [2.4K]  trim_reads_with_fastp.SRR4036060_Bud.log
│   ├── [2.5K]  trim_reads_with_fastp.SRR4036061_ZTF_G1.log
│   └── [2.5K]  trim_reads_with_fastp.SRR4036062_ZTF_S.log
├── [   0]  merged/
│   └── [4.8M]  Siefert2017_counts.txt
├── [   0]  processed_rse/
│   └── [ 32M]  Siefert2017_processed_rse.rds
├── [   0]  qc/
│   └── [   0]  fastp_reports/
│       ├── [452K]  SRR4036047_28hpf_G1.html
│       ├── [117K]  SRR4036047_28hpf_G1.json
│       ├── [452K]  SRR4036048_28hpf_S.html
│       ├── [117K]  SRR4036048_28hpf_S.json
│       ├── [452K]  SRR4036049_28hpf_G1.html
│       ├── [117K]  SRR4036049_28hpf_G1.json
│       ├── [453K]  SRR4036050_28hpf_S.html
│       ├── [117K]  SRR4036050_28hpf_S.json
│       ├── [454K]  SRR4036051_28hpf_G1.html
│       ├── [118K]  SRR4036051_28hpf_G1.json
│       ├── [454K]  SRR4036052_28hpf_S.html
│       ├── [119K]  SRR4036052_28hpf_S.json
│       ├── [451K]  SRR4036059_Bud.html
│       ├── [115K]  SRR4036059_Bud.json
│       ├── [451K]  SRR4036060_Bud.html
│       ├── [116K]  SRR4036060_Bud.json
│       ├── [450K]  SRR4036061_ZTF_G1.html
│       ├── [116K]  SRR4036061_ZTF_G1.json
│       ├── [451K]  SRR4036062_ZTF_S.html
│       └── [116K]  SRR4036062_ZTF_S.json
├── [   0]  rse/
│   └── [2.0M]  Siefert2017_rse.rds
└── [   0]  trimmed/
    ├── [1.4G]  SRR4036047_28hpf_G1_trimmed_R1.fastq.gz
    ├── [1.4G]  SRR4036047_28hpf_G1_trimmed_R2.fastq.gz
    ├── [1.6G]  SRR4036048_28hpf_S_trimmed_R1.fastq.gz
    ├── [1.6G]  SRR4036048_28hpf_S_trimmed_R2.fastq.gz
    ├── [1.7G]  SRR4036049_28hpf_G1_trimmed_R1.fastq.gz
    ├── [1.8G]  SRR4036049_28hpf_G1_trimmed_R2.fastq.gz
    ├── [2.2G]  SRR4036050_28hpf_S_trimmed_R1.fastq.gz
    ├── [2.2G]  SRR4036050_28hpf_S_trimmed_R2.fastq.gz
    ├── [4.4G]  SRR4036051_28hpf_G1_trimmed_R1.fastq.gz
    ├── [4.4G]  SRR4036051_28hpf_G1_trimmed_R2.fastq.gz
    ├── [5.0G]  SRR4036052_28hpf_S_trimmed_R1.fastq.gz
    ├── [5.0G]  SRR4036052_28hpf_S_trimmed_R2.fastq.gz
    ├── [720M]  SRR4036059_Bud_trimmed_R1.fastq.gz
    ├── [720M]  SRR4036059_Bud_trimmed_R2.fastq.gz
    ├── [910M]  SRR4036060_Bud_trimmed_R1.fastq.gz
    ├── [907M]  SRR4036060_Bud_trimmed_R2.fastq.gz
    ├── [2.2G]  SRR4036061_ZTF_G1_trimmed_R1.fastq.gz
    ├── [2.3G]  SRR4036061_ZTF_G1_trimmed_R2.fastq.gz
    ├── [2.2G]  SRR4036062_ZTF_S_trimmed_R1.fastq.gz
    └── [2.2G]  SRR4036062_ZTF_S_trimmed_R2.fastq.gz

15 directories, 283 files
```
