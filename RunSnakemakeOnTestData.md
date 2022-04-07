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

### 5A. Run pipeline with installed modules (Alternative 2)
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
