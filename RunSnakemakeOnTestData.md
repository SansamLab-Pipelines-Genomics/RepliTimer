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
![Config File Image](https://github.com/SansamLab/Process_RepTiming_Snakemake/blob/make-test-vignette/resources/configFileImage.png)


#### 4A. Modify the config/config.yml file
The config.yml file should be preconfigured for the test data set. Check it to be sure:

### 4. Run pipeline with test data
```
sbatch --wrap="\
snakemake \
-R \
-j 999 \
--cluster-config config/cluster_config.yml \
--cluster '\
sbatch \
--constraint=westmere \
-A {cluster.account} \
-p {cluster.partition} \
--cpus-per-task {cluster.cpus-per-task}  \
--mem {cluster.mem} \
--output {cluster.output}'"
```
