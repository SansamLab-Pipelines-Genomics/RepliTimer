### 1.  Load slurm and miniconda
Note. The commands to do this will be different on your machine. These commands are specific to an HPC using slurm with these modules installed.

```bash
ml slurm
ml miniconda
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
### 3A.  FIRST TIME ONLY:  Setup conda environment
```bash
# -f is the location of the environment .yml file. 
## The relative path assumes that you are in the root directory of this repository.
# -p is the path where you want to install this environment
conda env create -f workflow/envs/RTEnv.yml -p /s/sansam-lab/RT_Environment 
```

### 3B.  Activate conda environment
```bash
conda activate /s/sansam-lab/RT_Environment
```

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