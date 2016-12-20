# Conda for Bioinformatic Pipelines

This document describes some simple commands to get the conda package up and running on 64-bit Linux machines (and likely MacOS machines).

## Install Conda package manager

There exist two versions of the conda package manager: **Miniconda** and **Anaconda**.

The difference is that Anaconda includes 150+ open source packages. If you are mostly doing bioinformatic pipelines I suggest you install **Miniconda** as we have to install most tools anyway.

### Install Miniconda

Run the following two commands to install conda. Follow the on-screen installation instructions.
The default location for conda is `~/miniconda3`.

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Remember to prepend the miniconda path to PATHS. Either, during install or by running the following:

Bash: `export PATH=~/miniconda3/bin:$PATH`

Fish: `set -U fish_user_paths ~/miniconda3/bin $fish_user_paths`

**Restart the shell for the changes to take effect.**

### Install Anaconda

To be done...

## Search for packages to install

Conda searches the anaconda repository as default

```bash
conda search <package_name>
```

**Example:**
```bash
conda search python
```

### Search BioConda

Most of the tools we are interested in will be located in the bioconda repository. To search this, we use:
```bash
conda search -c bioconda <package_name>
```

**Example:**
```bash
conda search -c bioconda samtools
```

## Install packages

In general:

```bash
conda install <package_name>
```

**Example:**
```bash
conda install python
```

### Install from BioConda

Installing from the bioconda repository is as easy as searching it:
```bash
conda install -c bioconda <package_name>
```

**Example:**
```bash
conda install -c bioconda samtools
conda install -c bioconda bwa
conda install -c bioconda picard
conda install -c bioconda gatk
conda install -c bioconda snakemake
```

## Backup and Restore conda environments

For reprodibility it can sometimes be necessary rerun analysis. Therefore it can be crucial to run the analysis exactly as it was done the first time with same set of packages and versions.

### Backup

To backup the current conda environment to a file use:

```bash
conda list --export > <file_name.txt>
```

**Example:**
```bash
conda list --export > conda_packages.txt
```

### Restore

To restore a conda environment without touching the current installation use the following command which creates a new conda environment with the specific packages defined from the backup file:

```bash
conda create --name <env> --file <file_name.txt>
```

Now, this environment can be activated:

```bash
source activate <env>
```

or deactivated:

```bash
source deactivate
```

**Example:**

```bash
conda create --name old_env --file conda_packages.txt
source activate old_env
```
