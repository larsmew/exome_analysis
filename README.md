# exome_analysis
A repository for analyzing exome data using snakemake (makefile) to combine all the steps of the analysis

## Snakefile
The Snakefile currently describes a pipeline for analyzing exomes

**(To be updated with steps)**

## Update_pipeline_info.py

A script to automaize pipeline versions and updates for easy and systematic documentation of programs, files and tools used for the analysis.

Version 1 includes the following functionality:

- Update pipeline_version.txt describing the pipeline files, Snakefile version and conda packages.
- Extracting the pipeline information from a pipeline_version file.
- Update conda packages
- Backup conda packages list
- Archiving pipeline_version files and conda_packages files in pipeline_archive folder for later restoration of a specific pipeline version

Functionality planned for version 2:

- Easy restore of a specific pipeline version into seperate folder and conda environment for a conflict free restore without touching the current setup.
