__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "05/12/2016"
__version__ = "1.0"

import sys, os

updated_snakefile = False

###############################################################################
####                                                                       ####
####                            Current version                            ####
####                                                                       ####
###############################################################################
pipeline_file = "pipeline_version.txt"

with open(pipeline_file, "r") as f:
	for line in f:
		if line.startswith("Pipeline"):
			pipeline_version = line.split(":")[1].strip()
			print("current pipeline version:", pipeline_version)
		if line.startswith("Snakefile"):
			snakefile_version = line.split(":")[1].strip()
			print("current Snakefile version:", snakefile_version)
		if line.startswith("File"):
			file_version = line.split(":")[1].strip()
			print("current File version:", file_version)
		if line.startswith("Conda"):
			conda_version = line.split(":")[1].strip()
			print("current Conda version:", conda_version)

###############################################################################
####                                                                       ####
####                               Snakefile                               ####
####                                                                       ####
###############################################################################
def get_snakefile_date(snakefile):
	with open(snakefile, "r") as f:
		for line in f:
			if line.startswith("__date__"):
				date = line.split('"')[1]
				print(date)
				return date
			else:
				continue

def get_snakefile_version(snakefile):
	with open(snakefile, "r") as f:
		for line in f:
			if line.startswith("__version__"):
				version = line.split('"')[1]
				print(version)
				return float(version)
			else:
				continue

#def is_snakefile_updated(snakefile_version):
	

snakefile = "Snakefile"
new_snakefile_version = get_snakefile_version(snakefile)
new_snakefile_date = get_snakefile_date(snakefile)


###############################################################################
####                                                                       ####
####                             File versions                             ####
####                                                                       ####
###############################################################################
'''
These are updated manually at the moment
'''
Reference = "human_g1k_v37_decoy"
dbsnp = "dbsnp_138.b37.vcf"
Mills_and_1000G = "Mills_and_1000G_gold_standard.indels.b37.vcf"
interval = "NimbleGen_ExomeV3_UTR_CustomTargetRegion.interval"
bed = "NimbleGen_ExomeV3_UTR_CustomTargetRegion_padding100bp.bed"