__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "05/12/2016"
__version__ = "1.0"

import sys, os, shutil, time

###############################################################################
####                                                                       ####
####                       Current pipeline version                        ####
####                                                                       ####
###############################################################################
pipeline_file = "pipeline_version.txt"
archive = "pipeline_archive/"

# Dirty hack to extract dates, as date lines are similar.
snakefile_flag = False
conda_flag = False

# Scan file and extract relevant information
with open(pipeline_file, "r") as f:
	for line in f:
		if line.startswith("Pipeline"):
			pipeline_version = line.split(":")[1].strip()
			print("current pipeline version:", pipeline_version)
		elif line.startswith("Snakefile"):
			snakefile_version = line.split(":")[1].strip()
			print("current Snakefile version:", snakefile_version)
			snakefile_flag = True
		elif snakefile_flag:
			snakefile_date = line.split(":")[1].strip()
			print("current Snakefile date:", snakefile_date)
			snakefile_flag = False
		# To be updated manually - see 'File versions' section
		# elif line.startswith("File"):
		# 	file_version = line.split(":")[1].strip()
		# 	print("current File version:", file_version)
		elif line.startswith("Conda Packages"):
			conda_version = line.split(":")[1].strip()
			print("current Conda version:", conda_version)
			conda_flag = True
		elif conda_flag:
			conda_date = line.split(":")[1].strip()
			print("current Conda date:", conda_date)
			conda_flag = False

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
				# print(date)
				return date
			else:
				continue

def get_snakefile_version(snakefile):
	with open(snakefile, "r") as f:
		for line in f:
			if line.startswith("__version__"):
				version = line.split('"')[1]
				# print(version)
				return version
			else:
				continue

def is_snakefile_updated(snakefile_version, new_snakefile_version):
	# print(snakefile_version, new_snakefile_version)
	if snakefile_version == new_snakefile_version:
		return False
	else:
		return True

snakefile = "Snakefile"

# Check for snakefile updates
if is_snakefile_updated(snakefile_version, get_snakefile_version(snakefile)):
	new_snakefile_version = get_snakefile_version(snakefile)
	new_snakefile_date = get_snakefile_date(snakefile)
else:
	new_snakefile_version = snakefile_version
	new_snakefile_date = snakefile_date

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

# Remember to update this, if any of the above file changes!
file_version = '1'
file_date = "06/10/2016"


###############################################################################
####                                                                       ####
####                            Conda packages                             ####
####                                                                       ####
###############################################################################
'''
TBD
'''
def is_conda_updated():
	return False

if is_conda_updated():
	new_conda_version = str(int(conda_version) + 1)
	new_conda_date = time.strftime("%d/%m/%Y")
else:
	new_conda_version = conda_version
	new_conda_date = conda_date

###############################################################################
####                                                                       ####
####                            Update pipeline                            ####
####                                                                       ####
###############################################################################
'''

'''
def is_pipeline_updated(pipeline_version, new_pipeline_version):
	if pipeline_version == new_pipeline_version:
		return False
	else:
		return True

def write_pipeline_version_content(new_pipeline_version):
	with open("test_file.txt", "w") as f:
		f.write("Pipeline version: " + new_pipeline_version + "\n")
		f.write("Date: " + time.strftime("%d/%m/%Y") + "\n\n")

		f.write("Snakefile version: " + new_snakefile_version + "\n")
		f.write("Updated: " + new_snakefile_date + "\n\n")

		f.write("File version: " + file_version + "\n")
		f.write("Updated: " + file_date + "\n\n")

		f.write("Reference: " + Reference + "\n")
		f.write("dbsnp: " + dbsnp + "\n")
		f.write("Mills_and_1000G: " + Mills_and_1000G + "\n")
		f.write("interval: " + interval + "\n")
		f.write("bed: " + bed + "\n\n")

		f.write("Conda packages - version: " + new_conda_version + "\n")
		f.write("Updated: " + new_conda_date + "\n\n")

		with open("conda_packages_v"+new_conda_version+".txt", "r") as cp:
			for line in cp:
				f.write(line)

def archive_and_update_pipeline():
	
	# Copy current pipeline file to archive, if not already there
	current_pipeline_file = pipeline_file.split(".")[0]+"_"+pipeline_version+".txt"
	print(current_pipeline_file)
	if not os.path.isfile(archive+current_pipeline_file):
		shutil.copy2(pipeline_file, archive+current_pipeline_file)
		print("Archived", archive+current_pipeline_file)
	
	# Copy current conda packages file to archive
	current_conda_file = "conda_packages_v"+conda_version+".txt"
	print(current_conda_file)
	if not os.path.isfile(archive+current_conda_file):
		shutil.copy2(current_conda_file, archive+current_conda_file)
		print("Archived", archive+current_conda_file)
	
	new_pipeline_version = new_snakefile_version + "." + file_version \
						 + "." + new_conda_version
	
	if is_pipeline_updated:
		write_pipeline_version_content(new_pipeline_version)
		print("Pipeline updated to version:", new_pipeline_version)
	else:
		print("Pipeline is already up-to-date")

archive_and_update_pipeline()
