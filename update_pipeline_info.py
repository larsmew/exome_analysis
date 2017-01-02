__author__ = "Lars Andersen <larsmew@gmail.com>"
__version__ = "1.0"
__date__ = "02/01/2017"

from subprocess import check_call
import sys, os, shutil, time, filecmp

# Files to be uploaded, if updated.
archive_updated = False
snakefile_updated = False
conda_updated = False


###############################################################################
####                                                                       ####
####                       Current pipeline version                        ####
####                                                                       ####
###############################################################################
'''
Extracts the information from the current pipeline version
'''

# Describe pipeline_version infor file and archive folder
pipeline_file = "pipeline_version.txt"
archive = "pipeline_archive/"

# Dirty hack to extract dates, as date lines are similar.
snakefile_flag = False
conda_flag = False

# Scan file and extract relevant information
with open(pipeline_file, "r") as f:
	for line in f:
		# Extract pipeline version
		if line.startswith("Pipeline"):
			pipeline_version = line.split(":")[1].strip()
			print("current pipeline version:", pipeline_version)
		# Extract Snakefile version
		elif line.startswith("Snakefile"):
			snakefile_version = line.split(":")[1].strip()
			print("current Snakefile version:", snakefile_version)
			snakefile_flag = True
		# Extract Snakefile last update date
		elif snakefile_flag:
			snakefile_date = line.split(":")[1].strip()
			print("current Snakefile date:", snakefile_date)
			snakefile_flag = False
		# To be updated manually - see 'File versions' section
		'''
		elif line.startswith("File"):
			file_version = line.split(":")[1].strip()
			print("current File version:", file_version)
		'''
		# Extract Conda package version
		if line.startswith("Conda Packages"):
			conda_version = line.split(":")[1].strip()
			print("current Conda version:", conda_version)
			conda_flag = True
		# Extract Conda last update date
		elif conda_flag:
			conda_date = line.split(":")[1].strip()
			print("current Conda date:", conda_date)
			conda_flag = False

###############################################################################
####                                                                       ####
####                               Snakefile                               ####
####                                                                       ####
###############################################################################
'''
This section checks if the Snakefile is updated since last pipeline update.
'''

# Extracts date of newest Snakefile
def get_snakefile_date(snakefile):
	with open(snakefile, "r") as f:
		for line in f:
			if line.startswith("__date__"):
				date = line.split('"')[1]
				# print(date)
				return date
			else:
				continue

# Extracts version of newest Snakefile
def get_snakefile_version(snakefile):
	with open(snakefile, "r") as f:
		for line in f:
			if line.startswith("__version__"):
				version = line.split('"')[1]
				# print(version)
				return version
			else:
				continue

# Method to check if Snakefile updated since last pipeline version
def is_snakefile_updated(snakefile_version, new_snakefile_version):
	# print(snakefile_version, new_snakefile_version)
	if snakefile_version == new_snakefile_version:
		return False
	else:
		return True

# Name of Snakefile
snakefile = "Snakefile"

# Check for snakefile updates
new_snakefile_version = get_snakefile_version(snakefile)
if is_snakefile_updated(snakefile_version, new_snakefile_version):
	# If Snakefile updated
	new_snakefile_date = get_snakefile_date(snakefile)
	snakefile_updated = True
else:
	# If Snakefile not updated
	new_snakefile_version = snakefile_version
	new_snakefile_date = snakefile_date

###############################################################################
####                                                                       ####
####                             File versions                             ####
####                                                                       ####
###############################################################################
'''
This section describes the files used for the pipeline.

Note: These are updated manually in this file at the moment
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
These functions initiates conda update system
'''

# Define file names
current_conda_file = "conda_packages_v"+conda_version+".txt"
new_conda_version = str(int(conda_version) + 1)
new_conda_file = "conda_packages_v"+new_conda_version+".txt"

# Method to check if conda packages updated since last pipeline version
def update_conda():
	# Update all conda packages
	try:
	    check_call(["conda", "update", "--all", "-y"])
	# Throw error if fails
	except subprocess.CalledProcessError:
	    # handle errors in the called executable
		print("ERROR: in 'conda update' command")
		sys.exit(1)
	except OSError:
		# executable not found
		print("ERROR: 'conda update' command not found")
		sys.exit(1)

	# Export updated packages to new conda file
	with open(new_conda_file, "w") as f:
		# Run export command
		try:
		    check_call(["conda", "list", "--export"], stdout=f)
		# Throw error if fails
		except subprocess.CalledProcessError:
		    # handle errors in the called executable
			print("ERROR: in 'conda list' command")
			sys.exit(1)
		except OSError:
			# executable not found
			print("ERROR: 'conda list' command not found")
			sys.exit(1)

# Checks if any conda packages have been updated
def is_conda_updated(current_conda_file, new_conda_file):
	print("compares:", current_conda_file, "and", new_conda_file)
	if filecmp.cmp(current_conda_file, new_conda_file):
		os.remove(new_conda_file)
		print("No conda packages updated")
		return False
	else:
		print("conda packages updated")
		return True

# Run functions for conda packages
update_conda()

if is_conda_updated(current_conda_file, new_conda_file):
	new_conda_date = time.strftime("%d/%m/%Y")
	conda_updated = True
else:
	new_conda_version = conda_version
	new_conda_date = conda_date

###############################################################################
####                                                                       ####
####                            Update pipeline                            ####
####                                                                       ####
###############################################################################
'''
Creates the new pipeline_version file and commits and pushes it to git
'''
# Checks if pipeline updated
def is_pipeline_updated(pipeline_version, new_pipeline_version):
	if pipeline_version == new_pipeline_version:
		return False
	else:
		return True

# The information to write in new pipeline_version file
def write_pipeline_version_content(new_pipeline_version):
	with open("pipeline_version.txt", "w") as f:
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

		f.write("Conda Packages - version: " + new_conda_version + "\n")
		f.write("Updated: " + new_conda_date + "\n\n")

		with open("conda_packages_v"+new_conda_version+".txt", "r") as cp:
			for line in cp:
				f.write(line)

# Commits to git and pushes updates to GitHub
def commit_and_push_to_github(new_pipeline_version):
	# Add achive to git if updated
	if archive_updated:
		check_call(["git", "add", archive])
		print("Added achive to commit")

	# Obs: Snakefile is added to github when script updates
	# if snakefile_updated:
	# 	os.system("git add Snakefile")

	# Add achive to git if updated
	if conda_updated:
		check_call(["git", "add", new_conda_file])
		check_call(["git", "rm", current_conda_file])
		print("Added new conda_packages file and removed old one")

	# Add pipeline file to git
	check_call(["git", "add", pipeline_file])
	print("Added new pipeline_version file")

	# Commit changes to git
	check_call(["git", "commit", "-m",
			   "'Updated pipeline version to "+new_pipeline_version+"'"])
	print("Commited pipeline version:", new_pipeline_version)

	# Push updates to GitHub
	s = check_call(["git", "push"])
	if s == 0:
		print("Pushed changes to remote git (GitHub)")
	else:
		print("Failed to push changes to remote git (GitHub)")

# Copy current pipeline file to archive, if not already there
def archive_pipeline(pipe_vers):
	# Graps current pipeline version file
	current_pipeline_file = pipeline_file.split(".")[0]+"_"+pipe_vers+".txt"
	print(current_pipeline_file)

	# Check if already in archive, if not, then archive it
	if not os.path.isfile(archive+current_pipeline_file):
		# Uses copy2 to preserve metadata
		shutil.copy2(pipeline_file, archive+current_pipeline_file)
		# Flag update
		print("Archived", archive+current_pipeline_file)
		archive_updated = True

# Copy current conda packages file to archive, if not already there
def archive_conda(conda_file):
	print(conda_file)
	# Check if already in archive, if not, then archive it
	if not os.path.isfile(archive+conda_file):
		# Uses copy2 to preserve metadata
		shutil.copy2(conda_file, archive+conda_file)
		# Flag update
		print("Archived", archive+conda_file)
		archive_updated = True

# Method to archive, update and commit pipeline
def archive_and_update_pipeline():

	# Archive old versions
	archive_pipeline(pipeline_version)
	archive_conda(current_conda_file)

	# Get the newest pipeline version
	new_pipeline_version = new_snakefile_version + "." + file_version \
						 + "." + new_conda_version

	# If pipeline is updated
	if is_pipeline_updated(pipeline_version, new_pipeline_version):

		# Create the new pipeline_version file
		write_pipeline_version_content(new_pipeline_version)
		print("Pipeline updated to version:", new_pipeline_version)

		# Archive new versions
		archive_pipeline(new_pipeline_version)
		if conda_updated:
			archive_conda(new_conda_file)

			# Remove old version from working dir
			os.remove(current_conda_file)
			print("Removed old conda file from top folder")

		# Commit changes to git and push to github
		commit_and_push_to_github(new_pipeline_version)
	else:
		# If no updates to pipeline
		print("Pipeline already up-to-date")

archive_and_update_pipeline()

# Method to fix a failed commit, only for bugfixing
def fix_failed_commit():
	new_pipeline_version = "1.0.1.3"
	new_conda_file = "conda_packages_v3.txt"
	current_conda_file = "conda_packages_v2.txt"
	archive_updated = True
	conda_updated = True
	commit_and_push_to_github(new_pipeline_version)

# Fixing failed commit (uncomment to run)
# fix_failed_commit()