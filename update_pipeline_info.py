__author__ = "Lars Andersen <larsmew@gmail.com>"
__version__ = "1.0"
__date__ = "02/01/2017"

from subprocess import check_call
import sys, os, shutil, time, filecmp, subprocess

# Files to be uploaded, if updated.
archive_updated = False
snakefile_updated = False
conda_updated = False
py2_updated = False
py3_updated = False
r_updated = False
condaenv_py2_updated = False


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
python2_flag = False
python3_flag = False
r_flag = False
condaenv_ver_flag = False
condaenv_date_flag = False

# list of conda environments
conda_env = []


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
		# Extract python2 version
		if line.startswith("Python2 packages"):
			py2_version = line.split(":")[1].strip()
			print("current Python2 version:", py2_version)
			python2_flag = True
		# Extract python2 date
		elif python2_flag:
			py2_date = line.split(":")[1].strip()
			print("current Python2 date:", py2_date)
			python2_flag = False
		# Extract python3 version
		if line.startswith("Python3 packages"):
			py3_version = line.split(":")[1].strip()
			print("current Python3 version:", py3_version)
			python3_flag = True
		# Extract python3 date
		elif python3_flag:
			py3_date = line.split(":")[1].strip()
			print("current Python3 date:", py3_date)
			python3_flag = False
		# Extract R version
		if line.startswith("R packages"):
			r_version = line.split(":")[1].strip()
			print("current R version:", r_version)
			r_flag = True
		# Extract python3 date
		elif r_flag:
			r_date = line.split(":")[1].strip()
			print("current R date:", r_date)
			r_flag = False
		# Extract conda environment name
		elif line.startswith("Conda environment:"):
			condaenv_name = line.split(":")[1].strip()
			print("Conda environment:", condaenv_name)
			condaenv_ver_flag = True
		# Extract conda environment version
		elif condaenv_ver_flag:
			condaenv_ver = line.split(":")[1].strip()
			print("Conda environment version:", condaenv_ver)
			condaenv_ver_flag = False
			condaenv_date_flag = True
		# Extract conda environment date
		elif condaenv_date_flag:
			condaenv_date = line.split(":")[1].strip()
			print("Conda environment updated:", condaenv_date)
			condaenv_date_flag = False
			conda_env.append((condaenv_name, condaenv_ver, condaenv_date))


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

Note: These are updated manually in this file as they changes rarely
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
new_conda_date = time.strftime("%d/%m/%Y")

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
	
	# Checks if any conda packages updated or new installed 
	if is_conda_updated(current_conda_file, new_conda_file):
		# set updated date and flags
		conda_updated = True
	else:
		# No updates, use current version
		new_conda_version = conda_version
	

# Checks if any conda packages have been updated
def is_conda_updated(current_conda_file, new_conda_file):
	print("compares:", current_conda_file, "and", new_conda_file)
	if filecmp.cmp(current_conda_file, new_conda_file):
		os.remove(new_conda_file) # remove temp file
		print("No conda packages updated")
		return False
	else:
		print("conda packages updated")
		return True

# Run functions for conda packages
# update_conda()

###############################################################################
####                                                                       ####
####                       Conda environements                             ####
####                                                                       ####
###############################################################################
'''
Packages in conda environments other than default are updated manually 
as there are older packages installed for some tools to work properly
e.g. bedtools 2.23.00
'''

new_conda_env = []
new_condaenv_py2_file = ""
current_condaenv_py2_file = ""

def update_conda_env(conda_env):
	for env, ver, date in conda_env:
		print(env, ver, date)
		
		current_conda_env_file = "conda_env_"+env+"_v"+ver+".txt"
		new_ver = str(int(ver) + 1)
		new_conda_env_file = "conda_env_"+env+"_v"+new_ver+".txt"

		os.system("source activate "+env+"; conda list --export > "+new_conda_env_file)
		
		# Export updated packages to new conda file
		# with open(new_conda_env_file, "w") as f:
		# 	# Run export command
		# 	try:
		# 	    check_call(["source", "activate", env+";", "conda", "list", "--export"], stdout=f)
		# 	# Throw error if fails
		# 	except subprocess.CalledProcessError:
		# 	    # handle errors in the called executable
		# 		print("ERROR: in 'conda list' command")
		# 		sys.exit(1)
		# 	except OSError:
		# 		# executable not found
		# 		print("ERROR: 'conda list' command not found")
		# 		sys.exit(1)
		
		# Checks if any conda packages updated or new installed 
		if is_conda_updated(current_conda_env_file, new_conda_env_file):
			# set updated date and flags
			new_date = time.strftime("%d/%m/%Y")
			new_env = (env, new_ver, new_date)
			new_conda_env.append(new_env)
			if env == "py2":
				new_condaenv_py2_file = new_conda_env_file
				current_condaenv_py2_file = current_conda_env_file
				condaenv_py2_updated = True
		else:
			# No updates, use current version
			new_date = time.strftime("%d/%m/%Y")
			old_env = (env, ver, new_date)
			new_conda_env.append(old_env)


def is_conda_env_updated(current_conda_env_file, new_conda_env_file):
	print("compares:", current_conda_env_file, "and", new_conda_env_file)
	if filecmp.cmp(current_conda_env_file, new_conda_env_file):
		os.remove(new_conda_env_file) # remove temp file
		print("Conda environment not updated")
		return False
	else:
		print("Conda environment updated")
		return True

#update_conda_env(conda_env)

###############################################################################
####                                                                       ####
####                           python packages                             ####
####                                                                       ####
###############################################################################
'''
These functions initiates the python update system
'''

current_py2_file = "conda_env_py2_v"+py2_version+".txt"
new_py2_version = str(int(py2_version) + 1)
new_py2_file = "conda_env_py2_v"+new_py2_version+".txt"
new_py2_date = time.strftime("%d/%m/%Y")

def update_python2():
	
	# Update python2 packages
	os.system("source activate py2; \
			   pip freeze --local | grep -v '^\-e' | cut -d = -f 1  | xargs -n1 pip install -U; \
			   pip freeze > " + new_py2_file)
			   
	# Checks if any python2 packages updated or new installed
	print("compares:", current_py2_file, "and", new_py2_file)
	if filecmp.cmp(current_py2_file, new_py2_file):
		print("Python2 packages not updated") # Files the same
		os.remove(new_py2_file) # remove temp file
		new_py2_version = py2_version
	else:
		print("Python2 packages are updated")
		py2_updated = True


current_py3_file = "conda_env_py3_v"+py3_version+".txt"
new_py3_version = str(int(py3_version) + 1)
new_py3_file = "conda_env_py3_v"+new_py3_version+".txt"
new_py3_date = time.strftime("%d/%m/%Y")

def update_python3():
	
	# Update python3 packages
	os.system("pip freeze --local | grep -v '^\-e' | cut -d = -f 1  | xargs -n1 pip install -U; \
			   pip freeze > " + new_py3_file)
			   
	# Checks if any python3 packages updated or new installed
	print("compares:", current_py3_file, "and", new_py3_file)
	if filecmp.cmp(current_py3_file, new_py3_file):
		print("Python3 packages already up-to-date") # Files the same
		os.remove(new_py3_file) # remove temp file
		new_py3_version = py3_version
	else:
		print("Python3 packages are updated")
		py3_updated = True
		


#update_python2()
#update_python3()

###############################################################################
####                                                                       ####
####                               R packages                              ####
####                                                                       ####
###############################################################################
'''
These functions initiates the R package update / backup system


'''

# from rpy2 import *
#
# installedpkgs = ro.r('installed.packages()[,c("Package", "Version")]')
# ro.r('write.table(installedpkgs, file=filename, row.names = F, col.names = T, sep = "=", quote = F)')

current_r_file = "R_packages_v"+r_version+".txt"
new_r_version = str(int(r_version) + 1)
new_r_name = "R_packages_v"+new_r_version
new_r_file = new_r_name+".txt"
new_r_date = time.strftime("%d/%m/%Y")

def update_R():
	os.system("Rscript updateR.R " + new_r_name)
	
	if filecmp.cmp(current_r_file, new_r_file):
		print("R packages already up-to-date") # Files the same
		os.remove(new_r_file) # remove temp file
		os.remove(new_r_name+".rda")
		new_r_version = r_version
	else:
		print("R packages are updated")
		r_updated = True
	
#update_R()

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
		
		f.write("Python2 packages - version:" + new_py2_version + "\n")
		f.write("Updated:" + new_py2_date + "\n\n")

		f.write("Python3 packages - version:" + new_py3_version + "\n")
		f.write("Updated:" + new_py3_date + "\n\n")

		f.write("R packages - version:" + new_r_version + "\n")
		f.write("Updated:" + new_R_date + "\n\n")
		
		f.write("# Conda environments (manually updated) \n")
		for env, ver, date in new_conda_env:
			f.write("Conda environment: " + env + "\n")
			f.write("Version:" + ver + "\n")
			f.write("Updated:" + date + "\n")
		
		# NO LONGER PRINTS ALL PACKAGES AT END - TOO MUCH (WASTED) SPACE
		# with open("conda_packages_v"+new_conda_version+".txt", "r") as cp:
		# 	for line in cp:
		# 		f.write(line)

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
		archive_file(new_conda_file)
		check_call(["git", "rm", current_conda_file])
		os.remove(current_conda_file)
		print("Added new conda_packages file and removed old one")
	
	# Add achive to git if updated
	if py2_updated:
		check_call(["git", "add", new_py2_file])
		archive_file(new_py2_file)
		check_call(["git", "rm", current_py2_file])
		os.remove(current_py2_file)
		print("Added new py2_packages file and removed old one")

	# Add achive to git if updated
	if py3_updated:
		check_call(["git", "add", new_py3_file])
		archive_file(new_py3_file)
		check_call(["git", "rm", current_py3_file])
		os.remove(current_py3_file)
		print("Added new py3_packages file and removed old one")

	# Add achive to git if updated
	if r_updated:
		check_call(["git", "add", new_r_file])
		archive_file(new_r_file)
		check_call(["git", "rm", current_r_file])
		os.remove(current_r_file)
		print("Added new r_packages file and removed old one")
	
	# Add achive to git if updated
	if condaenv_py2_updated:
		check_call(["git", "add", new_condaenv_py2_file])
		archive_file(new_condaenv_py2_file)
		check_call(["git", "rm", current_condaenv_py2_file])
		os.remove(current_condaenv_py2_file)
		print("Added new condaenv_py2_packages file and removed old one")

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
		print("Changes already commited!")
		print("Check internet connection or if GitHub is down")

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
def archive_file(file):
	print(file)
	# Check if already in archive, if not, then archive it
	if not os.path.isfile(archive+file):
		# Uses copy2 to preserve metadata
		shutil.copy2(file, archive+file)
		# Flag update
		print("Archived", archive+file)
		archive_updated = True

# Method to archive, update and commit pipeline
def archive_and_update_pipeline():

	# Archive old versions
	archive_pipeline(pipeline_version)
	archive_file(current_conda_file)
	archive_file(current_condaenv_py2_file)
	archive_file(current_py2_file)
	archive_file(current_py3_file)
	archive_file(current_r_file)

	# Get the newest pipeline version
	new_pipeline_version = new_snakefile_version + "." + file_version \
						 + "." + new_conda_version + "." + new_py2_version \
						 + "." + new_py3_version + "." + new_r_version

	# If pipeline is updated
	if is_pipeline_updated(pipeline_version, new_pipeline_version):

		# Create the new pipeline_version file
		write_pipeline_version_content(new_pipeline_version)
		print("Pipeline updated to version:", new_pipeline_version)

		# Archive new pipeline version
		archive_pipeline(new_pipeline_version)

		# Commit changes to git and push to github
		commit_and_push_to_github(new_pipeline_version)
	else:
		# If no updates to pipeline
		print("Pipeline already up-to-date")

#archive_and_update_pipeline()

###############################################################################
####                                                                       ####
####                             Run commands                              ####
####                                                                       ####
###############################################################################

update_conda()
update_conda_env(conda_env)
update_python2()
update_python3()
update_R()
archive_and_update_pipeline()

###############################################################################
####                                                                       ####
####                               Bugfixing                               ####
####                                                                       ####
###############################################################################
'''
Methods for fixing bugs and the state of pipeline files.
'''
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
