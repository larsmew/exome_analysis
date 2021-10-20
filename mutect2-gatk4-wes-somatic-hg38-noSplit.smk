__title__ = "Pipeline for Somatic Variant Calling with Mutect2 - Kristina's project"
__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "10/09/2019"
__version__ = "1.0"

import time, os, sys

#########################################################
####                       Input                     ####
#########################################################
# Matched tumour-normal samples information
configfile: "somatic_matched_samples.yaml"

# Panel of normals location
pon_location = "somaticPanelOfNormals/"

# Sample information
# FAMNAME = glob_wildcards("{sample}.bam")[0]
# list_fam = FAMNAME
# print(list_fam)
# FAMNAME = set(['-'.join(fam.split("-",2)[0:2]) for fam in glob_wildcards("data/{sample}.bam")[0]])
FAMNAME, = glob_wildcards("data/{sample}_tumor_tagseq-medexome.recalibrated.bam")
FAMNAME = [fam+"_tumor_tagseq-medexome" for fam in FAMNAME]
print(FAMNAME)

NORMALS, = glob_wildcards(pon_location+"{normal}.bam")
# NORMALS, = glob_wildcards("{normal}_normal_*.bam")
print(NORMALS)
# sys.exit()

# Explicit paths for external input files
ref = "/work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta"
# gnomead = "/work/sduvarcall/knownSNPs/gnomead/af-only-gnomad.raw.sites.b37.vcf.gz"
gnomad = "/work/sdukoldby/resources/hg38/af-only-gnomad.hg38.vcf.gz"
interval_list = "/work/sdukoldby/resources/hg38/MedExome_hg38_capture_targets.interval_list"
common_variants = "/work/sdukoldby/resources/hg38/small_exac_common_3.hg38.vcf.gz"


#########################################################
####                      Output                     ####
#########################################################
log_file = "log_file_somatic.txt"
output_somatic = "somatic_variants"


#########################################################
####                       Setup                     ####
#########################################################
# Timing
totim = time.time()
timeFormat = "%Y_%m_%d:%X" # year, month, day, time H:M:S

# Memory
# mem = "-Xmx12g" # login nodes - be careful not running too many jobs at once!
# mem = "-Xmx24g" # slim nodes
# mem = "-Xmx32g" # Fat nodes
mem = "-Xmx64g"

# Define chromosomes - perhaps implement function to extract from reference
# CHROM = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
#          "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT",
# 		 "hs37d5"]
#
# PATCHES = ["GL000207.1", "GL000226.1", "GL000229.1", "GL000231.1", "GL000210.1",
# 		   "GL000239.1", "GL000235.1", "GL000201.1", "GL000247.1", "GL000245.1",
# 		   "GL000197.1", "GL000203.1", "GL000246.1", "GL000249.1", "GL000196.1",
# 		   "GL000248.1", "GL000244.1", "GL000238.1", "GL000202.1", "GL000234.1",
# 		   "GL000232.1", "GL000206.1", "GL000240.1", "GL000236.1", "GL000241.1",
# 		   "GL000243.1", "GL000242.1", "GL000230.1", "GL000237.1", "GL000233.1",
# 		   "GL000204.1", "GL000198.1", "GL000208.1", "GL000191.1", "GL000227.1",
# 		   "GL000228.1", "GL000214.1", "GL000221.1", "GL000209.1", "GL000218.1",
# 		   "GL000220.1", "GL000213.1", "GL000211.1", "GL000199.1", "GL000217.1",
# 		   "GL000216.1", "GL000215.1", "GL000205.1", "GL000219.1", "GL000224.1",
# 		   "GL000223.1", "GL000195.1", "GL000212.1", "GL000222.1", "GL000200.1",
# 		   "GL000193.1", "GL000194.1", "GL000225.1", "GL000192.1", "NC_007605"]

# CHROM = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]

# CHROM = ["chr19", "chr20"]

# PATCHES = [l.strip('\n') for l in open("contigs_hg38.txt", "r")]
# print(PATCHES[0:10])
# print(len(PATCHES))

# sys.exit()


#########################################################
####  Define workflow start, stop and error actions  ####
#########################################################
onstart:
	# shell("echo $(head /work/sduvarcall/exome_analysis/pipeline_version.txt -n 1) >> {log_file}")
	# shell("echo 'Started execution of pipeline:' $(date +'%Y-%m-%d %H:%M:%S') >> {log_file}")
	shell("mkdir -p "+output_somatic)

# onsuccess:
# 	fiTime = 'Total time:', str((time.time()-totim) / 60), 'minutes'
# 	shell("echo 'Finished execution on' $(date +'%Y-%m-%d %H:%M:%S') >> {log_file}")
# 	shell("echo Execution time: {fiTime} >> {log_file}")
# 	shell("cat {log_file} >> log_file_somatic_success.txt")
#
# onerror:
# 	fiTime = 'Total time:', str((time.time()-totim) / 60), 'minutes'
# 	shell("echo 'Finished execution on' $(date +'%Y-%m-%d %H:%M:%S') >> {log_file}")
# 	shell("echo Execution time: {fiTime} >> {log_file}")
# 	shell("echo 'ERROR OCCURED, PLEASE REFER TO SLURM LOGFILE FOR DETAILS' >> {log_file}")
# 	shell("cp {log_file} log_file_somatic_error.txt")


#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''
# Matched normal-tumor samples
# if FAMNAME in config:
# 	PAIRS = config[FAMNAME]
# 	print(PAIRS)

print([expand(output_somatic+"/{tumor}_vs_{normal}_somatic_variants_mutect2_filtered.vcf",
			normal=config[fam]["normal"],
			tumor=config[fam]["tumor"]) for fam in FAMNAME])
# sys.exit()
	
rule all_pairs:
	input:
		[expand(output_somatic+"/{tumor}_vs_{normal}_somatic_variants_mutect2_filtered.vcf",
			normal=config[fam]["normal"],
			tumor=config[fam]["tumor"]) for fam in FAMNAME]


#########################################################
####       Create Somatic Panel of Normals           ####
#########################################################
'''
Run Mutect2 in tumor-only mode for each normal sample
'''
rule Mutect2_tumor_only_pon:
	input:
		bam=pon_location+"{sample}.bam"
	output:
		vcf=pon_location+"{sample}_for_pon.vcf",
		idx=pon_location+"{sample}_for_pon.vcf.idx"
	threads: 24
	shell:
		"""
		gatk --java-options {mem} Mutect2 \
		-R {ref} \
		-I {input} \
		-max-mnp-distance 0 \
		--native-pair-hmm-threads {threads} \
		-L {interval_list} \
		-O {output.vcf}
		"""
		# -tumor ECV2-35-020118-blod_tagseq-medexome_HC3NJBGX9_S5 \

'''
Create GenomicsDB from the normal Mutect2 calls
'''
# rule GenomicsDB:
# 	input:
# 		vcf=expand(pon_location+"{sample}_for_pon.vcf", sample=NORMALS)
# 	output:
# 		db=pon_location+"pon_db"
# 	params:
# 		vcf=expand("-V "+pon_location+"/{sample}_for_pon.vcf", sample=NORMALS)
# 	shell:
# 		"""
# 		gatk GenomicsDBImport \
# 		-R {ref} \
# 		-L {interval_list} \
# 		--genomicsdb-workspace-path {output} \
# 		{params.vcf}
# 		"""


'''
Combine the normal calls using CreateSomaticPanelOfNormals.
'''
rule CreateSomaticPanelOfNormals:
	input:
		# db=pon_location+"pon_db"
		vcf=expand(pon_location+"{sample}_for_pon.vcf", sample=NORMALS)
	output:
		pon=pon_location+"pon.vcf.gz"
	params:
		vcf=expand("-V "+pon_location+"/{sample}_for_pon.vcf", sample=NORMALS)
	shell:
		"""
		gatk --java-options {mem} CreateSomaticPanelOfNormals \
		-R {ref} \
		{params.vcf} \
		-O {output}
		"""
		# -V gendb://{input} \

##########################################################
####  Call Somatic Variants using Mutect2 on matched  ####
####   Tumor-Normal samples on per chromesome basis   ####
##########################################################
'''
Mutect2 on matched Tumor-Normal samples
'''
rule Mutect2_matched:
	input:
		normal="data/{normal}.recalibrated.bam",
		tumor="data/{tumor}.recalibrated.bam",
		pon=pon_location+"pon.vcf.gz"
	output:
		vcf=output_somatic+"/{tumor}_vs_{normal}_somatic_variants_mutect2.vcf",
		idx=output_somatic+"/{tumor}_vs_{normal}_somatic_variants_mutect2.vcf.idx",
		f1r2="{tumor}_vs_{normal}_f1r2.tar.gz"
	threads: 24
	shell:
		"""
		gatk --java-options {mem} Mutect2 \
		-R {ref} \
		-I {input.tumor} \
		-I {input.normal} \
		-pon {input.pon} \
		--germline-resource {gnomad} \
		--af-of-alleles-not-in-resource 0.0000025 \
		--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
		--native-pair-hmm-threads {threads} \
		--f1r2-tar-gz {output.f1r2} \
		-L {interval_list} \
		-O {output.vcf} \
		-bamout {wildcards.tumor}_mutect2.bam
		"""
		# -tumor {wildcards.tumor} \
		# -normal {wildcards.normal} \


#########################################################
####          Learn Read Orientation Bias            ####
#########################################################
rule learnReadOrientationModel:
	input:
		"{tumor}_vs_{normal}_f1r2.tar.gz"
	output:
		"{tumor}_vs_{normal}_read-orientation-model.tar.gz"
	shell:
		"""
		gatk LearnReadOrientationModel \
		-I {input} \
		-O {output}
		"""

#########################################################
####           Create Contamination table            ####
#########################################################
### Sæt sammen til en metode, som tager højde for med og uden normal... (Se gatk4 scripts)
rule GetPileupSummaries_normal:
	input:
		bam="data/{normal}.recalibrated.bam"
	output:
		pileup=output_somatic+"/{normal}_normal_pileup.table"
	shell:
		"""
		gatk --java-options {mem} GetPileupSummaries \
		-I {input.bam} \
		-V {common_variants} \
		-L {common_variants} \
		-O {output}
		"""

rule GetPileupSummaries_tumor:
	input:
		bam="data/{tumor}.recalibrated.bam"
	output:
		pileup=output_somatic+"/{tumor}_tumor_pileup.table"
	shell:
		"""
		gatk --java-options {mem} GetPileupSummaries \
		-I {input.bam} \
		-V {common_variants} \
		-L {common_variants} \
		-O {output}
		"""


rule CalculateContamination:
	input:
		normal=output_somatic+"/{normal}_normal_pileup.table",
		tumor=output_somatic+"/{tumor}_tumor_pileup.table"
	output:
		contamination=output_somatic+"/{tumor}_vs_{normal}_contamination.table"
	shell:
		"""
		gatk --java-options {mem} CalculateContamination \
		-I {input.tumor} \
		-matched {input.normal} \
		-O {output}
		"""


#########################################################
####              Filter Mutect2 Calls               ####
#########################################################
rule FilterMutectCalls:
	input:
		vcf=output_somatic+"/{tumor}_vs_{normal}_somatic_variants_mutect2.vcf",
		idx=output_somatic+"/{tumor}_vs_{normal}_somatic_variants_mutect2.vcf.idx",
		contamination=output_somatic+"/{tumor}_vs_{normal}_contamination.table",
		read_orientation="{tumor}_vs_{normal}_read-orientation-model.tar.gz"
	output:
		vcf=output_somatic+"/{tumor}_vs_{normal}_somatic_variants_mutect2_filtered.vcf",
		idx=output_somatic+"/{tumor}_vs_{normal}_somatic_variants_mutect2_filtered.vcf.idx"
	shell:
		"""
		gatk --java-options {mem} FilterMutectCalls \
		-V {input.vcf} \
		--contamination-table {input.contamination} \
		--orientation-bias-artifact-priors {input.read_orientation} \
		-L {interval_list} \
		-O {output.vcf}
		"""