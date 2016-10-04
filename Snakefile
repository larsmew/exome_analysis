import time

ref = "../bwa-0.7.13/reference/human_g1k_v37_decoy.fasta"
timeFormat = "%Y_%m_%d:%X" # year, month, day, time H:M:S

'''
Create log file, containing:
	- Programs used and their version info
	- Start and stop time of complete workflow
	- Start and stop time of complete workflow
	- sample(s)
'''
rule log_file:
	output:
		log="log_file.txt"
	shell:
		"touch {input} \ "
		"&& echo $(date +'%Y-%m-%d:%H:%M:%S')"
		

'''
Map reads to reference genome with bwa (adding read groups with -R, random placeholder text)
Uses samtools to sort mappings and convert to bam file
'''
rule bwa_map:
	input:
		ref={ref},
		faq=expand("{{sample}}_{pair}_001.fastq.gz", pair=["R1","R2"])
	output:
		"{sample}_addedGroup.bam"
	threads: 48
	shell:
		"bwa mem \ -t {threads} -R '@RG\tID:bwa\tLB:lib1\tPL:illumina\tSM:20\tPU:unit1' \ "
		"{input} | samtools sort - > {output}"


'''
Remove duplicate mappings - for faster processing later
'''
rule remove_duplicates:
	input:
		"{sample}_addedGroup.bam"
	output:
		"{sample}_reDup_addedGroup.bam"
	shell:
		"picard -Xmx4g MarkDuplicates I={input} O={output} REMOVE_DUPLICATES=true METRICS_FILE=duplicates.txt"


'''
Creates index bai file with Samtools
'''
rule index_bam:
	input:
		"{sample}_reDup_addedGroup.bam"
	output:
		"{sample}_reDup_addedGroup.bam.bai"
	shell:
		"samtools index {input}"


'''
Downloads list of known SNPs if not already in the directory
'''
rule get_dbsnp_vcf_file:
	params:
		link="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/"
	output:
		"dbsnp_138.b37.vcf"
	shell:
		"wget {params}{output}.gz | gunzip {output}.gz"


'''
Downloads list of golden standard indels if not already in the directory
'''
rule get_gold_indels_vcf_file:
	params:
		link="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/"
	output:
		"Mills_and_1000G_gold_standard.indels.b37.vcf"
	shell:
		"wget {params}{output}.gz | gunzip {output}.gz"


'''
Obtain recalibration information
'''
rule patterns_covariation_pass:
	input:
		bam="{sample}_reDup_addedGroup.bam",
		gold_indels="Mills_and_1000G_gold_standard.indels.b37.vcf",
		dbsnp="dbsnp_138.b37.vcf"
	output:
		"{sample}_recal_data.table"
	shell:
		"GenomeAnalysisTK \ "
		"-T BaseRecalibrator \ "
		"-R {ref} \ "
		"-I {input.bam} \ "
		"-L 20 \ "
		"-knownSites {input.dbsnp} \ "
		"-knownSites {input.gold_indels} \ "
		"-cov ReadGroupCovariate \ "
		"-cov QualityScoreCovariate \ "
		# "-cov CycleCovariate \ "
		# "-cov ContextCovariate \ "
		# "-rf BadCigar \ "
		"-o {output}"


'''
Apply Recalibration
'''
rule apply_recalibration:
	input:
		bam="{sample}_reDup_addedGroup.bam",
		table="{sample}_recal_data.table"
	output:
		"{sample}_recal_reads.bam"
	shell:
		"GenomeAnalysisTK \ "
		"-T PrintReads \ "
		"-R {ref} \ "
		"-I {input.bam} \ "
		"-L 20 \ "
		"-BQSR {input.table} \ "
		# "-rf BadCigar \ "
		"-o {output}"


'''
Rule to create directory for VCF files		
'''
rule create_vcf_dir:
	output:
		"vcf_files"
	shell:
		"mkdir {output}"

'''
HaplotypeCaller
'''
rule haplotypeCaller:
	input:
		bam="{sample}_recal_reads.bam",
		# dbsnp="dbsnp_138.b37.vcf",
		vcfdir="vcf_files"
	output:
		gvcf="{input.vcfdir}/{sample}_raw_variants.g.vcf",
		#gvcf_list="{sample}_gvcf_files.txt"
	shell:
		"GenomeAnalysisTK \ "
		"-T HaplotypeCaller \ "
		"-R {ref} \ "
		"-I {input.bam} \ "
		"-L 20 \ "
		# "--dbsnp {input.dbsnp} \ "
		"--genotyping_mode DISCOVERY \ "
		"-stand_emit_conf 10 \ "
		"-stand_call_conf 30 \ "
		"--emitRefConfidence GVCF \ "
		"--variant_index_type LINEAR \ "
		"--variant_index_parameter 128000 \ "
		"-o {output} \ "
		# "| bgzip -c > {output} \ "
		# "&& tabix -p vcf {output}"
		"&& echo '{output}' >> gvcf_files.txt"


'''
GenotypeGVCFs
'''
rule genotypeGVCFs:
	input:
		gvcfs="gvcf_files.txt",
		dbsnp="dbsnp_138.b37.vcf",
		vcfdir="vcf_files"
	output:
		"{sample}_variants.vcf"
	shell:
		"GenomeAnalysisTK \ "
		"-T GenotypeGVCFs \ "
		"-R {ref} \ "
		"-V {input.gvcfs} \ "
		"--dbsnp {input.dbsnp} \ "
		"-L 20 \ "
		"-stand_emit_conf 10 \ "
		"-stand_call_conf 30 \ "
		"-o {input.vcfdir}/{output}"


'''

Maybe useful rules:

rule sort_bam:
	input:
		"{sample}_reDup.bam"
	output:
		"{sample}_reDup_sorted.bam"
	shell:
		"picard SortSam SO=coordinate I={input} O={output}"

rule add_read_groups:
	input:
		"{sample}_reDup_sorted.bam"
	output:
		"{sample}_reDup_sorted_addedGroups.bam"
	shell:
		"picard -Xmx4g AddOrReplaceReadGroups I={input} O={output} RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20"
		
		
rule realignerTargetCreator:
	input:
		ref={ref},
		bam="{sample}_reDup_addedGroup.bam",
		bai="{sample}_reDup_addedGroup.bam.bai"
	output:
		"{sample}_output.intervals"
	shell:
		"GenomeAnalysisTK -T RealignerTargetCreator -R {input.ref} -o {output} -I {input.bam}"

-----
Test af recalibration
		
rule patterns_covariation_pass2:
	input:
		bam="{sample}_reDup_addedGroup.bam",
		gold_indels="Mills_and_1000G_gold_standard.indels.b37.vcf",
		dbsnp="dbsnp_138.b37.vcf",
		table="{sample}_recal_data.table"
	output:
		"{sample}_post_recal_data.table"
	shell:
		"GenomeAnalysisTK -T BaseRecalibrator -R {ref} -I {input.bam} -L 20 -knownSites {input.dbsnp} -knownSites {input.gold_indels} -BQSR {input.table} -o {output}"

rule generate_recalibration_plots:
	input:
		before="{sample}_recal_data.table",
		after="{sample}_post_recal_data.table"
	output:
		"{sample}_recalibration_plots.pdf"
	shell:
		"GenomeAnalysisTK -T AnalyzeCovariates -R {ref} -L 20 -before {input.before} -after {input.after} -plots {output}"
-----

'''