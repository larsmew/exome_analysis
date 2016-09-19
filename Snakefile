ref = "../bwa-0.7.13/reference/human_g1k_v37_decoy.fasta"

'''
Map reads to reference genome with bwa (adding read groups with -R, random placeholder text)
Use samtools to sort mappings and convert to bam file
'''
rule bwa_map:
	input:
		ref={ref},
		faq=expand("{{sample}}_{pair}_001.fastq.gz", pair=["R1","R2"])
	output:
		"{sample}_addedGroup.bam"
	threads: 48
	shell:
		"bwa mem -t {threads} -R '@RG\tID:bwa\tLB:lib1\tPL:illumina\tSM:20\tPU:unit1' {input} | samtools sort - > {output}"


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

'''
rule get_dbsnp_vcf_file:
	params:
		link="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/"
	output:
		"dbsnp_138.b37.vcf"
	shell:
		"wget {params}{output}.gz | gunzip {output}.gz"

rule get_gold_indels_vcf_file:
	params:
		link="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/"
	output:
		"Mills_and_1000G_gold_standard.indels.b37.vcf"
	shell:
		"wget {params}{output}.gz | gunzip {output}.gz"

rule patterns_covariation_pass1:
	input:
		bam="{sample}_reDup_addedGroup.bam",
		gold_indels="Mills_and_1000G_gold_standard.indels.b37.vcf",
		dbsnp="dbsnp_138.b37.vcf"
	output:
		"{sample}_recal_data.table"
	shell:
		"GenomeAnalysisTK -T BaseRecalibrator -R {ref} -I {input.bam} -L 20 -knownSites {input.dbsnp} -knownSites {input.gold_indels} -o {output}"

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

rule apply_recalibration:
	input:
		bam="{sample}_reDup_addedGroup.bam",
		table="{sample}_recal_data.table"
	output:
		"{sample}_recal_reads.b"
	shell:
		"GenomeAnalysisTK -T PrintReads -R {ref} -I {input.bam} -L 20 -BQSR {input.table} -o {output}"





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

'''