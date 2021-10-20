__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "10/01/2020"
__version__ = "1.0"

import time, os

configfile: samples.yaml

SAMPLES, = glob_wildcards("{sample}_R1.fastq.gz")
# FAMNAME = os.getcwd().rsplit("/",1)[1]

totim = time.time()
timeFormat = "%Y_%m_%d:%X" # year, month, day, time H:M:S
log_file = "log_gatk4-wes-germline-hofte-dysplasi_" + timeFormat + ".txt"

mem = "-Xmx12g" # login nodes
# mem = "-Xmx24g" # slim nodes
# mem = "-Xmx32g" # Fat nodes
# mem = "-Xmx50g"

# Resources
ref = "/work/sduvarcall/bwa-0.7.13/reference/human_g1k_v37_decoy.fasta"
bed = ""
interval = ""
dbsnp = "/work/sduvarcall/knownSNPs/dbsnp_150.b37.vcf.gz"
mills_1000G = "/work/sduvarcall/knownSNPs/Mills_and_1000G_gold_standard.indels.b37.vcf"
cosmic = "/work/sduvarcall/cosmic/Cosmic-combined_v81_b37.vcf"
hapmap = "/work/sduvarcall/knownSNPs/hapmap_3.3.b37.vcf"
omni = "/work/sduvarcall/knownSNPs/1000G_omni2.5.b37.vcf"
phase1_1000G = "/work/sduvarcall/knownSNPs/1000G_phase1.indels.b37.vcf"


onstart:
	shell("echo 'Started execution of pipeline:' $(date +'%Y-%m-%d %H:%M:%S') >> {log_file}")
	shell("mkdir -p Metrics; mkdir -p bam; mkdir -p gvcf_files")

onsuccess:
	fiTime = 'Total time:', str((time.time()-totim) / 60), 'minutes'
	shell("echo 'Finished execution on' $(date +'%Y-%m-%d %H:%M:%S') >> {log_file}")
	shell("echo {fiTime} >> {log_file}")

onerror:
	fiTime = 'Total time:', str((time.time()-totim) / 60), 'minutes'
	shell("echo 'Finished execution on' $(date +'%Y-%m-%d %H:%M:%S') >> {log_file}")
	shell("echo {fiTime} >> {log_file}")
	shell("echo 'ERROR OCCURED, PLEASE REFER TO SLURM LOGFILE FOR DETAILS' >> {log_file}")


###############################################################################
### Rule all                                                                ###
###############################################################################
rule all:
	input:
		# BAM file
		expand("bam/{sample}_recal.bam", sample=SAMPLES),
		
		# VCF files
        expand("{fam_name}_variants_recal_gatk.vcf.gz", fam_name=FAMNAME),
		#expand("{fam_name}_variants_recal_gatk.vcf.gz", fam_name=FAMNAME),
		#expand("{fam_name}_variants_combined_freebayes.vcf", fam_name=FAMNAME),
		
		# Metrics and Statistics
        expand("Metrics/{sample}_post_recalibration.grp", sample=SAMPLES),
		expand("Metrics/{sample}.HS_Metrics.txt", sample=SAMPLES),
		expand("Metrics/{sample}_unsorted.quality_distribution_metrics", sample=SAMPLES),
		expand("Metrics/{sample}.alignment_summary_metrics", sample=SAMPLES),
		expand("Metrics/{sample}_aggregation.alignment_summary_metrics", sample=SAMPLES)



###############################################################################
### Create recalibrated bam file                                            ###
###############################################################################
'''
Map reads to reference genome with bwa
'''
rule MapAndSort:
    input:
		f1 = "{sampleid}_{protocol}_{flowcell}_R1.fastq.gz",
		f2 = "{sampleid}_{protocol}_{flowcell}_R2.fastq.gz"
    output:
        bam = temp("{sampleid}_{protocol}_{flowcell}_sorted.bam")
	params:
    	rgid = "{sampleid}_{protocol}_{flowcell}",
    	rglb = "{protocol}",
    	rgsm = "{sampleid}",
    	rgpl = "Illumina",
    	rgpu = "{flowcell}",
    threads: 24
    shell:
        """
        bwa mem -M -t {threads} \
        -R \"@RG\\tID:{params.rgid}\\tLB:{params.rglb}\\tSM:{params.rgsm}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\" \
        {ref} {input} | \
        gatk --java-options {mem} SortSam \
        --INPUT=/dev/stdin \
        --OUTPUT={output.bam} \
        --VALIDATION_STRINGENCY=LENIENT \
        --SORT_ORDER=coordinate \
        --CREATE_INDEX=TRUE
        """


'''
Remove duplicate reads
'''
rule MarkDuplicates:
	input:
		bam = "{sample}_sorted.bam",
		bai = "{sample}_sorted.bai"
	output:
		bam = temp("{sample}_dedup.bam"),
		bai = temp("{sample}_dedup.bai"),
		met = "Metrics/{sample}_duplicate_metrics.txt"
	shell:
		"""
		gatk --java-options {mem} MarkDuplicates \
		--INPUT={input.bam} \
		--OUTPUT={output.bam} \
		--METRICS_FILE={output.met} \
		--OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
		--CREATE_INDEX=true \
		"""


'''
Obtain recalibration information
'''
rule BaseRecalibrator:
	input:
		bam = "{sample}_dedup.bam",
		bai = "{sample}_dedup.bai",
	output:
		table = "Metrics/{sample}_pre_recalibration.grp"
	shell:
		"""
		gatk --java-options {mem} BaseRecalibrator \
		-R={ref} \
		-I={input.bam} \
		--known-sites={dbsnp} \
		--known-sites={mills_1000G} \
        --known-sites={phase1_1000G} \
		-O={output} \
		"""


'''
Apply Recalibration
'''
rule ApplyRecalibration:
	input:
		bam="{sample}_dedup.bam",
		bai="{sample}_dedup.bai",
		table = "Metrics/{sample}_pre_recalibration.grp"
	output:
		"bam/{sample}_recal.bam"
	shell:
		"""
		gatk --java-options {mem} ApplyBQSR \
		-R={ref} \
		-I={input.bam} \
		--bqsr-recal-file={input.table} \
		--add-output-sam-program-record \
		-O={output} \
		"""

###############################################################################
### Collect Statistics                                                      ###
###############################################################################
'''
Collect recalibration metrics
'''
rule RecalibrationMetrics:
	input:
		bam = "bam/{sample}_recal.bam",
	output:
		table = "Metrics/{sample}_post_recalibration.grp"
	shell:
		"""
		gatk --java-options {mem} BaseRecalibrator \
		-R={ref} \
		-I={input.bam} \
		--known-sites={dbsnp} \
		--known-sites={mills_1000G} \
        --known-sites={phase1_1000G} \
		-O={output} \
		"""

'''
Collect Hybrid Selection (HS) metrics
'''
rule collectHsMetrics:
	input:
		bam="bam/{sample}_recal.bam"
	output:
		"Metrics/{sample}.HS_Metrics.txt"
	shell:
		"""
		gatk --java-options {mem} CollectHsMetrics \
		--INPUT={input.bam} \
		--REFERENCE_SEQUENCE={ref} \
		--OUTPUT={output} \
		--BAIT_INTERVALS={interval} \
		--TARGET_INTERVALS={interval} \
		"""



rule CollectUnsortedReadgroupBamQualityMetrics:
	input:
		bam="bam/{sample}_recal.bam"
	output:
		"Metrics/{sample}_unsorted.quality_distribution_metrics"
		# And many more, not written here for simplicity
	params:
		out="Metrics/{sample}_unsorted"
	shell:
		"""
		gatk --java-options {mem} CollectMultipleMetrics \
		--INPUT={input.bam} \
		--OUTPUT={params} \
		--ASSUME_SORTED=true \
		--PROGRAM="null" \
		--PROGRAM="CollectBaseDistributionByCycle" \
		--PROGRAM="CollectInsertSizeMetrics" \
		--PROGRAM="MeanQualityByCycle" \
		--PROGRAM="QualityScoreDistribution" \
		--METRIC_ACCUMULATION_LEVEL="null" \
		--METRIC_ACCUMULATION_LEVEL="ALL_READS" \
		"""


rule CollectReadgroupBamQualityMetrics:
	input:
		bam="bam/{sample}_recal.bam"
	output:
		"Metrics/{sample}.alignment_summary_metrics"
		# And many more, not written here for simplicity
	params:
		out="Metrics/{sample}"
	shell:
		"""
		gatk --java-options {mem} CollectMultipleMetrics \
		--INPUT={input.bam} \
		--REFERENCE_SEQUENCE={ref} \
		--OUTPUT={params} \
		--ASSUME_SORTED=true \
		--PROGRAM="null" \
		--PROGRAM="CollectAlignmentSummaryMetrics" \
		--PROGRAM="CollectGcBiasMetrics" \
		--METRIC_ACCUMULATION_LEVEL="null" \
		--METRIC_ACCUMULATION_LEVEL="READ_GROUP" \
		"""
	
rule CollectAggregationMetrics:
	input:
		bam="bam/{sample}_recal.bam"
	output:
		"Metrics/{sample}_aggregation.alignment_summary_metrics"
		# And many more, not written here for simplicity
	params:
		out="Metrics/{sample}_aggregation"
	shell:
		"""
		gatk --java-options {mem} CollectMultipleMetrics \
		--INPUT={input.bam} \
		--REFERENCE_SEQUENCE={ref} \
		--OUTPUT={params} \
		--ASSUME_SORTED=true \
		--PROGRAM="null" \
		--PROGRAM="CollectAlignmentSummaryMetrics" \
		--PROGRAM="CollectInsertSizeMetrics" \
		--PROGRAM="CollectGcBiasMetrics" \
		--PROGRAM="CollectSequencingArtifactMetrics" \
		--PROGRAM="QualityScoreDistribution" \
		--METRIC_ACCUMULATION_LEVEL="null" \
		--METRIC_ACCUMULATION_LEVEL="SAMPLE" \
		--METRIC_ACCUMULATION_LEVEL="LIBRARY"
		"""
        

###############################################################################
### Call germline variants                                                  ###
###############################################################################
'''
HaplotypeCaller
'''
rule HaplotypeCaller:
	input:
		bam="bam/{sample}_recal.bam"
	output:
		gvcf="gvcf_files/{sample}_variants.g.vcf.gz"
	shell:
		"""
		gatk --java-options {mem} HaplotypeCaller \
		-R={ref} \
		-I={input.bam} \
		-O={output.gvcf} \
		-ERC=GVCF \
		--dbsnp={dbsnp} \
		-A Coverage \
		-A DepthPerAlleleBySample \
		-A BaseQualityRankSumTest
		"""
        # Coverage - Total depth of coverage per sample and over all samples
        # DepthPerAlleleBySample - Unfiltered alternative allele depth (AD)
        # BaseQualityRankSumTest - Rank Sum Test of REF versus ALT base quality scores

'''
Create single multi-sample g.vcf file    
'''
rule CombineGVCFs:
	input:
		gvcfs=temp(expand("gvcf_files/{sample}_variants.g.vcf.gz", sample=config["{fam_name}"]))
	output:
		gvcf="gvcf_files/{fam_name}_variants_combined.g.vcf.gz"
	params:
		gvcfs=expand("-V=gvcf_files/{sample}_variants.g.vcf.gz", sample=config["{fam_name}"]),
	shell:
		"""
		gatk --java-options {mem} CombineGVCFs \
		-R {ref} \
		{params.gvcfs} \
		-O {output.gvcf}
		"""

'''
GenotypeGVCFs
'''
rule GenotypeGVCFs:
	input:
		gvcf="gvcf_files/{fam_name}_variants_combined.g.vcf.gz"
	output:
		vcf="{fam_name}_variants.vcf.gz"
	shell:
		"""
		gatk --java-options {mem} GenotypeGVCFs \
		-R={ref} \
		-V={input.gvcf} \
		--dbsnp={dbsnp} \
		-O={output} \
		"""


# '''
# Recalibrate Variant Quality Score
# '''
# rule VariantRecalibrator:
#     input:
#         vcf="{fam_name}_variants.vcf.gz"
#     output:
#         recal="{fam_name}_variants.recal",
#         tranches="{fam_name}_variants.tranches",
#         rscript="{fam_name}_variants.plots.R"
#     shell:
#         """
#         gatk --java-options {mem} VariantRecalibrator \
#         -R {ref} \
#         -V {input} \
#         -O {output.recal} \
#         --resource hapmap,known=false,training=true,truth=true,prior=15.0:{hapmap} \
#         --resource omni,known=false,training=true,truth=false,prior=12.0:{omni} \
#         --resource 1000G,known=false,training=true,truth=false,prior=10.0:{phase1_1000G} \
#         --resource dbsnp,known=true,training=false,truth=false,prior=2.0:{dbsnp} \
#         -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
#         -mode BOTH \
#         --tranches-file {output.tranches} \
#         --rscript-file {output.rscript}
#         """
#
#
# '''
# Apply Variant Recalibration
# '''
# rule ApplyVQSR:
#     input:
#         vcf="{fam_name}_variants.vcf.gz",
#         recal="{fam_name}_variants.recal",
#         tranches="{fam_name}_variants.tranches"
#     output:
#         vcf="{fam_name}_variants_recal_gatk.vcf.gz"
#     shell:
#         """
#         gatk --java-options {mem} ApplyVQSR \
#         -R {ref} \
#         -V {input.vcf} \
#         -O {output} \
#         --tranches-file {input.tranches} \
#         --recal-file {input.recal} \
#         --mode BOTH
#         """
