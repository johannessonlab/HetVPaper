### SNPpop: Variant calling and population structure in *Podospora anserina*
#############################################################################
# Variant calling for all the sequenced Wageningen collection

### Sources:
# https://software.broadinstitute.org/gatk/documentation/quickstart.php
# First steps: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
# Germline short variant discovery (SNPs + Indels): https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145
# HaplotypeCaller: https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
# GenomicsDBImport: https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.1.0/org_broadinstitute_hellbender_tools_genomicsdb_GenomicsDBImport.php
# GenotypeGVCFs: https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.1.0/org_broadinstitute_hellbender_tools_walkers_GenotypeGVCFs.php
# GATK hard filtering: https://software.broadinstitute.org/gatk/documentation/article?id=23216#2
# VCF INFO annotations: https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_annotator_ReadPosRankSumTest.php
# http://corearray.sourceforge.net/tutorials/SNPRelate/

#############################################################################
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/04/10-06/07
# ---------------
# +++++++++++++++++++++++++++++++++++++++++++++++++
# Version 1


# -------------------------------------------------
# Samples
samples = config["SampleIDs"]

# Illumina reads path
Illumina = config["Illumina"]

# The reference genome
Podan2file = config["Podan2file"]

# Path to repeated elements library
TElib = config["TElib"]

# Plotting script
rscript = config["rscript"]

# Data for plotting script
metadata = config["metadata"]
matingdata = config["matingdata"]

# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: indexbwa, snpsvcf, snpsvcfnomiss, bedtoolsTEs, bgzip_tabix, VariantFiltration, getPASSsites, bgzip_tabix2, snpsvcfnomiss, makeintervals
# ----------

rule all:
	input:
		# Some numbers on the number of sites filtered
		"results/filteringstats.txt",
		# Figures
		"results/Figures/Fig2_MatingVsSNPs.pdf", # The main figure
		"results/Figures/FigS4_PCAs.pdf" # The supplementary figure
		
# ---------------------------------

rule indexbwa:
	""" Index genome with BWA """
	input:
		genome = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa"
	output:
		index = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa.bwt"
	version: "1"
	shell:
		"""
		bwa index {input.genome}
		"""

rule bwa_mem:
	""" Map Illumina reads with BWA """
	input:
		genome = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa",
		index = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa.bwt",
		read1 = "data/Illumina/{sample}_postQC.1.fq.gz",
		read2 = "data/Illumina/{sample}_postQC.2.fq.gz",
	output:
		bwaoutput = temp("mapping/{sample}/{sample}-to-Podan2.bam.sorted"),
		log = "logs/bwa_mem/{sample}.log"
	params:
		time = "3:30:00",
		threads = 10, 
		refbase = "Podan2",
		rg = "@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina",
	version: "1"
	shell:
		"""
		(bwa mem {input.genome} {input.read1} {input.read2} -t {params.threads} -R '{params.rg}' -M | samtools view -Su - | samtools sort -l 5 -O bam -T {wildcards.sample}'-to-'{params.refbase} -@ {params.threads} > {output.bwaoutput}) 2> {output.log}
		# -l 5 following Doug
		"""

rule markduplicates:
	""" Mark duplicates in BAM """
	# https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
	input:
		bwaoutput = "mapping/{sample}/{sample}-to-Podan2.bam.sorted"
	output:
		mdoutput = "mapping/{sample}/{sample}-to-Podan2.sorted.debup.bam",
		mdmetrics = "mapping/{sample}/{sample}-to-Podan2.sorted.metrics.txt"
	params:
		time = "1:30:00",
		threads = 3,
	version: "1"
	shell:
		"""
		# Using normal Picard
		picard MarkDuplicates I={input.bwaoutput} O={output.mdoutput} M={output.mdmetrics} ASSUME_SORT_ORDER=coordinate CREATE_INDEX=true TMP_DIR="temp"

		"""	
		# # VALIDATION_STRINGENCY=ValidationStringency
		# #                               Validation stringency for all SAM files read by this program.  Setting stringency to
		# #                               SILENT can improve performance when processing a BAM file in which variable-length data
		# #                               (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. This
		# #                               option can be set to 'null' to clear the default value. Possible values: STRICT,
		# #                               LENIENT, SILENT
		# # CREATE_INDEX=Boolean          Whether to create a BAM index when writing a coordinate-sorted BAM file.  Default value:
		# #                               false. This option can be set to 'null' to clear the default value. Possible values:
		# #                               true, false
		# # TMP_DIR (File)  Default value: null. This option may be specified 0 or more times.

# -------------- GATK4 -----------------
rule indexsanddict:
	""" Index reference for GATK """ 
	input:
		genome = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa",
	output:
		indexsamtools = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa.fai",
		diction = "data/Podan2/Podan2_AssemblyScaffoldsmt.dict"
	params:
		time = "1:00:00",
		threads = 1,
	version: "1.1"
	shell:
		"""
		# Make a reference index
		samtools faidx {input.genome}

		# Make a reference dictionary
		picard CreateSequenceDictionary R={input.genome} O={output.diction}
		"""	

rule HaplotypeCaller:
	""" Produce a GVCF file from BAM - haploid """
	# https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
	input:
		bam = "mapping/{sample}/{sample}-to-Podan2.sorted.debup.bam",
		ref = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa",
		indexsamtools = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa.fai",
		diction = "data/Podan2/Podan2_AssemblyScaffoldsmt.dict"
	output:
		gvcf = "gvcfs/{sample}.g.vcf",
		newbam = "realignedbams/{sample}.newhaplo.bam",
	params:
		time = "15:00:00",
		threads = 1,
		JavaMem = int(1 * 6.8),
		ploidy = 1
	shell:
		"""
		gatk --java-options "-Xmx{params.JavaMem}G" HaplotypeCaller \\
		-I {input.bam} -R {input.ref} \\
		-O {output.gvcf} \\
		-ploidy {params.ploidy} \\
		-ERC GVCF \\
		--bam-output {output.newbam}

		# -nct {params.threads} \\
		# --bam-output,-bamout:String   File to which assembled haplotypes should be written  Default value: null.
		# --annotateNDA 	Annotate number of alleles observed
		# --useNewAFCalculator	Use new AF model instead of the so-called exact model
		# --emitRefConfidence GVCF 	Mode for emitting reference confidence scores
		"""

rule makeintervals:
	""" Make an interval file for GenomicsDBImport """
	input:
		"data/Podan2/Podan2_AssemblyScaffoldsmt.fa.fai"
	output:
		"data/Podan2/Podan2.intervals"
	shell:
		"cat {input} | cut -f1 > {output}"

rule GenomicsDBImport:
	""" Import single-sample GVCFs into GenomicsDB before joint genotyping """
	# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.1.1/org_broadinstitute_hellbender_tools_genomicsdb_GenomicsDBImport.php
	input:
		"data/Podan2/Podan2.intervals",
		expand("gvcfs/{sample}.g.vcf", sample = samples)
	output:
		directory("genomicsdb")
	params:
		time = "2:00:00",
		threads = 1,
		JavaMem = int(8 * 6.8),
		ploidy = 1
	run:
		# Create a string in the format --variant path/to/gvcf/sample1 --variant path/to/gvcf/sample2 etc...
		variantlist = ""
		for sample in input[1:]: # Ignore the reference in the end
			variantlist += "--variant " + sample + " "
		
		gatkcommand = 'gatk --java-options "-Xmx{}G" GenomicsDBImport {} --genomicsdb-workspace-path {} -L {}'.format(params.JavaMem, variantlist, output[0], input[0]) 
		shell(gatkcommand)
		# Notice GATK will create the output directory

rule GenotypeGVCFs:
	""" Perform joint genotyping """
	# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.1.0/org_broadinstitute_hellbender_tools_walkers_GenotypeGVCFs.php
	input:
		my_database = "genomicsdb",
		ref = "data/Podan2/Podan2_AssemblyScaffoldsmt.fa",
	output:
		rawvcf = "results/PodoPop.vcf.gz"
	params:
		time = "2:00:00",
		threads = 1,
		JavaMem = int(8 * 6.8),
		ploidy = 1
	shell:
		"""
		gatk --java-options "-Xmx{params.JavaMem}G" GenotypeGVCFs \\
		-R {input.ref} \\
		-V gendb://{input.my_database} \\
		-O {output.rawvcf} \\
		-ploidy {params.ploidy} \\
		--create-output-variant-index 

		# --create-output-variant-index	If true, create a VCF index when writing a coordinate-sorted VCF file.
		# --use-new-qual-calculator / -new-qual 	Use the new AF model instead of the so-called exact model. Default: true
		# By default, GATK HaplotypeCaller and GenotypeGVCFs do not emit variants with QUAL < 10, controlled with -stand-call-conf
		"""

# ------- RepeatMasking -------

rule repeatmasker:
	""" Use RepeatMasker to find regions that should be filtered out """
	input:
		"data/Podan2/Podan2_AssemblyScaffoldsmt.fa"
	output:
		"RepeatMasker/Podan2_AssemblyScaffoldsmt.fa.out.gff"
	params:
		time = "2:00:00",
		threads = 16,
		TElib = TElib
	shell:
		""" 
		RepeatMasker -pa {params.threads} -a -xsmall -gccalc -gff -excln -lib {params.TElib} -dir RepeatMasker {input}
		"""

# ------- Filtering -------

rule snpsvcf:
	""" Filter resulting vcf file to get only SNPs """
	input:
		rawvcf = "results/PodoPop.vcf.gz"
	output:
		snpsvcf = "results/PodoPop-snps.vcf.gz",
	shell:
		"vcftools --gzvcf {input.rawvcf} --remove-indels --recode --recode-INFO-all --not-chr 'PaMt_NC_001329.3' --stdout | bgzip > {output.snpsvcf}" # Notice I filtered out the mitochondria too

rule bedtoolsTEs:
	""" Filter vcf with the repeats from RepeatMasker """
	input:
		vcf = "results/PodoPop-snps.vcf.gz",
		gfffile = "RepeatMasker/Podan2_AssemblyScaffoldsmt.fa.out.gff"
	output:
		filteredvcf = "results/PodoPop-snps-NoTEs.vcf",
	shell:
		"""
		# Get the header
		bcftools view -h {input.vcf} > {output.filteredvcf}
		
		# Filter out the repeats
		bedtools intersect -a {input.vcf} -b {input.gfffile} -v >> {output.filteredvcf}
		"""

rule bgzip_tabix:
	""" Compress and index vcf file """
	input:
		"results/PodoPop-snps-NoTEs.vcf"
	output:
		"results/PodoPop-snps-NoTEs.vcf.gz"
	shell:
		"bgzip {input}; "
		"tabix -p vcf {output}"


rule VariantFiltration:
	""" Use GATK's VariantFiltration to mark out sites by INFO annotations """
	input:
		vcf = "results/PodoPop-snps-NoTEs.vcf.gz",
	output:
		filteredvcf = temp("results/PodoPop-snps-NoTEs-gatk.vcf.gz"),
	shell:
		"""
		gatk VariantFiltration \\
		-V {input.vcf} \\
		-O {output.filteredvcf} \\
		--filter-expression "QD < 2.0" --filter-name "QD2" \\
		--filter-expression "FS > 60.0" --filter-name "FS60" \\
		--filter-expression "MQ < 40.0" --filter-name "MQ40" \\
		--filter-expression "QUAL < 30.0" --filter-name "QUAL30" \\
		--filter-expression "SOR > 3.0" --filter-name "SOR3" \\
		--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
		
		##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
		##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
		##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
		##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
		##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
		##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
		# RankSum annotations can only be calculated for REF/ALT heterozygous sites and therefore will be absent from records that do not present read counts towards heterozygous genotypes.
		##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
		##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">
		##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">	
		# In practice, we only filter out low negative values when evaluating
		# variant quality because the idea is to filter out variants for which the
		# quality of the data supporting the alternate allele is comparatively low.
		# The reverse case, where it is the quality of data supporting the reference
		# allele that is lower (resulting in positive ranksum scores), is not really
		# informative for filtering variants.
		##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
		"""

rule getPASSsites:
	""" Remove sites rejected in VariantFiltration """
	input:
		vcf = "results/PodoPop-snps-NoTEs-gatk.vcf.gz",
	output:
		vcf = "results/PodoPop-snps-NoTEs-gatkPASS.vcf"
	shell:
		"""
		# Get the header
		bcftools view -h {input.vcf} > {output.vcf}
		
		# Filter out the repeats
		bcftools view -H {input.vcf} | grep 'PASS' >> {output.vcf}
		"""

rule bgzip_tabix2:
	""" Compress and index vcf file """
	input:
		"results/PodoPop-snps-NoTEs-gatkPASS.vcf"
	output:
		"results/PodoPop-snps-NoTEs-gatkPASS.vcf.gz"
	shell:
		"bgzip {input}; "
		"tabix -p vcf {output}"

rule snpsvcfnomiss:
	""" Filter out sites with missing data """
	input:
		vcf = "results/PodoPop-snps-NoTEs-gatkPASS.vcf.gz"
	output:
		filteredvcf = "results/PodoPop-snps-NoTEs-gatkPASS-miss1.vcf.gz"
	shell:
		"vcftools --gzvcf {input.vcf} --max-missing 1 --recode --recode-INFO-all --stdout | bgzip > {output.filteredvcf}; "
		"tabix -p vcf {output.filteredvcf}"


rule filteringreport:
	""" Get some stats of the filtering """
	input:
		rawvcf = "results/PodoPop.vcf.gz",
		snpsvcf = "results/PodoPop-snps.vcf.gz",
		tevcf = "results/PodoPop-snps-NoTEs.vcf.gz",
		gatkvcf = "results/PodoPop-snps-NoTEs-gatkPASS.vcf.gz",
		miss1vcf = "results/PodoPop-snps-NoTEs-gatkPASS-miss1.vcf.gz",
	output:
		"results/filteringstats.txt"
	params:
		time = "30:00",
		threads = 1,
	shell:
		""" 
		printf "Total number of variants: " > {output}
		bcftools view -H {input.rawvcf} | wc -l >> {output}

		printf "Number of SNPs: " >> {output}
		bcftools view -H {input.snpsvcf} | wc -l >> {output}

		printf "Number of SNPs excluding TEs: " >> {output}
		bcftools view -H {input.tevcf} | wc -l >> {output}

		printf "Number of SNPs excluding TEs and bad sites: " >> {output}
		bcftools view -H {input.gatkvcf} | wc -l >> {output}

		printf "Number of SNPs excluding TEs, bad sites, no siblings and no missing data: " >> {output}
		bcftools view -H {input.miss1vcf} | wc -l >> {output}
		"""

# ------- Plotting the final figure for the paper -------

rule plotinR:
	""" Make PCAs figures """
	input:  # In this order!
		metadata,
		matingdata,
		"results/PodoPop-snps-NoTEs-gatkPASS-miss1.vcf.gz"
	output:
		"temp/PodoPop-snps-NoTEs-gatkPASS-NoSibs-miss1.gds",
		"results/Figures/PaPCA_all.pdf", # SNPs from all the chromosomes
		"results/Figures/PaPCA_corr.pdf", # Not used in the paper
		"results/Figures/PaPCA12.pdf",
		"results/Figures/PaPCA13.pdf",
		"results/Figures/PaPCA23.pdf",
		"results/Figures/Fig2_MatingVsSNPs.pdf",
		"results/Figures/FigS4_PCAs.pdf",
	conda: 
		"envs/snprelate.yaml"
	params:
		time = "45:00",
		threads = 2,
	script:
		rscript

