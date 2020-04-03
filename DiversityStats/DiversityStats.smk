# -*- snakemake -*-

### DiversityStats: Calculate population genetic diversity statistics from Podospora data
#############################################################################

#############################################################################
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/07/04
# ---------------
# +++++++++++++++++++++++++++++++++++++++++++++++++
# Version 1

# -------------------------------------------------
# # List of samples to analyze
samples = config["SampleIDs"]
# The variants file
vcf = config["vcf"]
# The gff of sites to replace for missing data
TEgff = config["TEgff"]
# Path to the BAM files
path2BAM = config["path2BAM"]

# Scripts
badsites2vcf = config["badsites2vcf"]
totalcovergff = config["totalcovergff"]
DiversityStats_vcfR_plotter = config["DiversityStats_vcfR_plotter"]
DiversityStatsFstInterv = config["DiversityStatsFstInterv"]
DiversityStatsCalc = config["DiversityStatsCalc"]
DiversityStatsCalcPlot = config["DiversityStatsCalcPlot"]

# Metadata file
# metadata = config["metadata"]
hetgenes = config["hetgenes"]
# Number of chromosomes (assuming vcf has contigs names like "chromosome_1")
chrs = list(range(1, config["NCHR"] + 1))

# -------------------------------------------------

# Get name of file
input_base = os.path.basename(vcf) # Remove the path
input_name = input_base.split(".")[0] # Taking out the prefix of the file

# ----------
# Rules not submitted to a job
localrules: filtergenomecov, totalcovergff, removeoverlap, inputmissingsites, compressvcf, tabix, extractchrgff, extractchrcov, plotstats, getannotation
# ----------

wildcard_constraints: # So it doesn't get confused with the expansion of wildcards for chromosome 5
    nchr="\d+"

rule all:
	input:
		expand("results/FstReplicas/FstConfIntv_chr{nchr}.csv", nchr = chrs), 

		# Final plot of figures
		"results/figures/Fig1_PopData_Main_2Chrs.pdf", # MAIN figure 
		"results/figures/FigS2_PopData_Supp.pdf", # MAIN figure 
		"results/figures/Fig3_FstAll.pdf", # MAIN figure 
		"results/figures/FigS5_Dxy_supp.pdf", # MAIN figure 
		"results/figures/Fig4_HetVRgenotypes.pdf", # MAIN figure
		"results/stats/highTajima.csv", # Useful here

rule getannotation:
	""" Get Podan2 annotation """
	output:
		"data/Podan2_genes.gff"
	shell:
		"wget https://raw.githubusercontent.com/johannessonlab/SpokPaper/master/Fig1_3SppNetwork/references/Podan2_AssemblyScaffoldsmtGenesEd_gh.gff -O {output}"

# ------- Find sites that are not good -------

rule PlotCoverage:
	""" Plot coverage distribution for all samples and report """
	input:  # In this order!
		vcf
	output:
		"results/figures/Coverage.pdf",
		"data/cov/CoverageSummary.tab",
	params:
		time = "1:00:00",
		threads = 2, # It needs some memory, so I have to give it two
		lquantile = 0.25,
		uquantile = 0.985,
	conda: 
		"envs/statsplot.yaml"
	script:
		DiversityStats_vcfR_plotter

rule genomecov:
	""" genomecov makes BED with intervals of coverage """ 
	input:
		path2BAM + "/{sample}/{sample}-to-Podan2.sorted.debup.bam",
	output:
		"genomecov/{sample}_cov.bed"
	params:
		time = "15:00",
		threads = 1,
	shell:
		"bedtools genomecov -ibam {input} -bga > {output}"
		# bga - report coverage in intervals, including sites with 0 coverage

rule filtergenomecov:
	""" Filter the coverage sites based on the quantile of each sample's distribution """
	input:
		covsum = "data/cov/CoverageSummary.tab",
		bed = "genomecov/{sample}_cov.bed",
	output:
		"genomecov/{sample}_cov_filtered.bed"
	shell:
		"""
		lowquantile=$(grep {wildcards.sample} {input.covsum} | cut -f4)
		uppquantile=$(grep {wildcards.sample} {input.covsum} | cut -f5)

		awk -v low="$lowquantile" -v upp="$uppquantile" '$4 < low || $4 > upp' {input.bed} > {output}
		"""

rule totalcovergff:
	""" Make bed file of sites overlapping transposable elements """
	input:
		TEgff
	output:
		"data/cov/TEsites.bed"
	params:
		totalcovergff = totalcovergff
	shell:
		"{params.totalcovergff} {input} > {output}"

rule BEDOPS:
	""" Take all the BED files and overlap them into one single BED file"""
	# https://bedops.readthedocs.io/en/latest/content/overview.html#overview
	# https://bedops.readthedocs.io/en/latest/content/reference/set-operations/bedops.html#merge-m-merge
	input:
		"data/cov/TEsites.bed",
		expand("genomecov/{sample}_cov_filtered.bed", sample = samples) 
	output:
		temp("filter/excludedsites-raw.bed")
	params:
		time = "1:00:00",
		threads = 1,
	shell:
		"bedops --merge {input} > {output}"

# ----
def remove_overlap(ranges):
	""" Simplify a list of ranges; I got it from https://codereview.stackexchange.com/questions/21307/consolidate-list-of-ranges-that-overlap """
	result = []
	current_start = -1
	current_stop = -1 

	for start, stop in sorted(ranges):
		if start > current_stop:
			# this segment starts after the last segment stops
			# just add a new segment
			result.append( (start, stop) )
			current_start, current_stop = start, stop
		else:
			# current_start already guaranteed to be lower
			current_stop = max(current_stop, stop)
			# segments overlap, replace
			result[-1] = (current_start, current_stop) # SLAV: I modified this to update the stop too.
	return(result)
# ----

rule removeoverlap:
	""" Fuse overlapping ranges in a BED file """
	# Basically the same as the script totalcovergff.py but for bed files
	input:
		"filter/excludedsites-raw.bed"
	output:
		"filter/excludedsites.bed"
	run:
		# Make a dictionary of ranges
		beddic = {} # key: chromosome, value: (start, end)

		with open(input[0], 'r') as file: 
			for line in file:
				if '#' in line:
					pass
				elif line not in ['\n', '\r\n']: # Ignore empty lines
					cols = line.rstrip("\n").split("\t")

					contig = cols[0]
					start = int(cols[1])
					end = int(cols[2])

					if contig in list(beddic.keys()):
						beddic[contig].append([start, end])
					else: # contig is new
						beddic[contig] = [[start, end]]

			# Reduce the overlaps
			for ctg in beddic.keys():
				beddic[ctg] = remove_overlap(beddic[ctg])

			# Print them in a new file
			with open(output[0], 'w') as result:
				result.write('#Contig\tStart\tEnd\n') # header
				for ctg in beddic.keys():
					for interval in beddic[ctg]:
						result.write('{}\t{}\t{}\n'.format(ctg, interval[0], interval[1]))

rule inputmissingsites:
	""" Input the sites that should be excluded as explicit missing data sites """
	input:
		vcf = vcf,
		# vcf = "filter/" + input_name + "-NoSibs.vcf.gz",
		bed = "filter/excludedsites.bed",
	output:
		"filter/" + input_name + "-NoSibs-withNA.vcf"
	params:
		badsites2vcf = badsites2vcf
	shell:
		"{params.badsites2vcf} {input.vcf} {input.bed} --bed > {output}"

rule compressvcf:
	""" Compress and index vcf file """
	input:
		"filter/" + input_name + "-NoSibs-withNA.vcf"
	output:
		"filter/" + input_name + "-NoSibs-withNA.vcf.gz"
	shell:
		"bgzip {input}"

rule tabix:
	""" Sanity check, it shouldn't complain """
	input:
		"filter/" + input_name + "-NoSibs-withNA.vcf.gz"
	output:
		"filter/" + input_name + "-NoSibs-withNA.vcf.gz.tbi"
	shell:
		"tabix -p vcf {input}" # This is a sanity check, it shouldn't complain	

# ------- Prepare the files for PopGenome -------

rule extractchrvcf:
	""" Extract one individual chromosome from the vcf file """
	input:
		vcf = "filter/" + input_name + "-NoSibs-withNA.vcf.gz",
		tbi = "filter/" + input_name + "-NoSibs-withNA.vcf.gz.tbi" # just to force the previous rule to produce it
	output:
		"data/chr{nchr}/{input_name}-NoSibs-withNA-chr{nchr}.vcf"
	params:
		time = "15:00",
		threads = 1, 
	shell:
		"vcftools --gzvcf {input.vcf} --recode --recode-INFO-all --chr chromosome_{wildcards.nchr} --stdout > {output}" # Notice it has to be uncompressed for PopGenome


rule extractchrgff:
	""" Extract one individual chromosome from the gff file """
	input:
		"data/Podan2_genes.gff"
	output:
		"data/gffs/{input_name}-NoSibs-withNA-chr{nchr}.gff"
	shell:
		"""
		grep '#' {input} > {output}
		grep -E "$(printf 'chromosome_{wildcards.nchr}\\t')" {input} >> {output}
		""" 

rule extractchrcov:
	""" Extract one individual chromosome from the bed file """
	input:
		"filter/excludedsites.bed"
	output:
		"data/cov/excludedsites-chr{nchr}.bed"
	shell:
		"""
		head -n1 {input} > {output}
		grep 'chromosome_{wildcards.nchr}' {input} >> {output}
		""" 

rule removehybridschr5:
	""" Remove samples with recombinant haplotypes of het-v"""
	input:
		vcf = "data/chr5/{input_name}-NoSibs-withNA-chr5.vcf",
		gff = "data/gffs/{input_name}-NoSibs-withNA-chr5.gff",
	output:
		vcf = "data/chr5/{input_name}-NoSibs-withNA-chr5-nohybrids.vcf",
		gff = "data/gffs/{input_name}-NoSibs-withNA-chr5-nohybrids.gff",
	params:
		time = "15:00",
		threads = 1, 
	shell:
		"vcftools --gzvcf {input.vcf} --recode --recode-INFO-all --remove-indv PaWa60p --remove-indv PaWa61m --stdout > {output.vcf}; " # Identified based on the PCA
		"cp {input.gff} {output.gff}" # It has to be the same name as the vcf...

# ------- Calculate the statistics -------

rule FstConfIntervals:
	""" Calculate the confidence intervals through permutation of Fst """
	input: # In this order!
		hetgenes, #metadata,
		"data/chr{nchr}/" + input_name + "-NoSibs-withNA-chr{nchr}.vcf",
		"data/gffs/" + input_name + "-NoSibs-withNA-chr{nchr}.gff",
	output:
		"data/FstReplicas/FstReplicas_chr{nchr}.csv",
		"results/FstReplicas/FstConfIntv_chr{nchr}.csv",
	params:
		time = "3-00:00:00",
		threads = 1, # I need to give these many in case memory is needed	
	conda: 
		"envs/statsplot.yaml"
	script:
		DiversityStatsFstInterv

rule CalcStats:
	""" Calculate diversity statistics """
	input:  # In this order!
		hetgenes, #metadata,
		"data/cov/excludedsites-chr{nchr}.bed",
		"data/chr{nchr}/{input_name}-NoSibs-withNA-chr{nchr}.vcf",
		"data/gffs/{input_name}-NoSibs-withNA-chr{nchr}.gff",
	output:
		"results/stats/{input_name}-NoSibs-withNA-chr{nchr}-stats.csv",
	params:
		time = "1-00:00:00",
		threads = 6, # I need to give these many in case memory is needed
	script:
		DiversityStatsCalc

rule CalcStatsChr5:
	""" Calculate diversity statistics but for chromosome 5 without the recombinant samples"""
	input:  # In this order!
		hetgenes, #metadata,
		"data/cov/excludedsites-chr5.bed",
		"data/chr5/{input_name}-NoSibs-withNA-chr5-nohybrids.vcf",
		"data/gffs/{input_name}-NoSibs-withNA-chr5-nohybrids.gff",
	output:
		"results/stats/{input_name}-NoSibs-withNA-chr5-nohybrids-stats.csv",
	params:
		time = "1-00:00:00",
		threads = 6, # I need to give these many in case memory is needed
	script:
		DiversityStatsCalc

rule plotstats:
	""" Plot diversity statistics """
	input: # In this order!
		expand("results/stats/{input_name}-NoSibs-withNA-chr{nchr}-stats.csv", nchr = chrs, input_name = input_name),
		expand("results/FstReplicas/FstConfIntv_chr{nchr}.csv", nchr = chrs),
		hetgenes,
		expand("results/stats/{input_name}-NoSibs-withNA-chr5-nohybrids-stats.csv", input_name = input_name),
	output:
		main2chrs = "results/figures/Fig1_PopData_Main_2Chrs.pdf", # MAIN figure
		mainSupchrs = "results/figures/FigS2_PopData_Supp.pdf", # MAIN figure 
		hetrvtime = "results/figures/Fig4_HetVRgenotypes.pdf", # MAIN figure
		mainfst = "results/figures/Fig3_FstAll.pdf", # MAIN figure 
		supdxy = "results/figures/FigS5_Dxy_supp.pdf", # MAIN figure 
		extratajima = "results/figures/Tajima_supp.pdf", # Useful for me
		tajimatable = "results/stats/highTajima.csv", # Useful here
	conda: 
		"envs/statsplot.yaml"
	params:
		time = "15:00",
		threads = 1, 
	script:
		DiversityStatsCalcPlot



