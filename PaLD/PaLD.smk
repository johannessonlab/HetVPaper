# -*- snakemake -*-

### PaLD: A pipeline to estimate linkage disequilibrium in *Podospora anserina*
#############################################################################

# Pipeline based on the linkage analyses of Hench et al. (2019)
# Inter-chromosomal coupling between vision and pigmentation genes during
# genomic divergence, Nature Ecology & evolution 

# * 2.2.8.1.interchromLD_global_subset.sh
# * 2.2.8.1.linkage_decay.sh

#############################################################################
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/08/07
# ---------------
# +++++++++++++++++++++++++++++++++++++++++++++++++
# Version 1

from glob import glob
import itertools

# -------------------------------------------------

notWa = config["notWa"]
groupV = config["groupV"]
groupA = config["groupA"]
vcf = config["vcf"] # The variants file

chrs = list(range(1, config["NCHR"] + 1)) # Number of chromosomes (assuming vcf has contigs names like "chromosome_1")
wins = list(range(1, config["NWIN"] + 1)) # Number of randomly sampled windows from each chromosome
WINLEN = config["WINLEN"] # Length of LD windows for decay
THIN = config["THIN"] # Thining of the SNPs for the intra and interchromosomal LD calculation

THINhet = config["THINhet"] # Thining for het-v and het-r areas
MAF = config["MAF"]
MAFhet = config["MAFhet"]

hetvSTART = config["hetvSTART"]
hetvEND = config["hetvEND"]
hetrSTART = config["hetrSTART"]
hetrEND = config["hetrEND"]

# Scripts
vcfconcat4LD = config["vcfconcat4LD"] # Manipulate vcf files
PaLDPlot = config["PaLDPlot"] # Plotting script LD decay
PaLDheatmap = config["PaLDheatmap"] # Plotting script LD heatmap

# -------------------------------------------------
# ----------
# Rules not submitted to a job
localrules: tabix, makefakevcf
# ----------

combchrs = list(itertools.combinations(chrs, 2)) # Make combinations of pairs
combchrs = list(zip(*combchrs)) # get them in two separate lists, one with first elements, one with the second elements


wildcard_constraints: # So it doesn't get confused with the expansion of wildcards
    thin="\d+"

rule all:
	input:
		# LD decay results
		"results/figures/LDdecayAllFit.pdf",
		"results/figures/LDdecayFitMG.pdf",

		# LD heatmap
		expand("results/figures/InterchrVR-LD_thin{thin}_maf{maf}.png", thin = THINhet, maf = MAFhet),
		expand("results/figures/LDchr{nchr}_thin{thin}_maf{maf}.png", thin = THIN, nchr = chrs, maf = MAF), # Individual plots for me
		expand("results/figures/LDchrALL_thin{thin}_maf{maf}.png", thin = THIN, maf = MAF)

# **************************
# ****** Prepare data ******
# **************************

def makelist(sampleslist, action):
	newstring = ""
	for sample in sampleslist:
		newstring += action + " " + sample + " "
	return(newstring)

rule preparevcf:
	""" Remove non Wageningen samples, and break in samples """
	input:
		vcf
	output:
		"data/vcfs/HighQWageningen.vcf.gz",
		"data/vcfs/HighQWageningen-groupV.vcf.gz", # V
		"data/vcfs/HighQWageningen-groupA.vcf.gz", # V1
	params:
		time = "30:00",
		threads = 1,
	run:
		cmd = "vcftools --gzvcf " + input[0] + " --recode --recode-INFO-all --stdout " + makelist(notWa, "--remove-indv") + " | bgzip > " + output[0]
		shell(cmd)

		cmd = "vcftools --gzvcf " + input[0] + " --recode --recode-INFO-all --stdout " + makelist(groupV, "--indv") + " | bgzip > " + output[1]
		shell(cmd)

		cmd = "vcftools --gzvcf " + input[0] + " --recode --recode-INFO-all --stdout " + makelist(groupA, "--indv") + " | bgzip > " + output[2]
		shell(cmd)

rule tabix:
	""" Sanity check, it shouldn't complain, and BCFtools might need it """
	input:
		allwa = "data/vcfs/HighQWageningen.vcf.gz",
		grV = "data/vcfs/HighQWageningen-groupV.vcf.gz",
		grA = "data/vcfs/HighQWageningen-groupA.vcf.gz",
	output:
		allwa = "data/vcfs/HighQWageningen.vcf.gz.tbi",
		grV = "data/vcfs/HighQWageningen-groupV.vcf.gz.tbi",
		grA = "data/vcfs/HighQWageningen-groupA.vcf.gz.tbi",
	shell:
		"tabix -p vcf {input.allwa}; "
		"tabix -p vcf {input.grV}; "
		"tabix -p vcf {input.grA}; "


# **************************
# ****** LD decay ******
# **************************

rule randomWindLD:
	""" Calculate LD in randomly selected windows """
	# Modified 2.2.8.1.linkage_decay.sh from Hench et al (2019)
	# http://tldp.org/LDP/abs/html/randomvar.html
	# http://srufaculty.sru.edu/david.dailey/unix/random_numbers.htm
	input: 
		allwa = "data/vcfs/HighQWageningen.vcf.gz",
		grV = "data/vcfs/HighQWageningen-groupV.vcf.gz",
		grA = "data/vcfs/HighQWageningen-groupA.vcf.gz",
		allwa_index = "data/vcfs/HighQWageningen.vcf.gz.tbi",
		grV_index = "data/vcfs/HighQWageningen-groupV.vcf.gz.tbi",
		grA_index = "data/vcfs/HighQWageningen-groupA.vcf.gz.tbi",
	output:
		allwa = "LDdecay/dist-all.chr{nchr}.w{win}.hap.ld",
		grV = "LDdecay/dist-grV.chr{nchr}.w{win}.hap.ld",
		grA = "LDdecay/dist-grA.chr{nchr}.w{win}.hap.ld",
	params:
		winlen = WINLEN,
		time = "6:00:00",
		threads = 1,
	shell:
		"""
		# Get length of this chromosome 
		lenrhischr=$(bcftools view -h {input.allwa} | grep "##contig=<ID=chromosome_{wildcards.nchr}" | sed -r 's;.*(length=)([0-9]*)>;\\2;')

		RANDOM=$(date +%N) #set seed using the current time in nanoseconds
		# $RANDOM will be in the interval 0 to 32,767
		length=${{#lenrhischr}}
		
		## Make a random name much larger than that using the length of the chromosome as max (Hench didn't do this)
		# bc is a calculator, I get a floating point and then I remove the 0s to turn it into a integer
		RANDOM2=$( echo "scale=$length; $RANDOM/32767" | bc | sed -r 's/\\.//' | sed 's/0//g' ) #  scale is the number of digits after a decimal point

		RN1=$(( ( RANDOM2 % lenrhischr ) + 1)) # If you need a random int within a certain range, use the 'modulo' operator.
		RN2=$(( $RN1+{params.winlen} ))
		echo "** range $RN1 $RN2 in chromosome {wildcards.nchr} of length $lenrhischr **"

		# The minimum MAF is tailored to remove singletons
		vcftools --gzvcf {input.allwa} --chr chromosome_{wildcards.nchr} --from-bp $RN1 --to-bp $RN2 --hap-r2 --stdout --maf 0.01 > {output.allwa} # all samples together
		vcftools --gzvcf {input.grV} --chr chromosome_{wildcards.nchr} --from-bp $RN1 --to-bp $RN2 --hap-r2 --stdout --maf 0.025 > {output.grV} # Mating group V
		vcftools --gzvcf {input.grA} --chr chromosome_{wildcards.nchr} --from-bp $RN1 --to-bp $RN2 --hap-r2 --stdout --maf 0.025 > {output.grA} # Mating group A
		""" 

rule catLDdecayPerChrAll:
	""" Put together the LD files """
	input:
		expand("LDdecay/dist-all.chr{{nchr}}.w{win}.hap.ld", win = wins)
	output:
		"LDdecay/dist-all.chr{nchr}.hap.ld.gz"
	params:
		time = "10:00",
		threads = 1,
	shell:
		"cat {input} | bgzip > {output}"

rule catLDdecayPerChrgrV:
	""" Put together the LD files """
	input:
		expand("LDdecay/dist-grV.chr{{nchr}}.w{win}.hap.ld", win = wins) ## the double braces signals a wildcard
	output:
		"LDdecay/dist-grV.chr{nchr}.hap.ld.gz"
	params:
		time = "10:00",
		threads = 1,
	shell:
		"cat {input} | bgzip > {output}"

rule catLDdecayPerChrgrA:
	""" Put together the LD files """
	input:
		expand("LDdecay/dist-grA.chr{{nchr}}.w{win}.hap.ld", win = wins) ## the double braces signals a wildcard
	output:
		"LDdecay/dist-grA.chr{nchr}.hap.ld.gz"
	params:
		time = "10:00",
		threads = 1,
	shell:
		"cat {input} | bgzip > {output}"

rule plotLDdecay:
	""" Plot LD decay for all samples """
	input:
		expand("LDdecay/dist-all.chr{nchr}.hap.ld.gz", nchr = chrs), # there should be 7 objects in list (chromosomes)
		expand("LDdecay/dist-grV.chr{nchr}.hap.ld.gz", nchr = chrs),
		expand("LDdecay/dist-grA.chr{nchr}.hap.ld.gz", nchr = chrs),
	output:
		"results/figures/LDdecayAllFit.pdf", # LD decay curve with raw data like Hench et al (2019), the Remington's non-linear regression, and the average r2 in bins with a smoothing line 
		"results/figures/LDdecayFitMG.pdf", # LD decay curve with just the Remington's fit, but for all and for each mating group
		"results/stats/LDdecayMid.csv", # A little table with aprox LD values at around r^2 = 0.2 per chromosome
	params:
		time = "10:00",
		threads = 1,
	conda: 
		"envs/PaLDPlot.yaml"
	script:
		PaLDPlot

# **************************
# ****** LD heatmap ******
# **************************

rule filtervcf:
	""" Prepare a vcf file with selected SNPs for whole genome LD calculations (with thining and minimum MAF) """
	input:
		"data/vcfs/HighQWageningen.vcf.gz",
	output:
		"data/vcfs/HighQWageningen-thin{thin}.vcf.gz",
	params:
		time = "30:00",
		threads = 1,
		thin = THIN,
		maf = MAF,
	shell:
		"vcftools --gzvcf {input} --recode --recode-INFO-all --stdout --thin {params.thin} --maf {params.maf} | bgzip > {output}" 
		# The minimum MAF is tailored to remove singletons


rule extractVCFhets:
	""" Extract the SNPs from around het-v and het-r into a vcf file """
	input:
		"data/vcfs/HighQWageningen.vcf.gz",
	output:
		"data/vcfs/HighQWageningen-thin{thin}_maf{maf}-hetVR.vcf.gz",
	params:
		time = "30:00",
		threads = 1,
		thin = THINhet,
		hetvSTART = hetvSTART,
		hetvEND = hetvEND,
		hetrSTART = hetrSTART,
		hetrEND = hetrEND,
		maf = MAF,
	shell:
		""" 
		printf "chrom\\tchromStart\\tchromEnd\\n" > hetVR.bed
		printf "chromosome_2\\t{params.hetrSTART}\\t{params.hetrEND}\\n" >> hetVR.bed
		printf "chromosome_5\\t{params.hetvSTART}\\t{params.hetvEND}\\n" >> hetVR.bed

		vcftools --gzvcf {input} --recode --recode-INFO-all --stdout --thin {params.thin} --maf {params.maf} --bed hetVR.bed | bgzip > {output}

		rm hetVR.bed
		"""

rule makefakevcf:
	""" Create a fake vcf file where all SNPs are in a single fake contig """
	input:
		"data/vcfs/HighQWageningen-thin{thin}_maf{maf}-hetVR.vcf.gz",
	output:
		"data/vcfs/HighQWageningen-thin{thin}_maf{maf}-hetVR-fakecoords.vcf.gz",
	params:
		vcfconcat4LD = vcfconcat4LD
	shell:
		"{params.vcfconcat4LD} {input} --fakectg hetVR --buffer 100000 | bgzip > {output}"

rule LDhetVR:
	""" Calculate r2 with SNPs within and between the het-v and het-r loci """
	input:
		"data/vcfs/HighQWageningen-thin{thin}_maf{maf}-hetVR-fakecoords.vcf.gz",
	output:
		"LDinterchr/LDhetVR-thin{thin}_maf{maf}-fakecoords.ld.gz",
	params:
		time = "6:00:00",
		threads = 2, # For memory
	shell:
		"vcftools --gzvcf {input} --hap-r2 --stdout | bgzip > {output}" 

rule intrachromLD:
	""" Calculate r2 with SNPs sampled from within a chromosome """
	input:
		"data/vcfs/HighQWageningen-thin{thin}.vcf.gz",
	output:
		"LDinterchr/HighQWageningen-thin{thin}_maf{maf}.intrachrom{nchr}.hap.ld.gz"
	params:
		time = "2:00:00",
		threads = 2,
		maf = MAF,
	shell:
		"vcftools --gzvcf {input} --stdout --hap-r2 --chr chromosome_{wildcards.nchr} --maf {params.maf} | bgzip > {output}" 

rule plotLDheatmap:
	""" Plot LD heatmaps for intra and inter chromosomal comparisons """
	input:
		expand("LDinterchr/LDhetVR-thin{thin}_maf{maf}-fakecoords.ld.gz", thin = THINhet, maf = MAFhet),
		expand("LDinterchr/HighQWageningen-thin{thin}_maf{maf}.intrachrom{nchr}.hap.ld.gz", thin = THIN, nchr = chrs, maf = MAF)
	output:
		expand("results/figures/InterchrVR-LD_thin{thin}_maf{maf}.png", thin = THINhet, maf = MAFhet),
		expand("results/figures/LDchr{nchr}_thin{thin}_maf{maf}.png", thin = THIN, nchr = chrs, maf = MAF), # Individual plots for me
		expand("results/figures/LDchrALL_thin{thin}_maf{maf}.png", thin = THIN, maf = MAF)
	params:
		time = "2:00:00",
		threads = 5, # Give it more memory
	conda: 
		"envs/PaLDPlot.yaml"
	script:
		PaLDheatmap	

