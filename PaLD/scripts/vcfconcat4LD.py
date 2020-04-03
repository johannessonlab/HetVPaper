#!/usr/bin/env python
# encoding: utf-8

# ================== vcfconcat4LD =================

# In order to facilitate the calculation with VCFTools and the plotting of
# pairwise r2 (LD) between SNPs (i.e. an LD heatmap) between chromosomes, 
# this script modifies the input vcf file in such a way that all SNPs have 
# the same "fake" chromosome, and the coordinates become accumulative. So 
# the second chromosome coordinates continue where the first chromosome ended.

# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/06/25
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import argparse  # For the fancy options
import sys  # To exit the script, and to pipe out
import os  # For the input name
import gzip # For unzipping files
from datetime import datetime
# from Bio import SeqIO # To deal with vcf files compressed with bgzip
# from Bio import bgzf # To deal with vcf files compressed with bgzip
# ------------------------------------------------------

version = 1
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* badsites2vcf *")  # Create the object using class argparse

# Add options
parser.add_argument('vcf', help="Standard (sorted) vcf file compressed with bgzip that contains ONLY the sites desired")
parser.add_argument('--uncompressed', '-z', help="VCF file is not compressed", default=False, action='store_true')
parser.add_argument('--buffer', '-b', help="Fake space in between chromosomes in bp", default=200000, type=int)
parser.add_argument('--fakectg', '-f', help="Name of the fake new contig for all sites", default="contig")

parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	if args.uncompressed:
		vcfopen = open(args.vcf, 'r')
	else:
		vcfopen = gzip.open(args.vcf, 'rt') # 't' is text mode, to interpret the tabs and new lines

except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()
# ============================

# ============================
# Read vcf file
# ============================

contig = ''
header = False

for line in vcfopen:
	if "##" in line: # header
		if "#contig" in line: 
			if not header:
				sys.stdout.write('##contig=<ID=' + args.fakectg + ',length=100000000>\n') # Some fake big number
				header = True
		else:
			sys.stdout.write(line) 
	elif "#CHROM" in line: # column names
		sys.stdout.write('##FILTER=<ID=vcfconcat4LD,Description="The contigs no longer reflect the truth, instead all sites are fused into a single fake contig">\n')
		sys.stdout.write('##INFO=<ID=CHR,Number=.,Type=String,Description="The real chromosome/contig that the SNP belongs to">\n') # Some fake big number # NOT WORKING

		now = datetime.now()
		sys.stdout.write('##vcfconcat4LD.py=v.' + versiondisplay + ' on ' + now.strftime("%d/%m/%Y %H:%M:%S") + '\n')

		sys.stdout.write(line) 
		headnext = line.rstrip("\n").split('\t')
		samplesIDs = headnext[9:]
		numberofsamples = len(headnext[9:])

		# pass
	else:
		# Extract information from the line
		columns = line.rstrip("\n").split('\t')
		CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = columns[0:9]
		POS = int(POS)
		sites = columns[9:]

		# Add the real chromosome to the INFO field to recover it later if needed
		INFO += ";CHR=" + CHROM # not working

		# Is it the very first site and contig?
		if contig == '': 
			contig = CHROM
			lastPOS = 0
			lastsite = 0

		# Next sites
		if contig != CHROM: # we are in the next contig now
			contig = CHROM
			lastsite = lastPOS - POS

		newPOS = POS + lastsite + args.buffer

		# Let's assume that the last SNP is close to the border (so do not worry about the actual size of the chromosome)
		newsite = args.fakectg + '\t' + str(POS + lastsite) + '\t' + '\t'.join(columns[2:7]) + '\t' + INFO + '\t' + FORMAT + '\t' + '\t'.join(sites) + '\n'
		sys.stdout.write(newsite)
		
		lastPOS = newPOS # So it remembers the previous site

		# print(contig, POS, lastsite) # debugging

