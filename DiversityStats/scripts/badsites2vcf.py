#!/usr/bin/env python
# encoding: utf-8

# ================== badsites2vcf =================
# Script to incorporate sites that have missing data but are not part of
# the standard vcf file, and hence are treated as invariants traditionally. It
# can take a gff or a bed file as input to get ranges of missing data.

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

version = 1.2
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* badsites2vcf *")  # Create the object using class argparse

# Add options
parser.add_argument('vcf', help="Standard vcf file compressed with bgzip")
parser.add_argument('gff', help="Gff file with regions that should be considered missing data, eg. transposable elements in the reference")
parser.add_argument('--uncompressed', '-z', help="VCF file is not compressed", default=False, action='store_true')
parser.add_argument('--bed', '-b', help="The annotation file is a bed file, not a gff", default=False, action='store_true')
# parser.add_argument('--gzvcf', '-z', help="VCF file is compressed", default=False, action='store_true')

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
	gffopen = open(args.gff, 'r')

except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()
# ============================

# ============================
# Read gff with regions to be excluded
# ============================
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


# Save the ranges in a dictionary
badsites = {} # key: chromosome, value: (start, end)

if args.bed: # bed
	startindex = 1
	endindex = 2
else: # gff
	startindex = 3
	endindex = 4

for line in gffopen:
	if '#' in line:
		pass
	elif line not in ['\n', '\r\n']: # Ignore empty lines
		cols = line.rstrip("\n").split("\t")

		contig = cols[0]
		start = int(cols[startindex])
		end = int(cols[endindex])

		if contig in list(badsites.keys()):
			badsites[contig].append((start, end))
		else:
			badsites[contig] = [(start, end)]

# Reduce the overlaps
for ctg in badsites.keys():
	badsites[ctg] = remove_overlap(badsites[ctg])

# ============================
# Read vcf file
# ============================


for line in vcfopen:
	if "##" in line: # header
		sys.stdout.write(line) 
		# pass	
	elif "#CHROM" in line: # column names
		sys.stdout.write('##FILTER=<ID=badsites2vcf,Description="Sites overlapping with gff ' + args.gff + '">\n')

		now = datetime.now()
		sys.stdout.write('##badsites2vcf.py=v. ' + versiondisplay + ' on ' + now.strftime("%d/%m/%Y %H:%M:%S") + '\n')

		sys.stdout.write(line) 
		headnext = line.rstrip("\n").split('\t')
		samplesIDs = headnext[9:]
		numberofsamples = len(headnext[9:])

	else:
		# Extract information from the line
		columns = line.rstrip("\n").split('\t')
		CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = columns[0:9]
		POS = int(POS)
		sites = columns[9:]

		# *** Get format of site ***
		# --------------------------------------
		# This can be determined from the header, but for now I'm just gonna go with
		# the GATK format for simplicity

		# formfields = FORMAT.split(":")

		# Make a NA string out of the GATK format for that number of samples
		nastring = ".:0,0:0:.:0,0" + '\t'
		samplesnastrings = nastring * numberofsamples
		samplesnastrings = samplesnastrings.rstrip("\t")
		# --------------------------------------
		
		# *** Add the ranges of bad sites until you find the variant site ***
		if CHROM in list(badsites.keys()): # otherwise ignore
			if badsites[CHROM] == []: # The very last sites after all bad ranges (tuples)
				sys.stdout.write(line) # Print the site
				continue # Go to next variant
			else: # work to do
				lastend = int(badsites[CHROM][-1][1])

			if POS > lastend: # It's above all tuples so empty them out for this contig
				# If you haven't already, print all the ranges as missing data
				if len(badsites[CHROM]) > 0:
					for rangetuble in badsites[CHROM]:
						start = int(rangetuble[0])
						end = int(rangetuble[1])
						
						for site in range(start, end + 1): # The vcf is base 1
							newsite = CHROM + '\t' + str(site) + '\t.\tN\tN\t0\tbadsites2vcf\tDP=0\t' + FORMAT + "\t" + samplesnastrings + "\n" 
							sys.stdout.write(newsite)
							# print("...etc..."); break # debugging

					badsites[CHROM] = [] # erase the tuples because they have been printed
				sys.stdout.write(line) # Print the site

			else: # It's not the end of the file
				for rangetuble in badsites[CHROM]:
					start = int(rangetuble[0])
					end = int(rangetuble[1])

					# Record also the start of the next tuple
					rangetupleindex = badsites[CHROM].index(rangetuble)
					lenbadsiteschr = len(badsites[CHROM])

					if rangetupleindex + 1 < lenbadsiteschr:
						# print(badsites[CHROM], badsites[CHROM][rangetupleindex])
						nextstart = int(badsites[CHROM][rangetupleindex + 1][0])
						nextend = int(badsites[CHROM][rangetupleindex + 1][1])
					else: # This is the last tuple
						nextstart = start
						nextend = end

					# print("New tuple: pos ", POS, start, end, "\t", nextstart, nextend) # debugging

					## A
					if (POS < start): # and POS < end and POS < nextstart): # Site is before the tuple
						# print("A") # debugging
						sys.stdout.write(line) # Print the site as normal
						break # Go to next variable site

					## B
					elif (POS >= start and POS < end): # and POS < nextstart): # this site is within a bad track (tuple)
						# print("B") # debugging
						for site in range(start, POS): # The vcf is base 1
							newsite = CHROM + '\t' + str(site) + '\t.\tN\tN\t0\tbadsites2vcf\tDP=0\t' + FORMAT + "\t" + samplesnastrings + "\n" 
							sys.stdout.write(newsite)
							# print("...etc..."); break # debugging
							
						# The actual site
						FILTER = "badsites2vcf\t"
						# Construct the info for site		
						linestring = '\t'.join(columns[0:6]) + "\t" + FILTER + '\t'.join(columns[7:9]) + "\t" + samplesnastrings + "\n" 
						sys.stdout.write(linestring)

						# Update the range tuple	
						rangetupleindex = badsites[CHROM].index(rangetuble)
						badsites[CHROM][rangetupleindex] = (POS + 1, end)

						break # Go to next variable site			

					## C
					elif (POS > start and POS >= end and POS < nextstart): # The site is after the tuple 
						# print("C") #, badsites[CHROM]) # debugging
						for site in range(start, end + 1): # The vcf is base 1
							newsite = CHROM + '\t' + str(site) + '\t.\tN\tN\t0\tbadsites2vcf\tDP=0\t' + FORMAT + "\t" + samplesnastrings + "\n"  # The tuple should all be NA
							sys.stdout.write(newsite)
							# print("...etc..."); break # debugging

						sys.stdout.write(line) # Print the site as normal

						# This tuple can be removed now
						badsites[CHROM].remove(rangetuble)

						break # Go to next variable site

					## D + E
					elif (POS > end and POS >= nextstart): # This tuple should all be NA, but the variant may or may not overlap the next tuple
						# print("D+E") #, badsites[CHROM]) # debugging
						for site in range(start, end + 1): # The vcf is base 1
							newsite = CHROM + '\t' + str(site) + '\t.\tN\tN\t0\tbadsites2vcf\tDP=0\t' + FORMAT + "\t" + samplesnastrings + "\n"  # The tuple should all be NA
							sys.stdout.write(newsite)
							# print("...etc..."); break # debugging

						# All the tuples start from the next one
						rangetupleindex = badsites[CHROM].index(rangetuble)
						badsites[CHROM] = badsites[CHROM][rangetupleindex + 1:]

						# Is the next tuple also D (< nextend) or E?
						if (POS <= nextend) or (POS > nextend): 
							# finish the tuple
							continue # go to next tuple

						sys.stdout.write(line) # Print the site as normal
						break # Go to next variant

					## C but last tuple
					elif (POS > start and POS > end and start == nextstart): # There are no more tuples coming
						sys.stdout.write(line) # Print the site as normal

						break # Go to next variable site
