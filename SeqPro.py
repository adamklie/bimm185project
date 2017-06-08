#! /usr/bin/env python
import subprocess

# parse user arguments 
import argparse
import sys
parser = argparse.ArgumentParser()
parser.add_argument("-f1", "--fastq1", required=True, help="forward fastq reads from bacterial sample", type=str)
parser.add_argument("-f2", "--fastq2", required=False, help="reverse fastq reads from bacterial sample    ", type=str)
parser.add_argument("-rt", "--read_type", required=True, help="type of read library input, s=single end, p=paired end, m=mate pair")
parser.add_argument("-r", "--reference", help="optional reference fasta to perform assembly against", required=False, type=str)
parser.add_argument("-v", "--verbosity", action="store_true", help="increase output verbosity")
pair = False
args = parser.parse_args()
if args.verbosity:
	log = sys.stderr
else:	
	log = open("log.txt", "w")

#check read type
if args.read_type in "pm":
	pair =  True
	if args.fastq2 is None:
		print("read pair not specified, exiting...")
		sys.exit()

#spades assembly with no reference
if args.reference is None:
	print("No reference given, performing de novo assembly")
	if args.read_type == "s":
		print("single end type assembly")
	elif args.read_type == "p":
		print("paired end type assembly")
	elif args.read_type == "m":
		print("mate pair type assembly")
	else:
		print("unrecognizable pair type, exiting...")
		sys.exit()
"""
#reference based assembly with reference given
else:
	print("reference given, performing reference based assembly")
	subprocess.call(["bwa", "index", args.reference], stderr=log)
	samfile = open("alignment.sam", "w")

#alignment step
	if not pair:
		print("performing single end alingment")
		subprocess.call(["bwa", "mem", args.reference, args.fastq1], stdout=samfile, stderr=log)
	else:
		print("performing paired end alignment")
		subprocess.call(["bwa", "mem", args.reference, args.fastq1, args.fastq2], stdout=samfile, stderr=log)

#convert to bam and sort
bamfile = open("alignment.bam", "w")
subprocess.call(["samtools", "view", "-S", "-b", "alignment.sam"], stdout=bamfile, stderr=log)
subprocess.call("samtools sort alignment.bam alignment_sorted", shell=True, stderr=log)
subprocess.call(["samtools", "index", "alignment_sorted.bam"], stderr=log) 
"""

#pileup reads at positions
fastqfile = open("consensus_alignment.fq", "w")
p = subprocess.Popen("samtools mpileup -uf " + args.reference + " alignment_sorted.bam | bcftools view -cg - | vcfutils.pl vcf2fq", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
fastqfile.write(p.stdout.read())
log.write(p.stderr.read())

#subprocess.call("rm " + args.reference + ".*", shell=True)
