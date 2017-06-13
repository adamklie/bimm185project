#! /usr/bin/env python
import annotate_helper
import subprocess
import mediaQuery

# parse user arguments 
import argparse
import sys
parser = argparse.ArgumentParser()
parser.add_argument("-f1", "--fastq1", required=True, help="forward fastq reads from bacterial sample", type=str)
parser.add_argument("-f2", "--fastq2", required=False, help="reverse fastq reads from bacterial sample    ", type=str)
parser.add_argument("-rt", "--read_type", required=True, help="type of read library input, s=single end, p=paired end, m=mate pair", type=str)
parser.add_argument("-i", "--insert_size", required=False, help="insert size for paired end assembly", type=str)
parser.add_argument("-r", "--reference", help="optional reference fasta to perform assembly against", required=False, type=str)
parser.add_argument("-v", "--verbosity", action="store_true", help="increase output verbosity")
pair = False
args = parser.parse_args()

#make output directory for sequences
subprocess.call("mkdir SeqPro_out", shell=True)

#setting up log file
if args.verbosity:
	log = sys.stderr
else:	
	log = open("./SeqPro_out/log.txt", "w")

#check read type
if args.read_type in "pm":
	pair =  True
	if args.fastq2 is None:
		print("read pair not specified, exiting...")
		sys.exit()

#make output directory for assembly
subprocess.call("mkdir assembly_out", shell=True)

#velvet assembly with no reference
if args.reference is None: 
	print("No reference given, performing de novo assembly")
	if args.read_type == "s":
		print("Performing single end assembly with Velvet")
		subprocess.call(["velveth", "assembly_out", "19", "-fastq", "-short", args.fastq1, args.fastq2], stdout=log)
		subprocess.call(["velvetg", "assembly_out"], stdout=log)	
	elif args.read_type == "p":
		if args.insert_size is None:
			print("Must specify insert size for paired end assembly, exiting...")
			sys.exit()
		print("Performing paired end assembly with Velvet")
		subprocess.call(["velveth", "assembly_out", "19", "-fastq", "-shortPaired", args.fastq1, args.fastq2], stdout=log)
		subprocess.call(["velvetg", "assembly_out", "-ins_length", args.insert_size, "-exp_cov", "auto",], stdout=log)  
	elif args.read_type == "m":
		print("Don't have functionality built in yet, exiting...")
		sys.exit()
	else:
		print("Unrecognizable pair type, exiting...")
		sys.exit()
	subprocess.call("mv assembly_out/contigs.fa .", shell=True)

#reference based assembly with reference given
else:
	print("Reference given, performing reference based assembly")
	subprocess.call("mkdir align_out", shell=True)
	subprocess.call(["bwa", "index", args.reference], stderr=log)
	samfile = open("alignment.sam", "w")

	#alignment step
	if not pair:
		print("Performing single end alingment with bwa")
		subprocess.call(["bwa", "mem", args.reference, args.fastq1], stdout=samfile, stderr=log)
	else:
		print("Performing paired end alignment with bwa")
		subprocess.call(["bwa", "mem", args.reference, args.fastq1, args.fastq2], stdout=samfile, stderr=log)

	samfile.close()

	#convert to bam and sort
	bamfile = open("alignment.bam", "w")
	subprocess.call(["samtools", "view", "-S", "-b", "alignment.sam"], stdout=bamfile, stderr=log)
	bamfile.close()
	subprocess.call("samtools sort alignment.bam alignment_sorted", shell=True, stderr=log)
	subprocess.call(["samtools", "index", "alignment_sorted.bam"], stderr=log) 
	
	#pileup reads at positions and generate fasta, need intermediate fastq with this pipeline
	fastqfile = open("contigs.fq", "w")
	p = subprocess.Popen("samtools mpileup -uf " + args.reference + " alignment_sorted.bam | bcftools view -cg - | vcfutils.pl vcf2fq", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
	fastqfile.write(p.stdout.read())
	log.write(p.stderr.read())
	fastqfile.close()
	
	subprocess.call(["seqret", "-osformat", "fasta", "contigs.fq", "-out2", "contigs.fa"], stderr=log) 
	subprocess.call("mv alignment* align_out", shell=True)
	subprocess.call("mv contigs.fq assembly_out", shell=True)

#annotate the assemblies using Prokka
print("Generating annotations with Prokka")
subprocess.call(["prokka", "contigs.fa", "--rnammer"], stderr=log)
subprocess.call("cp PROKKA*/PROKKA*.ffn .", shell=True)
subprocess.call("mv PROKKA*.ffn annotated_seqs.ffn", shell=True)

#finding the closest relative
print("Finding the most closely related strain")
annotate_helper.blastn("./contigs.fa", "./databases/rnammer_16S.fsa")
rna_16S = annotate_helper.get16S("rnammer_16S.xml")
if rna_16S is not None:
	rna_16S_ids = rna_16S[0]

	#outputting the closest relative
	closely_related_file = open("./SeqPro_out/closest_relative.txt", "w")
	closely_related_file.write("Closest related bacteria: ")
	closely_related = rna_16S[1].rstrip().lstrip()
	closely_related_file.write(closely_related)
	print("Most closely related bacterial strain: " + closely_related)
	print("")
	closely_related_file.close

	#finding the best media
	print("Detecting suggested media")
	print("")
	mediafile = open("./SeqPro_out/suggested_media.txt", "w")
	mediafile = open("./output/suggested_media.txt", "w")
	foundMedia = False
	for strain_id in rna_16S_ids:
		media = mediaQuery.mediaQuery(strain_id, mediafile)
		if media == 1:
			foundMedia = True
			break	
	if foundMedia == False:
		mediafile.write("No media data for closest relatives found in MediaDB.\nNOTE: This has been caused by occasional incompatible naming schemes of closest relatives (e.g. K12 vs. K-12) and inconsistent results returned from the 16S rRNA database.\nWe are working to resolve this issue.")
	mediafile.close()
else:
	print("Could not find a closest relative, will not output media")

#finding antibiotic resistance genes
resistancefile = open("./SeqPro_out/resistance_sequences.tsv", "w")
resistancefile.write("Hit\tE-Value\tPROKKA ID\tPROKKA Annotation\n")
print("Detecting antibiotic resistance")
annotate_helper.blastn("./annotated_seqs.ffn", "./databases/all_resistance.fsa")
annotate_helper.parseBlastXML("all_resistance.xml", resistancefile)
resistancefile.close()

#finding virulence genes
virulencefile = open("./SeqPro_out/virulence_sequences.tsv", "w")
virulencefile.write("Hit\tE-Value\tPROKKA ID\tPROKKA Annotation\n")
print("Detecting virulence factors")
annotate_helper.blastn("./annotated_seqs.ffn", "./databases/all_virulence.fsa")
annotate_helper.parseBlastXML("all_virulence.xml", virulencefile)
virulencefile.close()

#finding plasmid sequences
plasmidfile = open("./SeqPro_out/plasmid_sequences.tsv", "w")
plasmidfile.write("Hit\tE-Value\tContig\tPROKKA Annotation\n")
print("Detecting plasmid sequences")
annotate_helper.blastn("./contigs.fa", "./databases/all_plasmid.fsa")
annotate_helper.parseBlastXML("all_plasmid.xml", plasmidfile)
plasmidfile.close()

#finding bacterial serotype
serotypefile = open("./SeqPro_out/serotype_sequences.tsv", "w")
serotypefile.write("Hit\tE-Value\tPROKKA ID\tPROKKA Annotation\n")
print("Detecting serotype")
annotate_helper.blastn("./annotated_seqs.ffn", "./databases/all_serotype.fsa")
annotate_helper.parseBlastXML("all_serotype.xml", serotypefile)
serotypefile.close()

#move output files into proper directories
subprocess.call("mv contigs.fa assembly_out", shell=True)
subprocess.call("mv PROKKA* Prokka_out", shell=True)

# delete the ish
if args.reference != None:
	subprocess.call("rm " + args.reference + ".*", shell=True)
subprocess.call("rm *.xml", shell=True)
subprocess.call("rm annotated_seqs.ffn", shell=True)
subprocess.call("rm *.pyc", shell=True)
