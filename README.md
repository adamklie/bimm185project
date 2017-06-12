# bimm185project
Microbial Profile Generation

# Raw Ecoli X. reads used to test ProSeq can be found here; SRR292678 - paired end reads, insert size of 470 bp
foward reads: https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292678sub_S1_L001_R1_001.fastq.gz  
reverse reads: https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292678sub_S1_L001_R2_001.fastq.gz

# Dependencies
SeqPro requires:  
-Velvet  
-BWA-mem  
-BLAST  
-MySQL loaded with MediaDB  
-CGE's ResFinder, VirulenceFinder, PlasmidFinder, SerotypeFinder Databases  
-CBA's 16S rRNA Database  

# Setting Up SeqPro
python configure.py


# Running SeqPro
usage: SeqPro.py [-h] -f1 FASTQ1 [-f2 FASTQ2] -rt READ_TYPE [-r REFERENCE] [-v]  

optional arguments:  
  -h, --help            show this help message and exit  
  -f1 FASTQ1, --fastq1 FASTQ1  
                        forward fastq reads from bacterial sample  
  -f2 FASTQ2, --fastq2 FASTQ2  
                        reverse fastq reads from bacterial sample  
  -rt READ_TYPE, --read_type READ_TYPE  
                        type of read library input, s=single end, p=paired  
                        end, m=mate pair  
  -r REFERENCE, --reference REFERENCE  
                        optional reference fasta to perform assembly against  
  -v, --verbosity       increase output verbosity  
