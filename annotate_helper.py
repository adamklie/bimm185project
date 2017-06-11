import os
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW

def blastn(query_filepath, db_filepath):
	db_out = db_filepath.split("/")[-1].split(".")[0]
	blastOut = db_out + ".xml"	#xml file for parsing
	blastOut2 = db_out + ".txt"	#txt file for viewing, remove when confident in functionality
	queryBlastTXT = 'blastn -query ' + query_filepath + ' -db ' + db_filepath + ' > ' + blastOut2
	queryBlastXML = 'blastn -query ' + query_filepath + ' -db ' + db_filepath + ' -outfmt 5 ' + ' > ' + blastOut
	os.system(queryBlastXML)
	os.system(queryBlastTXT)
	return

def blastp(query_filepath, db_filepath):
        db_out = db_filepath.split("/")[-1].split(".")[0]
        blastOut = db_out + ".xml"      #xml file for parsing
        blastOut2 = db_out + ".txt"     #txt file for viewing, remove when confident in functionality
        queryBlastTXT = 'blastp -query ' + query_filepath + ' -db ' + db_filepath + ' > ' + blastOut2
        queryBlastXML = 'blastp -query ' + query_filepath + ' -db ' + db_filepath + ' -outfmt 5 ' + ' > ' + blastOut
        os.system(queryBlastXML)
        os.system(queryBlastTXT)
        return

def parseBlastXML(blast_out_filepath):
	result_handle = open(blast_out_filepath)
	blast_records = NCBIXML.parse(result_handle)
	E_VALUE_THRESH = 1
	#Loop through each protein query results
	for blast_record in blast_records:
		#Loop through the hits associated with particular sequence
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				#Hit must have e-value < threshold to be considered
				if hsp.expect < E_VALUE_THRESH:
					title = alignment.title	#title of hit
					print(title)
					print(blast_record.query)
			break
	return

def get16S(rnammer_out_filepath):
	"""
	result_handle = open(rnammer_out_filepath)
	blast_records = NCBIXML.parse(result_handle)
	max_score = float("-inf")
	best_hit_rna = ""
	for blast_record in blast_records:
		#database_name = 
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				title = alignment.title
				rna_type = title.split("/")[1].rstrip()
				if rna_type == "molecule=16s_rRNA":
					if hsp.score > max_score:
						max_score = hsp.score
						best_hit_rna = title
			break
	rna_title = best_hit_rna.split(" ")[1]
	result_handle.close()
     	result_handle = open(rnammer_out_filepath)
	blast_records = NCBIXML.parse(result_handle)
	database_name = list(blast_records)[0].database
	result_handle.close()
	print(database_name)
	rna_title = "rRNA_CU928145_4597992-4599521_DIR+"
	rna_database = open(database_name, "r")
	rna_seq = ""
	
	if rna_title != "":
		is_seq = False
		for line in rna_database: 
			if line[0] == ">":
				if is_seq:
					break
				elif rna_title in line:
					is_seq = True
			if is_seq:
				rna_seq += line

	rna_seq = ">rRNA_CU928145_4597992-4599521_DIR+ /molecule=16s_rRNA /score=1958.7\nAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAACAGCTTGCTGTTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTTGTTGGTGGGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCT"

	result = NCBIWWW.qblast("blastn", "refseq_genomic", rna_seq)
	rna_record = NCBIXML.read(result)
	for alignment in rna_record.alignments:
		for hsp in alignment.hsps:
			title = alignment.title
			#parse_title = title.split("|")[4].split(",")[0]
			print(title)
			#print(parse_title)
				
"""
	f = open("temp.txt", "r")
	titles = []
	for line in f:
		titles.append(line.lstrip())
	strain_ids = []
	for title in titles:
		parse_title = " ".join(title.replace(" complete genome", "").replace(" chromosome", "").replace(" whole genome shotgun sequence", "").replace(" DNA", "").split("|")[4].split(",")[0].split(" ")[3:]).rstrip("\n")
		if len(parse_title.split(" ")) == 1:
			strain_ids.append(parse_title)
		else:
			double_parse = parse_title.split(" ")
			if "str." in double_parse:
				ind = double_parse.index("str.")
				strain_ids.append(double_parse[ind+1])
			elif "strain" in double_parse:
				ind = double_parse.index("strain")
				strain_ids.append(double_parse[ind+1])
			else:
				strain_ids.append(" ".join(double_parse))
	
	for ids in strain_ids:
		ids = ids.lstrip(" ").rstrip(" ")
		ids = ids.replace("'", "").replace("\"", "")
	no_duplicates_ids = []
	ids_set = set()
	for ids in strain_ids:
		if ids in ids_set: 
			continue
		no_duplicates_ids.append(ids)
		ids_set.add(ids)
	for ids in no_duplicates_ids:
		print(ids)
#testing
blastn("./velvet_out/contigs.fa", "./databases/rnammer.fsa")
get16S("rnammer.xml")
