import os
from Bio.Blast import NCBIXML

def blastAgainstDatabase(query_filepath, db_filepath):
	db_out = db_filepath.split("/")[-1].split(".")[0]
	blastOut = db_out + ".xml"	#xml file for parsing
	blastOut2 = db_out + ".txt"	#txt file for viewing
	queryBlastTXT = 'blastn -query ' + query_filepath + ' -db ' + db_filepath + ' > ' + blastOut2
	queryBlastXML = 'blastn -query ' + query_filepath + ' -db ' + db_filepath + ' -outfmt 5 ' + ' > ' + blastOut
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
	return

#testing
blastAgainstDatabase("./tests/test_ecolix_annotation.ffn", "./databases/all_virulence.fsa")
parseBlastXML("all_virulence.xml")
