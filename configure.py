import os

n_type = "nucl" 
p_type = "prot"
resistance_db  = "makeblastdb -in ./databases/all_resistance.fsa -dbtype " + n_type
virulence_db = "makeblastdb -in ./databases/all_virulence.fsa -dbtype " + n_type
serotype_db = "makeblastdb -in ./databases/all_serotype.fsa -dbtype " + n_type
plasmid_db = "makeblastdb -in ./databases/all_plasmid.fsa -dbtype " + n_type
#rnammer_db = "makeblastdb -in ./databases/rnammer.fsa -dbtype " + n_type
db_16S = "makeblastdb -in ./databases/rnammer_16S.fsa -dbtype " + n_type

os.system(resistance_db)
os.system(virulence_db)
os.system(serotype_db)
os.system(plasmid_db)
#os.system(rnammer_db)
os.system(db_16S)
