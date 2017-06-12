import os

n_type = "nucl" 
p_type = "prot"
resistance_db  = "makeblastdb -in /home/adamklie/BIMM185/project/new_tool/bimm185project/databases/all_resistance.fsa -dbtype " + n_type
virulence_db = "makeblastdb -in /home/adamklie/BIMM185/project/new_tool/bimm185project/databases/all_virulence.fsa -dbtype " + n_type
serotype_db = "makeblastdb -in /home/adamklie/BIMM185/project/new_tool/bimm185project/databases/all_serotype.fsa -dbtype " + n_type
plasmid_db = "makeblastdb -in /home/adamklie/BIMM185/project/new_tool/bimm185project/databases/all_plasmid.fsa -dbtype " + n_type
rnammer_db = "makeblastdb -in /home/adamklie/BIMM185/project/new_tool/bimm185project/databases/rnammer.fsa -dbtype " + n_type

os.system(resistance_db)
os.system(virulence_db)
os.system(serotype_db)
os.system(plasmid_db)
os.system(rnammer_db)
