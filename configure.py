import os

db_type = "nucl" 
resistance_db  = "makeblastdb -in /home/adamklie/BIMM185/project/new_tool/bimm185project/databases/all_resistance.fsa -dbtype " + db_type
virulence_db = "makeblastdb -in /home/adamklie/BIMM185/project/new_tool/bimm185project/databases/all_virulence.fsa -dbtype " + db_type
serotype_db = "makeblastdb -in /home/adamklie/BIMM185/project/new_tool/bimm185project/databases/all_serotype.fsa -dbtype " + db_type
plasmid_db = "makeblastdb -in /home/adamklie/BIMM185/project/new_tool/bimm185project/databases/all_plasmid.fsa -dbtype " + db_type
os.system(resistance_db)
os.system(virulence_db)
os.system(serotype_db)
os.system(plasmid_db)
