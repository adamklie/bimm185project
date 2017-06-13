#!/usr/bin/python
import MySQLdb
from collections import defaultdict

def mediaQuery(strain, filename):
	strain = strain.replace("\'", "")

	db = MySQLdb.connect(host="localhost",    # your host, usually localhost
			     user="root",         # your username
			     passwd="doubleateam",  # your password
			     db="mediaDB")        # name of the data base


	# you must create a Cursor object. It will let
	#  you execute all the queries you need
	cur = db.cursor()


	# Use all the SQL you like
	cur.execute("SELECT strainID, Genus, Species, Strain, typeID, has_model, biomassID FROM organisms WHERE Strain = '{}';".format(strain))

	# If strain not found
	if not cur.rowcount:
		return(0)	
		

	# Holds all results
	strainID = ""
	genus = ""
	species = ""
	strain = ""
	hasModel = ""

	medIDList = []
	sourceIDList = []
	growth_RateList = []
	growth_UnitsList = []
	pHList = []
	temperature_CList = []
	additional_NotesList = []
	contributor_idList = []

	mediaNamesList = []
	isMinimalList = []

	sourceAuthorList = []
	sourceYearList = []
	sourceJournalList = []
	sourceLinkList = []
	sourcePubMedIDList = []


	# Save all results for the strain query
	if cur.rowcount > 1:
		print("WARNING: {} number of hits found for strain {}.\nIgnoring all except last hit".format(cur.rowcount, strain))

	for row in cur.fetchall():
		strainID = row[0]
		genus = row[1]
		species = row[2]
		strain = row[3]
		hasModel = row[5]


	# use strainID to query medID, Growth_Rate, Growth_Units, pH, Temperature_C, Additional_Notes, approved(?) from growth_data
	cur.execute("SELECT growthID, medID, sourceID, Growth_Rate, Growth_Units, pH, Temperature_C, Additional_Notes, contributor_id FROM growth_data WHERE strainID = {};".format((strainID)))


	# Save all results from the growth medium query
	for row in cur.fetchall():
		medIDList.append(row[1])
		sourceIDList.append(row[2])
		growth_RateList.append(row[3])
		growth_UnitsList.append(row[4])
		pHList.append(row[5])
		temperature_CList.append(row[6])
		additional_NotesList.append(row[7])
		contributor_idList.append(row[8])

	# use medID to query Media_name, is_Minimal from media_names 
	for id in medIDList:
		cur.execute("SELECT Media_name, Is_minimal FROM media_names WHERE medID = {};".format(id))

		for row in cur.fetchall():
			mediaNamesList.append(row[0])
			isMinimalList.append(row[1])
			

	# use medID to query compID from media_compounds
	compoundIDDict = defaultdict(list)
	compoundAmountDict = defaultdict(list)
	compoundNameDict = defaultdict(list)
	for id in medIDList:
		cur.execute("SELECT compID, Amount_mM FROM media_compounds WHERE medID = {}".format(id))
		
		for row in cur.fetchall():
			compoundIDDict[id].append(row[0])		
			compoundAmountDict[id].append(row[1])

		# use the compoundIDs for each mediaID to query the compound name
		for compID in compoundIDDict[id]:
			cur.execute("SELECT name FROM compounds WHERE compID = {};".format(compID))

			for row in cur.fetchall():
				compoundNameDict[id].append(row[0])


	# Get source info
	for id in sourceIDList:
		cur.execute("SELECT First_Author, Year, Link, pubmed_id, Journal FROM sources WHERE sourceID = {};".format(id))

		for row in cur.fetchall():
			sourceAuthorList.append(row[0])
			sourceYearList.append(row[1])
			sourceLinkList.append(row[2])
			sourcePubMedIDList.append(row[3])
			sourceJournalList.append(row[4])


	# Print results
	filename.write("Organism: {} {} {}\n".format(genus, species, strain))
	filename.write("Has Model? {}\n\n".format(hasModel))
	filename.write("{} Growth Conditions Found:\n\n".format(len(medIDList)))
	i = 0
	for id in medIDList:
		filename.write("Medium: {}\n".format(mediaNamesList[i]))
		if isMinimalList[i]:
			filename.write("Medium is minimal? Yes\n")
		else:
			filename.write("Medium is minimal? No\n")
		filename.write("In-depth Medium Details:\n")
		for compName, compAmount in zip(compoundNameDict[id], compoundAmountDict[id]):
			filename.write("\tCompound: {}\tAmount: {} mM\n".format(compName, compAmount))
		filename.write("Growth Data on Medium:\n")
		filename.write("\tGrowth Rate: {} {}\n".format(growth_RateList[i], growth_UnitsList[i]))
		filename.write("\tpH: {}\n".format(pHList[i]))
		filename.write("\tTemperature: {} C\n".format(temperature_CList[i]))
		filename.write("\tAdditional Notes: {}\n".format(additional_NotesList[i]))
		filename.write("Source Info:\n")
		filename.write("\tAuthor: {} Et al, {}\n".format(sourceAuthorList[i], sourceYearList[i]))
		filename.write("\tJournal: {}\n".format(sourceJournalList[i]))
		filename.write("\tPubmed ID: {}\n".format(sourcePubMedIDList[i]))
		filename.write("\tLink: {}\n".format(sourceLinkList[i]))
		filename.write('\n')

		i += 1

	db.close()
	return(1)
