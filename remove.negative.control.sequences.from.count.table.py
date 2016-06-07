import csv
import sys
import os
import processing as ps

OTU_TABLE_SAMPLE_ID_COL = 1
OTU_TABLE_OTU_START_COL = 3

if len(sys.argv) != 6 and len(sys.argv) != 7:
	print "Usage: python " + sys.argv[0] + " [.count_table File] [Metadata File] [Sample ID Col] [Sample Type Col] [Negative Control Sample Type]  [Batch Col - Optional]"
	sys.exit(1)
sharedFile = sys.argv[1]
metadataFile = sys.argv[2]
sampleIDCol = int(sys.argv[3]) - 1
diseaseCol = int(sys.argv[4]) - 1
negSamplePhrase = sys.argv[5]

isBatch = False
batchCol = -1
if len(sys.argv) == 7:
	isBatch = True
	batchCol = int(sys.argv[6]) - 1


print "Negative Control Sample Sequence Removal"
print "-----------------------------------"

sharedFile = sharedFile.replace('"', "")
otuTableSharedTSV = csv.reader(open(sharedFile), delimiter='\t', quotechar='|')
countTable = ps.processOTUMap(otuTableSharedTSV)

sampleIDMapTSV = csv.reader(open(metadataFile), delimiter='\t', quotechar='|')
metadata = ps.processOTUMap(sampleIDMapTSV)


# 
# Process
# 

newFilename = sharedFile.rsplit('.', 1)[0]
fileExt = sharedFile.rsplit('.', 1)[-1]
if isBatch:
	newFilename = newFilename + ".batch.negremoved." + fileExt
else:
	newFilename = newFilename + ".negremoved." + fileExt

if not isBatch:
	negativeIDs = {}
	for row in metadata:
		if row[diseaseCol] == negSamplePhrase:
			negativeIDs[row[sampleIDCol]] = 1

	numRemoved = 0
	negativeCols = []
	newCountTable = []
	i = 0
	while i < len(countTable):
		if i == 0:
			newCountTable.append(countTable[i])
			j = 2
			while j < len(countTable[i]):
				if countTable[i][j] in negativeIDs:
					negativeCols.append(j)
				j += 1
		else:
			total = 0
			for negativeCol in negativeCols:
				total += int(countTable[i][negativeCol])

			if total == 0:
				newCountTable.append(countTable[i])
			else:
				numRemoved += 1
		i += 1

	ps.exportToCSV(newCountTable, newFilename)
	print "Processed file at " + newFilename
	print "Number of sequences removed : " + str(numRemoved)

else:
	negativeBatchToSampleID = {}
	sampleIDToDisease = {}
	sampleIDToBatch = {}
	for row in metadata:
		sampleID = row[sampleIDCol]
		sampleIDToDisease[sampleID] = row[diseaseCol]
		sampleIDToBatch[sampleID] = row[batchCol]
		if row[diseaseCol] == negSamplePhrase:
			negativeBatchToSampleID[row[batchCol]] = sampleID

	sampleIDToCol = {}
	colToSampleID = {}
	i = 0
	while i < len(countTable):
		if i == 0:
			j = 2
			while j < len(countTable[i]):
				colToSampleID[j] = countTable[i][j]
				sampleIDToCol[countTable[i][j]] = j
				j += 1
		else:
			# For each representative sequence, get the negative control values per batch
			j = 2
			newSum = 0
			while j < len(countTable[i]):
				sampleID = colToSampleID[j]
				batch = sampleIDToBatch[sampleID]
				disease = sampleIDToDisease[sampleID]
				if disease != negSamplePhrase and batch in negativeBatchToSampleID:
					negativeBatchCol = sampleIDToCol[negativeBatchToSampleID[batch]]
					if float(countTable[i][negativeBatchCol]) > 0:
						countTable[i][j] = 0
				newSum += int(countTable[i][j])
				j += 1
			countTable[i][1] = newSum

		if i % 200 == 0:
			print "Processed " + str(i) + " / " + str(len(countTable))
		i += 1


	newCountTable = []
	numRemoved = 0
	i = 0
	while i < len(countTable):
		if i == 0:
			newCountTable.append(countTable[i])
		else:
			total = 0
			j = 2
			while j < len(countTable[i]):
				total += int(countTable[i][j])
				j += 1
			if total > 0:
				newCountTable.append(countTable[i])
			else:
				numRemoved += 1
		i += 1

	ps.exportToCSV(newCountTable, newFilename)
	print "Processed file at " + newFilename
	print "Number of sequences removed : " + str(numRemoved)
