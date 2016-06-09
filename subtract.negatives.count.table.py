import csv
import sys
import os
import processing as ps

COUNT_TABLE_SEQUENCE_COL = 0
COUNT_TABLE_SAMPLE_START_COL = 2

if len(sys.argv) != 6 and len(sys.argv) != 7:
	print "Usage: python " + sys.argv[0] + " [.count_table File] [Metadata File] [Sample ID Col] [Sample Type Col] [Negative Control Sample Type]  [Batch Col - Optional]"
	sys.exit(1)

print "\n"

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
	print "Negative Control Sample Sequence Subtraction [By Batch]\n"
else:
	print "Negative Control Sample Sequence Subtraction\n"
	
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
	newFilename = newFilename + ".batch.negsubtracted." + fileExt
else:
	newFilename = newFilename + ".negsubtracted." + fileExt


if not isBatch:
	negativeIDs = {}
	for row in metadata:
		if row[diseaseCol] == negSamplePhrase:
			negativeIDs[row[sampleIDCol]] = 1

	numRemoved = 0
	negativeCols = []
	rowsToRemove = {}
	i = 0
	while i < len(countTable):
		if i == 0:
			j = COUNT_TABLE_SAMPLE_START_COL
			while j < len(countTable[i]):
				if countTable[i][j] in negativeIDs:
					negativeCols.append(j)
				j += 1
		else:
			# Take the average of the negative control samples
			negTotal = 0
			for negativeCol in negativeCols:
				negTotal += float(countTable[i][negativeCol])
			negAvg = negTotal / float(len(negativeCols))

			# Subtract the "averaged" negative control sample values
			total = 0
			j = COUNT_TABLE_SAMPLE_START_COL
			while j < len(countTable[i]):
				if j not in negativeCols and countTable[i][j] != "":
					countTable[i][j] = float(countTable[i][j]) - negAvg
					if countTable[i][j] < 0:
						countTable[i][j] = 0
					total += countTable[i][j]
				j += 1

			countTable[i][1] = total

			if total == 0:
				rowsToRemove[i] = 1
				numRemoved += 1

		if i % 200 == 0 and i > 0:
			print "Processed " + str(i) + " / " + str(len(countTable))
		i += 1

	# Creates new count table excluding the negative samples columns and excluding the sequences with all zero counts
	r = 0
	newCountTable = []
	while r < len(countTable):
		if r not in rowsToRemove:
			newRow = []
			c = 0
			while c < len(countTable[r]):
				if c not in negativeCols and countTable[r][c] != "":
					newRow.append(countTable[r][c])
				c = c + 1
			newCountTable.append(newRow)
		r += 1
		
	ps.exportToCSV(newCountTable, newFilename)

	print "================================================="
	print "Finished subtracting negative sequence counts"
	if numRemoved > 0:
		print "Number of sequences removed : " + str(numRemoved)
	print "Output file at: " + newFilename






else:
	negativeBatchToSampleID = {}
	sampleIDToDisease = {}
	sampleIDToBatch = {}
	negativeIDs = {}
	for row in metadata:
		sampleID = row[sampleIDCol]
		sampleIDToDisease[sampleID] = row[diseaseCol]
		sampleIDToBatch[sampleID] = row[batchCol]
		if row[diseaseCol] == negSamplePhrase:
			negativeBatchToSampleID[row[batchCol]] = sampleID
			negativeIDs[sampleID] = 1

	sampleIDToCol = {}
	colToSampleID = {}
	negativeCols = []
	rowsToRemove = {}
	numRemoved = 0
	i = 0
	while i < len(countTable):
		if i == 0:
			j = COUNT_TABLE_SAMPLE_START_COL
			while j < len(countTable[i]):
				if countTable[i][j] != "":
					colToSampleID[j] = countTable[i][j]
					sampleIDToCol[countTable[i][j]] = j
					if countTable[i][j] in negativeIDs:
						negativeCols.append(j)
				j += 1
		else:
			# For each representative sequence, get the negative control values per batch
			j = COUNT_TABLE_SAMPLE_START_COL
			newSum = 0
			while j < len(countTable[i]):
				if countTable[i][j] != "":
					sampleID = colToSampleID[j]
					batch = sampleIDToBatch[sampleID]
					if j not in negativeCols and batch in negativeBatchToSampleID:
						negativeBatchCol = sampleIDToCol[negativeBatchToSampleID[batch]]
						countTable[i][j] = int(countTable[i][j]) - int(countTable[i][negativeBatchCol]) 
						if countTable[i][j] < 0:
							countTable[i][j] = 0
						newSum += int(countTable[i][j])
				j += 1
			countTable[i][1] = newSum

			if newSum == 0:
				rowsToRemove[i] = 1
				numRemoved += 1

		if i % 200 == 0 and i > 0:
			print "Processed " + str(i) + " / " + str(len(countTable))
		i += 1


	# Creates new count table excluding the negative samples columns and excluding the sequences with all zero counts
	r = 0
	newCountTable = []
	while r < len(countTable):
		if r not in rowsToRemove:
			newRow = []
			c = 0
			while c < len(countTable[r]):
				if c not in negativeCols:
					newRow.append(countTable[r][c])
				c = c + 1
			newCountTable.append(newRow)
		r += 1
		
	ps.exportToCSV(newCountTable, newFilename)

	print "================================================="
	print "Finished subtracting negative sequence counts by batch"
	if numRemoved > 0:
		print "Number of sequences removed : " + str(numRemoved)
	print "Output file at: " + newFilename