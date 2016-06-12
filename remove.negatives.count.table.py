from __future__ import absolute_import
from __future__ import print_function
import csv
import sys
import os
import processing as ps

COUNT_TABLE_SEQUENCE_COL = 0
COUNT_TABLE_SAMPLE_START_COL = 2

sampleIDCol = 0
diseaseCol = 1
batchCol = 2

if len(sys.argv) != 5:
	print("Usage: python " + sys.argv[0] + " [.count_table File] [Metadata File] [Negative Control Identifier] [Use Batch? - Y/N]")
	sys.exit(1)

print("\n")

sharedFile = sys.argv[1]
metadataFile = sys.argv[2]
negSamplePhrase = sys.argv[3]
useBatch = sys.argv[4]

isBatch = False
if useBatch.lower() == "y" or useBatch.lower() == "yes":
	isBatch = True
	print("Negative Control Sample Sequence Removal [By Batch]\n")
else:
	print("Negative Control Sample Sequence Removal\n")

sharedFile = sharedFile.replace('"', "")
countTable = ps.readInputFile(sharedFile)
metadata = ps.readInputFile(metadataFile)

# 
# Process
# 

newFilename = sharedFile.rsplit('.', 1)[0]
fileExt = sharedFile.rsplit('.', 1)[-1]
if isBatch:
	newFilename = newFilename + ".batch.negremoved." + fileExt
else:
	newFilename = newFilename + ".negremoved." + fileExt

accnosFile = sharedFile.rsplit('.', 1)[0]
if isBatch:
	accnosFile = accnosFile + ".batch.negremoved.accnos"
else:
	accnosFile = accnosFile + ".negremoved.accnos"


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
			isPos = False
			for negativeCol in negativeCols:
				if float(countTable[i][negativeCol]) > 0:
					isPos = True

			total = 0
			j = COUNT_TABLE_SAMPLE_START_COL
			while j < len(countTable[i]):
				if j not in negativeCols and countTable[i][j] != "":
					if isPos:
						countTable[i][j] = 0
					total += int(countTable[i][j])
				j += 1

			countTable[i][1] = total

			if total == 0:
				rowsToRemove[i] = 1
				numRemoved += 1

		if i % 200 == 0 and i > 0:
			print("Processed " + str(i) + " / " + str(len(countTable)))
		i += 1

	# Creates new count table excluding the negative samples columns and excluding the sequences with all zero counts
	removedSequences = []
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
		else:
			removedSequences.append([countTable[r][COUNT_TABLE_SEQUENCE_COL]])
		r += 1
		
	ps.exportToFile(newCountTable, newFilename)

	print("=================================================")
	print("Finished removing negative sequence counts")
	if numRemoved > 0:
		ps.exportToFile(removedSequences, accnosFile)
		print("Number of sequences removed : " + str(numRemoved))
		print("Accnos file at: " + accnosFile)
	print("Output file at: " + newFilename)




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
						if int(countTable[i][negativeBatchCol]) > 0:
							countTable[i][j] = 0
						newSum += int(countTable[i][j])
				j += 1
			countTable[i][1] = newSum

			if newSum == 0:
				rowsToRemove[i] = 1
				numRemoved += 1

		if i % 200 == 0 and i > 0:
			print("Processed " + str(i) + " / " + str(len(countTable)))
		i += 1


	# Creates new count table excluding the negative samples columns and excluding the sequences with all zero counts
	removedSequences = []
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
		else:
			removedSequences.append([countTable[r][COUNT_TABLE_SEQUENCE_COL]])
		r += 1
		
	ps.exportToFile(newCountTable, newFilename)

	print("=================================================")
	print("Finished removing negative sequence counts by batch")
	if numRemoved > 0:
		ps.exportToFile(removedSequences, accnosFile)
		print("Number of sequences removed : " + str(numRemoved))
		print("Accnos file at: " + accnosFile)
	print("Output file at: " + newFilename)