import csv
import sys
import os
import processing as ps

OTU_TABLE_SAMPLE_ID_COL = 1
OTU_TABLE_OTU_START_COL = 3

print "Negative Control Sample OTU Removal"
print "-----------------------------------"

sharedFile = raw_input('Enter the .shared file location (no spaces in filename): \n').strip()
sharedFile = str(sharedFile)
sharedFile = sharedFile.replace('"', "")

print "\n"

otuTableSharedTSV = csv.reader(open(sharedFile), delimiter='\t', quotechar='|')
base = ps.processOTUMap(otuTableSharedTSV)

metadataFile = raw_input('Enter the tab-delimited metadata file location (no spaces in filename): \n').strip()
# metadataFile = str(metadataFile)
# metadataFile = metadataFile.replace('"', "")
sampleIDMapTSV = csv.reader(open(metadataFile), delimiter='\t', quotechar='|')
metadata = ps.processOTUMap(sampleIDMapTSV)

print "\n"

sampleIDCol = raw_input('Which column (first column = 1) in the metadata file is the Sample ID? (eg. the column that corresponds to the Group column in the .shared file) \n')
sampleIDCol = int(sampleIDCol) - 1

print "\n"

diseaseCol = raw_input('Which column (first column = 1) in the metadata file is used to identify whether a sample is a negative control? \n')
diseaseCol = int(diseaseCol) - 1

print "\n"

negSamplePhrase = raw_input('In column ' + str(diseaseCol) + ' of the metadata file, which phrase is used to identify negative samples (eg. Negative) \n')

print "\n"

isBatchStr = raw_input('Are the negative controls using batches? [Yes/No] \n')
print "\n"

isBatch = False
batchCol = -1
if (isBatchStr.lower() == "yes"):
	isBatch = True
	batchCol = raw_input('Which column (first column = 1) in the metadata file is used to identify which batch a sample belongs to? \n')
	batchCol = int(batchCol) - 1

print "\n"


# 
# Process
# 

newFilename = sharedFile.rsplit('.', 1)[0]
fileExt = sharedFile.rsplit('.', 1)[-1]
if isBatch:
	newFilename = newFilename + ".batch.negremoved." + fileExt
else:
	newFilename = newFilename + ".negremoved." + fileExt

if isBatch:
	sampleIDToBatch = {}
	negBatchToSampleID = {}
	sampleIDToNegBatch = {}
	for row in metadata:
		if row[diseaseCol] == negSamplePhrase:
			negBatchToSampleID[row[batchCol]] = row[sampleIDCol]
			sampleIDToNegBatch[row[sampleIDCol]] = row[batchCol]
		sampleIDToBatch[row[sampleIDCol]] = row[batchCol]

	rowsToAvoid = {}
	negBatchToRow = {}
	i = 1
	while i < len(base):
		sampleID = base[i][OTU_TABLE_SAMPLE_ID_COL]
		if sampleID in sampleIDToNegBatch:
			negBatchToRow[sampleIDToBatch[sampleID]] = i
			rowsToAvoid[i] = 1
		i += 1

	c = OTU_TABLE_OTU_START_COL
	while c < len(base[0]):
		r = 1
		while r < len(base):
			sampleID = base[r][OTU_TABLE_SAMPLE_ID_COL]
			batch = sampleIDToBatch[sampleID]
			if r not in rowsToAvoid and batch in negBatchToRow:
				if float(base[negBatchToRow[batch]][c]) > 0:
					# The corresponding negative batch is positive for this OTU so set the sample to be 0
					base[r][c] = 0
			r += 1
		c = c + 1

	# Removes any zero OTUs
	numOTUsRemoved = 0
	colsToRemove = {}
	c = OTU_TABLE_OTU_START_COL
	while c < len(base[0]):
		allZeros = True
		r = 1
		while r < len(base):
			sampleID = base[r][OTU_TABLE_SAMPLE_ID_COL]
			if r not in rowsToAvoid:
				if float(base[r][c]) > 0:
					allZeros = False
			r += 1
		if allZeros:
			numOTUsRemoved += 1
			colsToRemove[c] = 1
		c = c + 1

	# Creates new OTU Table without negative rows and removed OTUs
	r = 0
	newBase = []
	while r < len(base):
		if r not in rowsToAvoid:
			newRow = []
			c = 0
			while c < len(base[r]):
				if c not in colsToRemove:
					newRow.append(base[r][c])
				c = c + 1
			newBase.append(newRow)
		r += 1

	ps.exportToCSV(newBase, newFilename)
	print "Removed " + str(numOTUsRemoved) + " OTUs"


else:
	negSampleIDs = {}
	for row in metadata:
		if row[diseaseCol] == negSamplePhrase:
			negSampleIDs[row[sampleIDCol]] = 1

	negSampleRows = []
	rowsToAvoid = {}
	i = 1
	while i < len(base):
		sampleID = base[i][OTU_TABLE_SAMPLE_ID_COL]
		if sampleID in negSampleIDs:
			rowsToAvoid[i] = 1
			negSampleRows.append(i)
		i += 1

	# Finds OTUs to remove
	numOTUsRemoved = 0
	colsToRemove = {}
	c = OTU_TABLE_OTU_START_COL
	while c < len(base[0]):
		toRemove = False
		r = 1
		while r < len(base):
			sampleID = base[r][OTU_TABLE_SAMPLE_ID_COL]
			for neg in negSampleRows:
				if float(base[neg][c]) > 0:
					# The corresponding negative sample is positive for this OTU so it will be removed
					toRemove = True
			r += 1
		if toRemove:
			colsToRemove[c] = 1
			numOTUsRemoved += 1
		c = c + 1

	# Creates new OTU Table without negative rows and removed OTUs
	r = 0
	newBase = []
	while r < len(base):
		if r not in rowsToAvoid:
			newRow = []
			c = 0
			while c < len(base[r]):
				if c not in colsToRemove:
					newRow.append(base[r][c])
				c = c + 1
			newBase.append(newRow)
		r += 1

	ps.exportToCSV(newBase, newFilename)
	print "Removed " + str(numOTUsRemoved) + " OTUs"