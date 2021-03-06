from __future__ import absolute_import
from __future__ import print_function
import csv
import sys
import os
import processing as ps

OTU_TABLE_SAMPLE_ID_COL = 1
OTU_TABLE_OTU_START_COL = 3

sampleIDCol = 0
diseaseCol = 1
batchCol = 2

if len(sys.argv) != 5:
	print("Usage: python " + sys.argv[0] + " [.shared File] [Metadata File] [Negative Control Sample Type] [Use Batch? - Y/N]")
	sys.exit(1)

print("\n")

sharedFile = sys.argv[1]
metadataFile = sys.argv[2]
negSamplePhrase = sys.argv[3]
useBatch = sys.argv[4]

isBatch = False
if useBatch.lower() == "y" or useBatch.lower() == "yes":
	isBatch = True
	print("Negative Control Sample OTU Subtraction [By Batch]\n")
else:
	print("Negative Control Sample OTU Subtraction\n")

sharedFile = sharedFile.replace('"', "")
otuTableSharedTSV = ps.readInputFile(sharedFile)
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
			numOTUs = 0
			newRow = []
			c = 0
			while c < len(base[r]):
				if c not in colsToRemove:
					newRow.append(base[r][c])
					if c >= OTU_TABLE_OTU_START_COL:
						numOTUs += 1
				c = c + 1

			if r > 0:
				newRow[2] = numOTUs
			newBase.append(newRow)
		r += 1

	ps.exportToFile(newBase, newFilename)

	print("=================================================")
	print("Finished!")
	print("Removed " + str(numOTUsRemoved) + " OTUs")
	print("Output file at: " + newFilename)

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
			numOTUs = 0
			c = 0
			while c < len(base[r]):
				if c not in colsToRemove:
					newRow.append(base[r][c])
					if c >= OTU_TABLE_OTU_START_COL:
						numOTUs += 1
				c = c + 1

			if r > 0:
				newRow[2] = numOTUs
			newBase.append(newRow)
		r += 1

	ps.exportToFile(newBase, newFilename)

	print("=================================================")
	print("Finished!")
	print("Removed " + str(numOTUsRemoved) + " OTUs")
	print("Output file at: " + newFilename)