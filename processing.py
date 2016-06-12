from __future__ import absolute_import
import csv
import sys
import os

# Reads in a tab-delimited file and returns a 2D array representation
def readInputFile(sharedFile):
	otuMap = []
	with open(sharedFile) as csvfile:
		sharedFileTSV = csv.reader(open(sharedFile), delimiter='\t', quotechar='|')
		otuMap = processOTUMap(sharedFileTSV)
	return otuMap

# Given a csv reader object, returns a 2D array representation
def processOTUMap(otuMapCSV):
	otuMap = []
	for o in otuMapCSV:
		if len(o) > 1:
			otuMap.append(o)
	return otuMap

# Given a 2D array representation, produces a tab-delimited file with the new name
def exportToFile(base, name):
	if sys.version_info[0] < 3:
		with open(name, 'wb') as csvfile:
			outputCSV = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			for row in base:
				outputCSV.writerow(row)
	else:
		with open(name, 'w', newline='') as csvfile:
			outputCSV = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			for row in base:
				outputCSV.writerow(row)