import csv
import sys
import os

def processOTUMap(otuMapCSV):
	otuMap = []
	for o in otuMapCSV:
		if len(o) > 1:
			otuMap.append(o)
	return otuMap

def exportToCSV(base, name):
	outputCSV = csv.writer(open(name, 'wb'), delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
	for row in base:
		outputCSV.writerow(row)
