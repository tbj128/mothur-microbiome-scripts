# mothur Negative Control Removal/Subtraction Tool
Useful scripts to modify data files as part of the mothur SOP

Usage:

Remove all sequences that occur in the negative controls from the count table
python remove.negatives.count.table.py [count file] [metadata file] [negative sample identifier] [is batch? y/n]

Subtracts negative control sequence reads from each sample in the count table
python subtract.negatives.count.table.py [count file] [metadata file] [negative sample identifier] [is batch? y/n]

Remove all sequences that occur in the negative controls from the shared file
python remove.negatives.shared.py [shared file] [metadata file] [negative sample identifier] [is batch? y/n]

Subtracts negative control sequence reads from each sample in the shared table
python subtract.negatives.shared.py [shared file] [metadata file] [negative sample identifier] [is batch? y/n]

Metadata File:
Tab-delimited file such that first column is sample ID, second column is sample type (how negative sample or treatment sample is identified), third optional column is the batch number (if the samples were processed in batches)

Negative Sample Identifier:
The single word phrase used in the second column of the metadata file that identifies a sample as being a negative control

Is Batch:
y/n - describes whether sequences should be removed or subtracted according to the corresponding negative control done as part of the same batch (may not be applicable to all experiments)