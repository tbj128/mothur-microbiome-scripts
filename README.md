# mothur-microbiome-scripts
Useful scripts to modify data files as part of the mothur SOP

Usage Examples:

python subtract.negatives.count.table.py Data/test.count_table Data/test.metadata.txt 1 2 Negative
python subtract.negatives.count.table.py Data/test.count_table Data/test.metadata.txt 1 2 Negative 7

python subtract.negatives.shared.py Data/test.shared Data/test.metadata.txt 1 2 Negative
python subtract.negatives.shared.py Data/test.shared Data/test.metadata.txt 1 2 Negative 7


python remove.negatives.count.table.py Data/test.count_table Data/test.metadata.txt 1 2 Negative
python remove.negatives.count.table.py Data/test.count_table Data/test.metadata.txt 1 2 Negative 7

python remove.negatives.shared.py Data/test.shared Data/test.metadata.txt 1 2 Negative
python remove.negatives.shared.py Data/test.shared Data/test.metadata.txt 1 2 Negative 7