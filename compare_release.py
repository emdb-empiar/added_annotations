import pathlib
import argparse
import sys
import os

latest_data_path = sys.argv[1]
previous_data_path = sys.argv[2]

def count_unique(filename):
	entries = set()
	with open(filename, 'r') as fr:
		for line in fr:
			emdb_id = line.split("\t")[0]
			entries.add(emdb_id)
	return len(entries)-1 #Subtract the index

message = "Resource    Number of entries (current)    Number of entries (previous)    Difference\n"

for latest_file_path in pathlib.Path(latest_data_path).glob('*.log'):
	filename = latest_file_path.parts[-1]
	resource = filename.split("_")[-1].split(".")[0].upper()
	count_latest = count_unique(str(latest_file_path))
	previous_file_path = previous_data_path + filename.replace('log','tab')
	if pathlib.Path(previous_file_path).is_file():
		count_previous = count_unique(str(previous_file_path))
		diff = count_latest-count_previous
		diff_str = f"+{diff}" if diff > 0 else str(diff)
		message += f"{resource}    {count_latest}    {count_previous}    {diff_str}\n"

os.system(f'mail -s "Summary of weekly annotations workflow" pdb_em@ebi.ac.uk <<< "{message}"')