import pysam
import re
import sys
import os
from tqdm import tqdm
import pandas as pd

def create_folder(folder_path):
	if not os.path.exists(folder_path):
		os.makedirs(folder_path)
		print(f"Folder '{folder_path}' created successfully.")
	else:
		print(f"Folder '{folder_path}' already exists.")


def create_folder_for_file(file_path):
	folder_path = os.path.dirname(file_path)
	if not os.path.exists(folder_path):
		os.makedirs(folder_path)
		print(f"Folder '{folder_path}' created successfully.")
	else:
		print(f"Folder '{folder_path}' already exists.")


def count_polyA_reads(input_bam_file):
	# Input BAM file path
	tmp_bam_file = pysam.AlignmentFile(input_bam_file, "rb")
	polyA_positions = []
	for read in tqdm(tmp_bam_file.fetch()):
		if read.is_reverse:
			if read.cigartuples[0][0] == pysam.CSOFT_CLIP and read.cigartuples[0][1] >= 6:
				clip_len=read.cigartuples[0][1]
				sub_seq=read.query_sequence[(clip_len-10):clip_len]
				retain_seq=read.query_sequence[clip_len:(clip_len+10)]
				if re.search(r".{0,2}?(T{3,}[^T]{0,2})*T{6,}([^T]{0,2}T{3,})*", sub_seq):
					if re.search(r".{0,2}?(T{3,}[^T]{0,2})*T{6,}([^T]{0,2}T{3,})*", retain_seq):
						polyA_positions.append((read.reference_name, "-", "IPE", str(read.reference_start+1)))
						# print(read.reference_name, "-", "IPE", str(read.reference_start))
					else:
						polyA_positions.append((read.reference_name, "-", "PAS", str(read.reference_start+1)))
						# print('\t'.join([read.reference_name, "-", "PAS", str(read.reference_start)]))
		else:
			if read.cigartuples[-1][0] == pysam.CSOFT_CLIP and read.cigartuples[-1][1] >= 6:
				clip_len=len(read.query_sequence)-read.cigartuples[-1][1]
				sub_seq=read.query_sequence[clip_len:(clip_len+10)]
				retain_seq=read.query_sequence[(clip_len-10):clip_len]
				if re.search(r"(A{3,}[^A]{0,2})*A{6,}([^A]{0,2}A{3,})*.{0,2}?", sub_seq):
					if re.search(r"(A{3,}[^A]{0,2})*A{6,}([^A]{0,2}A{3,})*.{0,2}?", retain_seq):
						polyA_positions.append((read.reference_name, "+", "IPE", str(read.reference_end)))
						# print('\t'.join([read.reference_name, "+", "IPE", str(read.reference_end)]))
					else:
						polyA_positions.append((read.reference_name, "+", "PAS", str(read.reference_end)))
						# print('\t'.join([read.reference_name, "+", "PAS", str(read.reference_end)]))
	tmp_bam_file.close()
	return(polyA_positions)

# Check the execution mode
if __name__ == "__main__":
	# Get the input file name
	if len(sys.argv) < 3:
		print("Please provide the study and batch path as an argument")
		sys.exit(1)


polyA_count = count_polyA_reads(sys.argv[2])

# Output the results to a table
out_file = "{0}/site_report/{1}_polyA_positions.csv".format(sys.argv[3], sys.argv[1])
create_folder_for_file(out_file)
df = pd.DataFrame(polyA_count, columns=['chr', 'strand', 'type', 'coord'])
count_df = df.groupby(['chr', 'strand', 'type', 'coord']).size().reset_index(name='count')
count_df.to_csv(out_file, index=False)
