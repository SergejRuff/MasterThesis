import subprocess

gff3_file = "top_contig.fasta.gff3"
edited_gff3 = "edited_top_contig.fasta.txt"

subprocess.run(['grep', 'ID=', gff3_file], stdout=open(edited_gff3, 'w'))

import re
import pandas as pd

# Open the text file and read its content
with open('edited_top_contig.fasta.txt', 'r') as file:
    lines = file.readlines()

# Regular expressions to match the required fields
sequence_id_re = re.compile(r'^(SRR\d+_cap3_Contig-\d+)\s')
orf_re = re.compile(r'\s+(\d+)\s+(\d+)\s+')
signature_desc_re = re.compile(r'signature_desc=([^;]+);')
interpro_accession_re = re.compile(r'Dbxref="InterPro:(IPR\d+)"')

# Dictionary to store sequence lengths for each sequence_id
sequence_lengths = {}

# List to store extracted data
extracted_data = []

# First pass: extract sequence lengths from polypeptide entries
for line in lines:
    # Extract the sequence_id (1st column)
    sequence_match = sequence_id_re.match(line)
    if sequence_match:
        sequence_id = sequence_match.group(1)

    # Check if the line contains a polypeptide entry and extract the ORF length
    if "polypeptide" in line:
        orf_match = orf_re.search(line)
        if orf_match:
            orf_start = int(orf_match.group(1))
            orf_end = int(orf_match.group(2))
            sequence_length = orf_end - orf_start + 1  # Calculate sequence length
            sequence_lengths[sequence_id] = sequence_length

# Second pass: extract other relevant information and add sequence length
for line in lines:
    # Extract the sequence_id (1st column)
    sequence_match = sequence_id_re.match(line)
    if sequence_match:
        sequence_id = sequence_match.group(1)

    # Extract the ORF start and end positions (even for non-polypeptide lines)
    orf_match = orf_re.search(line)
    if orf_match:
        orf_start = int(orf_match.group(1))
        orf_end = int(orf_match.group(2))
    else:
        orf_start = None
        orf_end = None

    # Extract signature description
    signature_desc_match = signature_desc_re.search(line)
    if signature_desc_match:
        signature_description = signature_desc_match.group(1)
    else:
        signature_description = None

    # Extract interpro accession
    interpro_accession_match = interpro_accession_re.search(line)
    if interpro_accession_match:
        interpro_accession = interpro_accession_match.group(1)
        interpro_description = signature_description  # Assuming description is the same as signature_desc
    else:
        interpro_accession = None
        interpro_description = None

    # Add the sequence length (if available) to the corresponding rows based on sequence_id
    sequence_length = sequence_lengths.get(sequence_id, None)  # Retrieve the sequence length based on sequence_id

    # If it's a row with signature or interpro data, append it
    if signature_description or interpro_accession:
        extracted_data.append({
            'sequence_id': sequence_id,
            'orf_start': orf_start,
            'orf_end': orf_end,
            'sequence_length': sequence_length,  # Add the sequence length from polypeptide entry
            'signature_description': signature_description,
            'interpro_accession': interpro_accession,
            'interpro_description': interpro_description
        })

# Create a pandas DataFrame from the extracted data
df = pd.DataFrame(extracted_data)
df.to_csv('interproscan_data.tsv', sep='\t', index=False)
df