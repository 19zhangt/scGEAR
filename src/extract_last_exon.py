import gzip
import sys

def parse_attributes(attr_str):
    attributes = {}
    for part in attr_str.split(';'):
        part = part.strip()
        if not part:
            continue
        # Split into key and value, considering possible spaces in values (which should be quoted)
        key_value = part.split(' ', 1)
        if len(key_value) != 2:
            continue
        key = key_value[0]
        value = key_value[1].strip('"')  # Remove surrounding quotes
        attributes[key] = value
    return attributes

def extract_last_exons(gtf_path, output_path):
    transcripts = {}  # Key: transcript_id, Value: {max_exon_num: int, line: str}

    # Open GTF file (handles gzipped files if ending with .gz)
    open_func = open if not gtf_path.endswith('.gz') else gzip.open
    with open_func(gtf_path, 'rt') as f:  # 'rt' mode for reading text
        for line in f:
            if line.startswith('#'):
                continue  # Skip comment lines
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue  # Malformed line
            feature_type = fields[2]
            if feature_type != 'exon':
                continue  # Only process exon features
            attributes = parse_attributes(fields[8])
            transcript_id = attributes.get('transcript_id')
            exon_number = attributes.get('exon_number')
            
            if not transcript_id or not exon_number:
                continue  # Skip if missing required attributes
            
            try:
                exon_number = int(exon_number)
            except ValueError:
                continue  # Invalid exon_number format
            
            # Update the transcript entry
            current = transcripts.get(transcript_id)
            if not current or exon_number > current['max_exon_num']:
                transcripts[transcript_id] = {
                    'max_exon_num': exon_number,
                    'line': line
                }

    # Write the last exons to the output file
    with open(output_path, 'w') as out_f:
        for data in transcripts.values():
            out_f.write(data['line'])

# Example usage
extract_last_exons(sys.argv[1], sys.argv[2])
