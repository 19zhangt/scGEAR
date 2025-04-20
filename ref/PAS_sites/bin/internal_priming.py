import subprocess
import sys

try:
    import Bio
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython"])

from Bio import SeqIO
import re
import sys

def merge_regions(regions):
    if not regions:
        return []
    sorted_regions = sorted(regions, key=lambda x: x[0])
    merged = [sorted_regions[0]]
    for current in sorted_regions[1:]:
        last = merged[-1]
        if current[0] <= last[1] + 1:
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    return merged

def sliding_window_analysis(fasta_file, window_size=10, step=1):
    """Sliding window detection of A enriched regions (≥7 A or ≥6 consecutive A)"""
    results = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        risk_windows = []

        # Sliding window traversal (step length = 1 to ensure full coverage)
        for i in range(len(seq) - 9):
            window = seq[i:i+10]
            a_count = window.count('A')
            consecutive_a = max( (len(match) for match in re.findall(r'A+', window)), default=0 )
            if a_count >= 7 or consecutive_a >= 6:
                risk_windows.append( (i+1, i+10) )

        merged_regions = merge_regions(risk_windows)
        results[record.id] = {
            "regions": merged_regions,
        }
    return results

# Example
results = sliding_window_analysis(sys.argv[1])

with open(sys.argv[2], "w") as res_info:
    for seq_id, data in results.items():
        for start, end in data['regions']:
            res_info.write(f"{seq_id}\t{start}\t{end}\tA_rich\n")
