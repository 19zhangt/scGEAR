
# https://www.biostars.org/p/461377/

from cgat import GTF
from cgatcore import iotools
import sys

input_file = iotools.open_file(sys.argv[1])

gtf_enteries = GTF.iterator(input_file)
transcripts = GTF.transcript_iterator(gtf_enteries)

for transcript in transcripts:
    exons = [e for e in transcript if e.feature == "exon"]
    if len(exons) == 0:
        continue
    if transcript[0].strand == "-":
        exons = sorted(exons, key=lambda x: -1*x.end)
    else:
        exons = sorted(exons, key=lambda x: x.start)
    last_exon = exons[-1]
    print (str(last_exon))
