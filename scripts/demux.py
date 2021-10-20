import pysam
from Bio import SeqIO
import os
import gzip
from pathlib import Path
from argparse import ArgumentParser
import sys


# parse args
parser = ArgumentParser(description='Read demultiplexer')
parser.add_argument('--read1', help='Path to read1')
parser.add_argument('--read2', help='Path to read2')
parser.add_argument('--read3', help='Path to read3')
parser.add_argument('--read4', help='Path to read4')
parser.add_argument('--sampleindex', help='Path to sample index FASTA file')
parser.add_argument('--tn5_i5', help='Path to Tn5 i5 index FASTA file')
parser.add_argument('--tn5_i7', help='Path to Tn5 i7 index FASTA file')
parser.add_argument('--min_barcode_len', help='Minimum length of cell barcode', default=16, type=int)
parser.add_argument('--output', help='Path to output directory')
args = parser.parse_args()


def open_fastq(f):
    if os.path.isfile(f):
        if f.endswith(".gz"):
            f_open = gzip.open(f, "rt") # rt read text mode (decodes binary gzip)
        else:
            f_open = open(f, "r")
        return(f_open)
    else:
        raise Exception("File not found")


def get_barcodes(inpath, maxlen=8):
    all_bc = dict()
    for i in SeqIO.parse(open(inpath),'fasta'):
        values = i.id.split("+")
        samplename = values[0]
        well = values[1]
        mark = values[2]
        barcode = str(i.seq[:maxlen])
        if samplename in all_bc.keys():
            all_bc[samplename][barcode] = (mark, well)
        else:
            all_bc[samplename] = {barcode: (mark, well)}
    return(all_bc)


def read_sample_barcodes(inpath, maxlen=8):
    all_bc = dict()
    for i in SeqIO.parse(open(inpath),'fasta'):
        all_bc[str(i.seq[:maxlen])] = i.id
    return(all_bc)


def get_entry(f):
    return([f.readline(), f.readline(), f.readline(), f.readline()])


def extract_barcodes(sequence, adapter, full_adapter_len, bc1_len, bc2_len):
    # return first 9 bp as tn5_bc, remove next 15 bp (adaptor), then return remaining as sample index
    tn5_bc = sequence[:bc1_len]
    trimmed = sequence[bc1_len+1:]
    idx_match = trimmed.find(adapter)
    if idx_match == -1:
        # not found
        si = trimmed[full_adapter_len:][:-1] # remove trailing \n
    else:
        si = trimmed[idx_match+len(adapter):][:-1] # remove trailing \n
        if len(si) > bc2_len:
            si = si[:bc2_len]
    return((tn5_bc, si))


tn5_i5_barcodes = get_barcodes(args.tn5_i5, maxlen=8) # using first 8 bases of barcode
tn5_i7_barcodes = get_barcodes(args.tn5_i7, maxlen=8)
sample_barcodes = read_sample_barcodes(args.sampleindex, maxlen=8)

r1 = open_fastq(f=args.read1)
r2 = open_fastq(f=args.read2)
r3 = open_fastq(f=args.read3)
r4 = open_fastq(f=args.read4)

# tn5_i5_barcodes = get_barcodes("./NTT_barcodes/tn5_i5.fa", maxlen=8)
# tn5_i7_barcodes = get_barcodes("./NTT_barcodes/tn5_i7.fa", maxlen=8)
# sample_barcodes = read_sample_barcodes("./NTT_barcodes/sampleindex.fa", maxlen=8)

# r1 = open_fastq(f="/home/stuartt/data/210924_NB552173_0145_AHVTNNBGXH/fastq/unassigned_R1.fastq.gz")
# r2 = open_fastq("/home/stuartt/data/210924_NB552173_0145_AHVTNNBGXH/fastq/unassigned_R2.fastq.gz")
# r3 = open_fastq("/home/stuartt/data/210924_NB552173_0145_AHVTNNBGXH/fastq/unassigned_R3.fastq.gz")
# r4 = open_fastq("/home/stuartt/data/210924_NB552173_0145_AHVTNNBGXH/fastq/unassigned_R4.fastq.gz")

# create dictionary with file handles
# key = barcode combination
# value = file handle
outf = dict()
# outpath = Path("./test")
outpath = Path(args.output)
if not outpath.exists():
    os.mkdir(outpath)

for sampleindex in tn5_i5_barcodes.keys():
    fpath = outpath / sampleindex
    if not fpath.exists():
        os.mkdir(fpath)
    unique_i5_marks = list(set([x[1][0] for x in tn5_i5_barcodes[sampleindex].items()]))
    unique_i7_marks = list(set([x[1][0] for x in tn5_i7_barcodes[sampleindex].items()]))
    unique_i5_marks.append('unknown')
    unique_i7_marks.append('unknown')
    for i5 in unique_i5_marks:
        for i7 in unique_i7_marks:
            bc =  i5 + "-" + i7
            # outputting uncompressed fastq is ~10x faster
            fname_1 = open(fpath / (bc + ".R1.fastq"), "w+")
            fname_2 = open(fpath / (bc + ".R2.fastq"), "w+")
            outf[bc + "-" + sampleindex] = (fname_1, fname_2)

# unknown sample index
unknown_path = outpath / "unknown"
if not unknown_path.exists():
    os.mkdir(unknown_path)
outf["unknown"] = (open(unknown_path / ("unknown.R1.fastq"), "w+"),
                   open(unknown_path / ("unknown.R2.fastq"), "w+"))

x = 0
while True:
    r1_entry = get_entry(r1)
    r2_entry = get_entry(r2)
    r3_entry = get_entry(r3)
    r4_entry = get_entry(r4)

    if r1_entry[0] == '':
        break

    # get barcodes
    i7_barcodes = extract_barcodes(sequence=r2_entry[1], adapter="CGAGAC", full_adapter_len=15, bc1_len=8, bc2_len=8)
    i5_barcodes = extract_barcodes(sequence=r3_entry[1], adapter="CGACGA", full_adapter_len=14, bc1_len=8, bc2_len=16)

    if i7_barcodes[1] in sample_barcodes:
        sample_index_name = sample_barcodes[i7_barcodes[1]]
        # look up mark names
        if i7_barcodes[0] in tn5_i7_barcodes[sample_index_name]:
            i7_mark, well = tn5_i7_barcodes[sample_index_name][i7_barcodes[0]]
        else:
            i7_mark, well = "unknown", "unknown"
        if i5_barcodes[0] in tn5_i5_barcodes[sample_index_name]:
            i5_mark, well = tn5_i5_barcodes[sample_index_name][i5_barcodes[0]]
        else:
            i5_mark, well = "unknown", "unknown"
    else:
        sample_index_name = "unknown"
        i7_mark, well = "unknown", "unknown"
        i5_mark, well = "unknown", "unknown"

    cell_barcode = i5_barcodes[1]
    if len(cell_barcode) >= args.min_barcode_len:
        # add barcodes to r1 and r2 genomic
        # need to include well information here
        bc_combination = "@" + cell_barcode + "-" + well + ":" + i5_barcodes[0] + "+" + i7_barcodes[0] + "+"
        r1_entry[0] = bc_combination + r1_entry[0][1:]
        r4_entry[0] = bc_combination + r4_entry[0][1:]
        
        if sample_index_name == "unknown":
            r1_outf = outf['unknown'][0]
            r2_outf = outf['unknown'][1]
        else:
            # write to correct output files based on barcode combination
            outfile = i5_mark + "-" + i7_mark + "-" + sample_index_name
            r1_outf = outf[outfile][0]
            r2_outf = outf[outfile][1]

        r1_outf.write("".join(r1_entry))
        r2_outf.write("".join(r4_entry))
    x += 1
    if x % 1e6 == 0:
        print("Processed " + str(int(x/1e6)) + " million reads", file=sys.stderr, end="\r")

# close all files
for i in outf.keys():
    outf[i][0].close()
    outf[i][1].close()
