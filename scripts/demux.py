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
        all_bc[str(i.seq[:maxlen])] = i.id
    return(all_bc)


def get_entry(f):
    return([f.readline(), f.readline(), f.readline(), f.readline()])


def extract_barcodes(sequence, adapter, full_adapter_len, bc2_len):
    # return first 9 bp as tn5_bc, remove next 15 bp (adaptor), then return remaining as sample index
    tn5_bc = sequence[:8]
    trimmed = sequence[9:]
    si = trimmed.split(adapter)
    if len(si) == 1:
        si = trimmed[full_adapter_len:][:-1] # remove trailing \n
    else:
        si = si[1][:-1]
    return((tn5_bc, si))



tn5_i5_barcodes = get_barcodes(args.tn5_i5, maxlen=8) # using first 8 bases of barcode
tn5_i7_barcodes = get_barcodes(args.tn5_i7, maxlen=8)
sample_barcodes = get_barcodes(args.sampleindex, maxlen=8)

r1 = open_fastq(f=args.read1)
r2 = open_fastq(f=args.read2)
r3 = open_fastq(f=args.read3)
r4 = open_fastq(f=args.read4)

i5_index_valid = list(set(tn5_i5_barcodes.values()))
i7_index_valid = list(set(tn5_i7_barcodes.values()))
sample_index_valid = list(set(sample_barcodes.values()))
i5_index_valid.append("unknown")
i7_index_valid.append("unknown")
sample_index_valid.append("unknown")

# create dictionary with file handles
# key = barcode combination
# value = file handle
outf = dict()
outpath = Path(args.output)
if not outpath.exists():
    os.mkdir(outpath)

for s5 in i5_index_valid:
    for n7 in i7_index_valid:
        for j in sample_index_valid:
            fpath = outpath / j
            if not fpath.exists():
                os.mkdir(fpath)
            bc =  s5 + "-" + n7
            # outputting uncompressed fastq is ~10x faster
            fname_1 = open(fpath / (bc + ".R1.fastq"), "w+")
            fname_2 = open(fpath / (bc + ".R2.fastq"), "w+")
            outf[bc + "-" + j] = (fname_1, fname_2)

x = 0
while True:
    r1_entry = get_entry(r1)
    r2_entry = get_entry(r2)
    r3_entry = get_entry(r3)
    r4_entry = get_entry(r4)
    
    if r1_entry[0] == '':
        break
    
    # get barcodes
    i7_barcodes = extract_barcodes(sequence=r2_entry[1], adapter="CGAGAC", full_adapter_len=15, bc2_len=8)
    i5_barcodes = extract_barcodes(sequence=r3_entry[1], adapter="CGACGA", full_adapter_len=14, bc2_len=16)
    
    # match to list of valid barcodes
    # if no match, classify as unknown
    if i7_barcodes[0] in tn5_i7_barcodes:
        i7_index_name = tn5_i7_barcodes[i7_barcodes[0]]
    else:
        i7_index_name = "unknown"

    if i5_barcodes[0] in tn5_i5_barcodes:
        i5_index_name = tn5_i5_barcodes[i5_barcodes[0]]
    else:
        i5_index_name = "unknown"
    
    if i7_barcodes[1] in sample_barcodes:
        sample_index_name = sample_barcodes[i7_barcodes[1]]
    else:
        sample_index_name = "unknown"

    cell_barcode = i5_barcodes[1]
    
    # add cell barcode to r1 and r2 genomic
    r1_entry[0] = "@" + cell_barcode + ":" + r1_entry[0][1:]
    r4_entry[0] = "@" + cell_barcode + ":" + r4_entry[0][1:]
    
    # write to correct output files based on barcode combination
    outfile = i5_index_name + "-" + i7_index_name + "-" + sample_index_name
    r1_outf = outf[outfile][0]
    r2_outf = outf[outfile][1]
    
    r1_outf.write(r1_entry[0])
    r1_outf.write(r1_entry[1])
    r1_outf.write(r1_entry[2])
    r1_outf.write(r1_entry[3])

    r2_outf.write(r4_entry[0])
    r2_outf.write(r4_entry[1])
    r2_outf.write(r4_entry[2])
    r2_outf.write(r4_entry[3])
    x += 1
    if x % 1e6 == 0:
        print("Processed " + str(int(x/1e6)) + " million reads", file=sys.stderr, end="\r")
    
# close all files
for i in outf.keys():
    outf[i][0].close()
    outf[i][1].close()

# another option is to only split by sample index, but add the tn5 barcode combination to the read name also
# can then map one set of fastqs per sample index, and extract the fragment information from the read name
# would need to modify the fragment function to get the extra information
# in the end we want one fragment file per condition anyway, unless we want to modify the fragment file format to add metadata column (would need to change a lot)