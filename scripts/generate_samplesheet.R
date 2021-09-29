
## Samplesheet and FASTA file generation
# Generate one samplesheet for the run
# Generate one set of FASTA files (i5 and i6) per sample index within the run
# FASTA files named {group_name}_i5.fa and {group_name}_i7.fa
# Optionally reverse complement sequences in the sample file

args <- commandArgs(trailingOnly = TRUE)
samplepath <- args[1]
samplename <- args[2]
rc <- args[3]

samples <- read.table(samplepath, sep = "\t", header = TRUE)

# barcodes can be variable length
# can have multiple sample barcodes

# # check if all sample barcodes the same length
# bc1 <- sapply(samples$i5_index, nchar)
# bc2 <- sapply(samples$i7_index, nchar)

if (rc == "True") {
  message("Reverse complementing barcodes")
  revcomp <- function(x) {
    return(as.character(Biostrings::reverseComplement(Biostrings::DNAString(x))))
  }
  samples$i5_index <- unname(sapply(samples$i5_index, revcomp))
  samples$i7_index <- unname(sapply(samples$i7_index, revcomp))
}

# create sample sheet containing sample index
header <- paste0("[Header]
Experiment Name,", samplename, "\n",
"Date,", Sys.Date(), "\n",
"Module,GenerateFASTQ - 2.0.0
Workflow,GenerateFASTQ
Library Prep Kit,Custom
Chemistry,Amplicon
[Data]
Sample_ID,Description,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project")

outdir <- paste0(samplename, "_barcodes")
if (!dir.exists(paths = outdir)) {
  dir.create(path = outdir)
}
samplesheet <- paste0(outdir, "/samplesheet.csv")
write(x = header, file = samplesheet, append = TRUE)

if (any(!is.na(x = samples$sample_index))) {
  # get unique groups
  grp.keep <- !duplicated(x = samples$sample_index)
  groups <- samples$sample_index[grp.keep]
  grn <- samples$group_name[grp.keep]
  
  # split by sample index
  for (i in seq_along(along.with = groups)) {
    group_split <- unlist(x = strsplit(x = groups[i], ","))
    for (j in group_split) {
      outstr <- paste(
        grn[i], grn[i], j, j,
        "","",
        sep = ","
      )
      write(x = outstr, file = samplesheet, append = TRUE)
    }
  }
} else {
  # no sample index, demux everything to unassigned
  outstr <- paste(
    "unassigned", "unassigned",
    "","","","",
    sep = ","
  )
  write(x = outstr, file = samplesheet, append = TRUE)
}

# create barcode fasta files
unique_groups <- unique(samples$group_name)
for (i in seq_along(along.with = unique_groups)) {
  # create barcode fasta files for cutadapt
  fa1 <- paste0(outdir, "/", unique_groups[i], "_i5.fa")
  fa2 <- paste0(outdir, "/", unique_groups[i], "_i7.fa")
  samples_use <- samples[samples$group_name == unique_groups[i], ]
  for (i in seq_len(length.out = nrow(x = samples_use))) {
    write(x = paste0(">", samples_use[i, "sample_name"]), file = fa1, append = TRUE)
    write(x = paste0(">", samples_use[i, "sample_name"]), file = fa2, append = TRUE)
    write(x = samples_use[i, "i5_index"], file = fa1, append = TRUE)
    write(x = samples_use[i, "i7_index"], file = fa2, append = TRUE)
  }
}
