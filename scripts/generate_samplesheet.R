
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
samplesheet <- paste0(outdir, "/samplesheet.csv")
if (!dir.exists(paths = outdir)) {
  dir.create(path = outdir)
} else {
  if (file.exists(samplesheet)) {
    file.remove(samplesheet)
  }
}
# demux everything to unassigned
write(x = header, file = samplesheet, append = TRUE)
outstr <- paste(
  "unassigned", "unassigned",
  "","","","",
  sep = ","
)
write(x = outstr, file = samplesheet, append = TRUE)

# create barcode fasta files
unique_groups <- unique(samples$group_name)

si_file <- paste0(outdir, "/sampleindex.fa")
i5_file <- paste0(outdir, "/tn5_i5.fa")
i7_file <- paste0(outdir, "/tn5_i7.fa")

# sample index
for (i in seq_along(along.with = unique_groups)) {
  # create barcode fasta files for cutadapt
  samples_use <- samples[samples$group_name == unique_groups[i], ]
  
  # write sample index fasta
  if (!length(x = unique(x = samples_use$sample_index))) {
    stop("Same samples listed with different sample barcodes. Check the configuration files.")
  }
  si_use <- samples_use$sample_index[[1]]
  si_use <- unlist(x = strsplit(x = si_use, split = ","))
  for (j in seq_along(along.with = si_use)) {
    write(x = paste0(">", unique_groups[i]), file = si_file, append = TRUE)
    write(x = si_use[[j]], file = si_file, append = TRUE)
  }
}

# tn5 barcodes
for (i in seq_len(length.out = nrow(x = samples))) {
  samples_use <- samples[i, ]
  write(x = paste0(">", samples_use$sample_name), file = i5_file, append = TRUE)
  write(x = samples_use$i5_index, file = i5_file, append = TRUE)
  write(x = paste0(">", samples_use$sample_name), file = i7_file, append = TRUE)
  write(x = samples_use$i7_index, file = i7_file, append = TRUE)
}