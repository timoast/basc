
# Generate one set of FASTA files (i5 and i6) per sample index within the run
# FASTA files named {group_name}_i5.fa and {group_name}_i7.fa
# Optionally reverse complement sequences in the sample file

args <- commandArgs(trailingOnly = TRUE)
samplepath <- args[1]
samplename <- args[2]
rc_i5 <- args[3]
rc_i7 <- args[4]

samples <- read.table(samplepath, header = TRUE)
outdir <- paste0(samplename, "_barcodes")

if (rc_i5 == "True") {
  message("Reverse complementing i5 barcodes")
  revcomp <- function(x) {
    return(as.character(Biostrings::reverseComplement(Biostrings::DNAString(x))))
  }
  samples$i5_index <- unname(sapply(samples$i5_index, revcomp))
}
if (rc_i7 == "True") {
  message("Reverse complementing i7 barcodes")
  revcomp <- function(x) {
    return(as.character(Biostrings::reverseComplement(Biostrings::DNAString(x))))
  }
  samples$i7_index <- unname(sapply(samples$i7_index, revcomp))
}

# create barcode fasta files
unique_groups <- unique(samples$sample_name)

si_file <- paste0(outdir, "/sampleindex.fa")
i5_file <- paste0(outdir, "/tn5_i5.fa")
i7_file <- paste0(outdir, "/tn5_i7.fa")

if (!dir.exists(paths = outdir)) {
  dir.create(path = outdir)
} else {
  if (file.exists(si_file)) {
    file.remove(si_file)
  }
  if (file.exists(i5_file)) {
    file.remove(i5_file)
  }
  if (file.exists(i7_file)) {
    file.remove(i7_file)
  }
}

# sample index
for (i in seq_along(along.with = unique_groups)) {
  samples_use <- samples[samples$sample_name == unique_groups[i], ]
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
# write mark and well combination
for (i in seq_len(length.out = nrow(x = samples))) {
  samples_use <- samples[i, ]
  comb <- paste0(">", samples_use$sample_name, "+", samples_use$well, "+", samples_use$mark)
  write(x = comb, file = i5_file, append = TRUE)
  write(x = samples_use$i5_index, file = i5_file, append = TRUE)
  write(x = comb, file = i7_file, append = TRUE)
  write(x = samples_use$i7_index, file = i7_file, append = TRUE)
}
