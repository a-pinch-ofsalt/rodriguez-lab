dirs <- c("A", "B", "C", "D", "E", "G")
filenames <- c("H_F", "H_T", "K_F", "K_T", "L_F", "L_T")

for (dir in dirs) {
  for (fname in filenames) {
    print(fname)
    filepath <- file.path(dir, paste0(fname, ".tsv"))
    df <- read.table(filepath, header = TRUE, sep = "\t")
    lengths <- nchar(df$cdr3_aa)
    print(length(lengths))
    out_file <- file.path(paste0(dir, "_", fname, "_cdr3_aa_length_hist.png"))
    png(out_file)
    hist(lengths, breaks = 35, main = paste(dir), xlab = "Length")
    dev.off()
  }
}
