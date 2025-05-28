dirs <- c("A", "B", "C", "D", "E", "G")
filenames <- c("H_F", "H_T", "K_F", "K_T", "L_F", "L_T")
# 
# for (dir in dirs) {
#   for (fname in filenames) {
#     print(fname)
#     filepath <- file.path(dir, paste0(fname, ".tsv"))
#     df <- read.table(filepath, header = TRUE, sep = "\t")
#     lengths <- nchar(df$cdr3_aa)
#     print(length(lengths))
#     out_file <- file.path(paste0(dir, "_", fname, "_cdr3_aa_length_hist.png"))
#     png(out_file)
#     hist(lengths, breaks = 35, main = paste(dir), xlab = "Length")
#     dev.off()
#   }
# }
# 
# 
# all_repertoires <- function(func) {
#   dirs <- c("A", "B", "C", "D", "E", "G")
#   filenames <- c("H_F", "H_T", "K_F", "K_T", "L_F", "L_T")
#   total <- list()
#   for (dir in dirs) {
#     for (fname in filenames) {
#       filepath <- file.path(dir, paste0(fname, ".tsv"))
#       repertoire <- fread(filepath)  # Faster + safer than read.table
#       label <- paste(dir, fname, sep = "_")
#       total[[label]] <- func(repertoire)
#     }
#   }
#   
#   df <- rbindlist(total, fill = TRUE)
#   df[is.na(df)] <- 0  # Fill any remaining NAs
#   df$sample <- names(total)  # Add the sample names as a column
#   setcolorder(df, c("sample", setdiff(names(df), "sample")))  # Reorder
#   return(df)
#   
#   df <- rbindlist(total, fill = TRUE)
#   filenames <- as.vector(t(outer(c("A", "B", "C", "D", "E", "G"), c("H_F", "H_T", "K_F", "K_T", "L_F", "L_T"), paste, sep = "_")))
#   rownames(df) <- filenames
#   df <- df[-8]
#   y <- total
#   
#   # fwrite(df, as.character(substitute(func)) + ".csv", row.names=TRUE)
#   
#   
# }
  
  

library(data.table)

# Generic wrapper to apply a function over all repertoire files
all_repertoires <- function(func) {
  dirs <- c("A", "B", "C", "D", "E", "G")
  filenames <- c("H_F", "H_T", "K_F", "K_T", "L_F", "L_T")
  total <- list()
  
  for (dir in dirs) {
    for (fname in filenames) {
      print(fname)
      filepath <- file.path(dir, paste0(fname, ".tsv"))
      repertoire <- fread(filepath)
      label <- paste(dir, fname, sep = "_")
      total[[label]] <- func(repertoire)
    }
  }
  
  df <- rbindlist(total, fill = TRUE, idcol = "sample")
  return(df)
}


count_kmers <- function(repertoire) {
  kmers <- table(nchar(unlist(strsplit(repertoire$cdr3_aa, "\\*"))))
  counts <- kmers[c("3", "4", "5", "6","7", "8", "9")]
  print(counts)
}

all_repertoires(count_kmers)

library(data.table)

dirs <- LETTERS[1:7]
files <- c("H_F", "H_T", "K_F", "K_T", "L_F", "L_T")

# Prepare an empty list to collect results
results <- list()

for (dir in dirs) {
  for (file in files) {
    path <- file.path(dir, paste0(file, ".tsv"))
    if (!file.exists(path)) {
      warning(paste("File not found:", path))
      next
    }
    dat <- fread(path)
    mean_length <- mean(nchar(dat$cdr3_aa))
    # print(file)
    # 
    # mean_length <- mean(dat$cdr3_end[valid_rows] - dat$cdr3_start[valid_rows])
    # print(mean_length)
    
    results[[length(results) + 1]] <- data.table(
      Directory = dir,
      File = file,
      MeanCDR3Length = mean_length
    )
  }
}

# Combine all results into one data table
all_results <- rbindlist(results)

# Write to CSV
fwrite(all_results, "cdr3_means.csv")

print("Saved means to cdr3_means.csv")



filepath <- file.path(dir, paste0(fname, ".tsv"))
repertoire <- read.table(filepath, header = TRUE, sep = "\t")
kmers <- table(nchar(unlist(strsplit(repertoire$cdr3_aa, "\\*"))))
counts <- kmers[c("3", "4", "5", "6","7", "8", "9")]
counts_df <- as.data.frame(as.list(counts))
print(counts_df)

rowname <- paste(dir, fname, sep = "_")
counts_df$File <- rowname
all_counts_list[[rowname]] <- counts_df



aa_freq <- read.csv("cdr3_character_frequencies.csv")
rownames(aa_freq) <- aa_freq[,1]
aa_freq <- aa_freq[-1] #aa = amino acid

#heavy
fetal <- aa_freq[grepl("H_F", rownames(aa_freq)),]
transitional <- aa_freq[grepl("H_T", rownames(aa_freq)),]
fetal
transitional
wilcox.test(unlist(fetal), unlist(transitional), paired = TRUE)$p.value

#kappa
fetal <- aa_freq[grepl("K_F", rownames(aa_freq)),]
transitional <- aa_freq[grepl("K_T", rownames(aa_freq)),]
wilcox.test(unlist(fetal), unlist(transitional), paired = TRUE)$p.value

#lambda
fetal <- aa_freq[grepl("L_F", rownames(aa_freq)),]
transitional <- aa_freq[grepl("L_T", rownames(aa_freq)),]
wilcox.test(unlist(fetal), unlist(transitional), paired = TRUE)$p.value

#all
fetal <- aa_freq[grepl("_F", rownames(aa_freq)),]
transitional <- aa_freq[grepl("_T", rownames(aa_freq)),]
aas <- intersect(names(fetal), names(transitional))
p_values <- sapply(aas, function(aa){
  wilcox.test(unlist(fetal[aa], use.names=FALSE), unlist(transitional[aa], use.names=FALSE), paired=TRUE)$p.value
})
sorted_p_values <- p_values[order(p_values, decreasing=TRUE)]
write.table(sorted_p_values, sep = "\t", quote = FALSE)


#Identify Kmers 3 - 9 and their frequency for every repertoire. Detained if kmer frequency differ between stages. So like determine whether the lengths change too?
fetal <- aa_freq[grepl("_F", rownames(aa_freq)),]

all_counts_list <- list()


for (dir in c("A", "B", "C", "D", "E", "G")) {
  for (fname in c("H_F", "H_T", "K_F", "K_T", "L_F", "L_T")) {
    filepath <- file.path(dir, paste0(fname, ".tsv"))
    repertoire <- read.table(filepath, header = TRUE, sep = "\t")
    kmers <- table(nchar(unlist(strsplit(repertoire$cdr3_aa, "\\*"))))
    counts <- kmers[c("3", "4", "5", "6","7", "8", "9")]
    counts_df <- as.data.frame(as.list(counts))
    print(counts_df)
    
    rowname <- paste(dir, fname, sep = "_")
    counts_df$File <- rowname
    all_counts_list[[rowname]] <- counts_df
  }
}

repertoire <- read.table("G/L_T.tsv", header = TRUE, sep = "\t")
kmers <- table(nchar(unlist(strsplit(repertoire$cdr3_aa, "\\*"))))
counts <- kmers[c("3", "4", "5", "6","7", "8", "9")]
counts_df <- as.data.frame(as.list(counts))
print(counts_df)

counts_df$File <- "G_L_T"
counts_df

all_counts_list[["G_L_T"]] <- counts_df

all_counts_list
combined_df <- rbindlist(all_counts_list, fill = TRUE)
combined_df
combined_
filenames <- as.vector(t(outer(c("A", "B", "C", "D", "E", "G"), c("H_F", "H_T", "K_F", "K_T", "L_F", "L_T"), paste, sep = "_")))
rownames(combined_df) <- filenames
combined_df <- combined_df[-8]
combined_df
combined_df <- a
yoyo <- as.data.frame(combined_df)
class(yoyo)
yoyo
filenames <- as.vector(t(outer(c("A", "B", "C", "D", "E", "G"), c("H_F", "H_T", "K_F", "K_T", "L_F", "L_T"), paste, sep = "_")))
filenames
rownames(yoyo) <- filesnames
yoyo
yoyo <- yoyo[-8]
yoyo
class(yoyo)
fwrite(yoyo, "kmer_frequencies.csv", row.names=TRUE)

#gpt did the write-csv part

library(data.table)

dirs <- LETTERS[1:7]
files <- c("H_F", "H_T", "K_F", "K_T", "L_F", "L_T")

# Prepare an empty list to collect results
results <- list()

for (dir in dirs) {
  for (file in files) {
    path <- file.path(dir, paste0(file, ".tsv"))
    if (!file.exists(path)) {
      warning(paste("File not found:", path))
      next
    }
    dat <- fread(path)
    mean_length <- mean(nchar(dat$cdr3_aa))
    # print(file)
    # 
    # mean_length <- mean(dat$cdr3_end[valid_rows] - dat$cdr3_start[valid_rows])
    # print(mean_length)
    
    results[[length(results) + 1]] <- data.table(
      Directory = dir,
      File = file,
      MeanCDR3Length = mean_length
    )
  }
}

# Combine all results into one data table
all_results <- rbindlist(results)

# Write to CSV
fwrite(all_results, "cdr3_means.csv")

print("Saved means to cdr3_means.csv")

kmer_freq <- read.csv("kmer_frequencies.csv")
rownames(kmer_freq) <- kmer_freq$X
kmer_freq <- kmer_freq[-1]
#all
fetal <- kmer_freq[grepl("_F", rownames(kmer_freq)),]
transitional <- kmer_freq[grepl("_T", rownames(kmer_freq)),]
ks <- intersect(names(fetal), names(transitional))
p_values <- sapply(ks, function(k){
  wilcox.test(unlist(fetal[k], use.names=FALSE), unlist(transitional[k], use.names=FALSE), paired=TRUE)$p.value
})
sorted_p_values <- p_values[order(p_values, decreasing=TRUE)]
write.table(sorted_p_values, sep = "\t", quote = FALSE)

fetal
transitional


library(data.table)

all_counts_list <- list()

for (dir in c("A", "B", "C", "D", "E", "G")) {
  for (fname in c("H_F", "H_T", "K_F", "K_T", "L_F", "L_T")) {
    filepath <- file.path(dir, paste0(fname, ".tsv"))
    repertoire <- read.table(filepath, header = TRUE, sep = "\t")
    aa_soup <- paste0(repertoire$cdr3_aa, collapse="")
    counts <- table(strsplit(aa_soup, ""))
    
    counts_df <- as.data.frame(as.list(counts))  # one row, named columns
    
    rowname <- paste(dir, fname, sep = "_")
    
    counts_df$File <- rowname
    
    all_counts_list[[rowname]] <- counts_df
    
    print(all_counts_list)
    # Combine all into one data.frame and fill NAs with 0
    all_counts_df <- rbindlist(all_counts_list, fill = TRUE)
    all_counts_df[is.na(all_counts_df)] <- 0
    
    # Add rownames as a column for clarity
    #all_counts_df <- cbind(File = rownames(all_counts_df), all_counts_df)
    rownames(all_counts_df) <- NULL
    
    # See it!
  }
}
# Make sure something was collected
if (length(all_counts_list) == 0) stop("No data was collected.")

combined_df <- rbindlist(all_counts_list, fill = TRUE)
combined_df[is.na(combined_df)] <- 0

# Now add filenames

# combined_df <- read.csv("2cdr3_character_frequencies.csv") this was to fix the old file
filenames <- as.vector(t(outer(c("A", "B", "C", "D", "E", "G"), c("H_F", "H_T", "K_F", "K_T", "L_F", "L_T"), paste, sep = "_")))

filenames <- filenames[1:nrow(combined_df)]
combined_df$File <- filenames
setcolorder(combined_df, c("File", setdiff(names(combined_df), "File")))

print(combined_df)
fwrite(combined_df, "cdr3_character_frequencies.csv")
