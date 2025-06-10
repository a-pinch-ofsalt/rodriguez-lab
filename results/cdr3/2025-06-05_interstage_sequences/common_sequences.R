library(data.table)

# Directories and chain types
prefixes <- c("A", "B", "C", "D", "E", "G")
chains <- c("L_T", "H_T", "K_T")

# Function to calculate pairwise stage similarity matrix
stage_similarity_matrix <- function(dt) {
  dt <- dt[!is.na(cdr3_aa) & cdr3_aa != ""]  # Clean data
  stages <- sort(unique(dt$stage))
  mat <- matrix(0, nrow = length(stages), ncol = length(stages),
                dimnames = list(paste0("Stage_", stages), paste0("Stage_", stages)))
  
  for (i in stages) {
    seqs_i <- unique(dt[stage == i, cdr3_aa])
    for (j in stages) {
      seqs_j <- unique(dt[stage == j, cdr3_aa])
      shared <- length(intersect(seqs_i, seqs_j))
      total <- length(unique(c(seqs_i, seqs_j)))
      mat[paste0("Stage_", i), paste0("Stage_", j)] <- shared / total * 100  # percent shared
    }
  }
  return(mat)
}

# Store matrices
similarity_matrices <- list()

# Main loop
for (prefix in prefixes) {
  for (chain in chains) {
    file_path <- file.path(prefix, paste0(chain, ".tsv"))
    if (file.exists(file_path)) {
      dt <- fread(file_path)
      mat <- stage_similarity_matrix(dt)
      similarity_matrices[[paste(prefix, chain, sep = "_")]] <- mat
    }
  }
}

# Example: view one of the matrices
print(similarity_matrices[["A_L_T"]])

