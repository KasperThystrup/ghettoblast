library(magrittr)

read_blast_columns <- function(column_file){
  message("INFO: Reading columns from file - ", column_file)
  raw <- readr::read_tsv(file = column_file, col_types = "c")
  
  dplyr::pull(raw, column)
}
  

read_metadata <- function(metadata_file){
  message("INFO: Reading metadata file - ", metadata_file)
  metadata <- readr::read_tsv(file = metadata_file, show_col_types = FALSE)
  
  if(ncol(metadata) < 2){
    warning("There was only one column detected in the metadata file.\n",
    "Metadata will be ignored.")
    return(FALSE)
  }
  
  return(metadata)
}


read_blast <- function(blast_file, col_nms){
  blast_raw <- readr::read_tsv(
    file = blast_file,
    col_names = col_nms,
    col_types = "diiiidciicciicc"
  )
  
  sample <- stringr::str_remove_all(string = blast_file, pattern = "\\.blast") %>%
    basename
  
  dplyr::mutate(blast_raw, sample = sample)
}


read_faidx <- function(fai_file){
  message("INFO: Reading fasta index file ", fai_file)
  faidx_raw <- readr::read_tsv(
    file = fai_file,
    col_names = c("query_acc", "query_length", "offset", "linebases", "linewidth"),
    col_types = "ciiii"
  )
  
  dplyr::select(faidx_raw, query_acc, query_length)
}


merge_tables <- function(blast, col_names){
  message("INFO: Reading and merging all blast files")
  blast_list <- lapply(X = blast, FUN = read_blast, col_nms = col_names)
  
  do.call(blast_list, what = rbind)
}


calculate_coverage <- function(blast_joined){
  
  message("INFO: Calculating coverage and rearranging columns")
  blast_table <- dplyr::mutate(
    blast_joined,
    percent_query_coverage = alignment_length / query_length * 100
  ) 
  
  dplyr::relocate(blast_table, sample) %>%
    dplyr::relocate(query_length, percent_query_coverage, .after = gap_open)
}


make_sequences <- function(blast_joined){
  message("INFO: Making fasta sequences")
  dplyr::mutate(
    blast_joined,
    fasta = paste(
      paste0(">", paste(subject_acc, query_acc, sep = "_")),
      subject_seq,
      sep = "\n"
    )
  )
}


create_fasta <- function(blast_fasta, blast_sequences){
  message("INFO: Creating fasta sequences")
  query_accs <- dplyr::pull(blast_fasta, query_acc) %>%
    unique
  
  sapply(X = query_accs, FUN = write_sequences, blast_fasta = blast_fasta, blast_sequences = blast_sequences)
}


write_sequences <- function(query_accession, blast_fasta, blast_sequences){
  fasta <- subset(blast_fasta, query_acc == query_accession) %>%
    dplyr::pull(fasta)
  query_accession <- stringr::str_remove(string = query_accession, pattern = '[^a-zA-Z\\_0-9]')  ### Makes query string convertable to file_names
  fasta_file <- paste0(blast_sequences, query_accession, ".fasta")
  
  message("INFO: Writing fasta file - ", fasta_file)
  readr::write_lines(x = fasta, file = fasta_file, append = TRUE)
  
  message("Done")
}


merge_metadata <- function(blast_table, metadata){
  dplyr::left_join(blast_table, metadata, by = "sample")
}


write_blast <- function(blast_merged, blast_results){
  message("INFO: Writing blast results table")
  readr::write_tsv(x = blast_merged, file = blast_results)
}


merge_blast <- function(blast, faidx, blast_columns, metadata_file, exclude_seqs, blast_results){
  
  col_names <- c(
    "percent_identity",
    "number_identical",
    "alignment_length",
    "mismatches",
    "gap_open",
    "bit_score",
    "query_acc",
    "query_start",
    "query_end",
    "query_seq",
    "subject_acc",
    "subject_start",
    "subject_end",
    "subject_seq"
  )
  
  blast_merged  <- merge_tables(blast, col_names)
  query_index <- read_faidx(fai_file = faidx)
  
  message("INFO: Joing query index")
  blast_joined <- dplyr::inner_join(x = blast_merged, y = query_index, by = "query_acc")
  blast_table <- calculate_coverage(blast_joined)
  
  blast_fasta <- make_sequences(blast_table)
  blast_sequences <- stringr::str_remove(string = blast_results, pattern = "blast.tsv")
  create_fasta(blast_fasta, blast_sequences)
  
  if(file.exists(metadata_file)){
    metadata <- read_metadata(metadata_file)
    blast_table <- merge_metadata(blast_table, metadata)
  }
  
  message("INFO: Cleaning output")
  if (exclude_seqs)
    blast_table <- dplyr::select(blast_table, -c(subject_seq, query_seq))
  write_blast(blast_table, blast_results)
    
  message("Done")
}

tmp_dir <- file.path(snakemake@params[["outdir"]], "tmp")
dir.create(tmp_dir)
message("DEBUG: Saving image to ", tmp_file <- file.path(tmp_dir, "merge_blast.RData"))
save.image(file = tmp_file)

message(paste("INFO: Merging", length(snakemake@input[["blast"]]), "samples"))

## Execute workflow
merge_blast(
  blast = snakemake@input[["blast"]],
  faidx = snakemake@input[["faidx"]],
  metadata_file = snakemake@params[["metadata"]],
  exclude_seqs = snakemake@params[["exclude_seqs"]],
  blast_results = snakemake@output[["blast_results"]]
)

