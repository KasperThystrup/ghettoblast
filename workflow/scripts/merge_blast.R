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


read_blast <- function(blast_file, col_names){
  blast_raw <- readr::read_tsv(
    file = blast_file,
    col_names = col_names,
    col_types = "diiiidciicciicc"
  )
  
  sample <- stringr::str_remove_all(string = blast_file, pattern = "\\.blast") %>%
    basename
  
  dplyr::mutate(blast_raw, sample = sample, .before = 1)
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
  tables <- purrr::map_dfr(.x = blast, .f = read_blast, col_names = col_names)
  
  dplyr::select(tables, dplyr::any_of(c("sample", col_names)))
}


recalculate_statistics <- function(blast_joined){
  
  message("INFO: Correcting identity and calculating coverage")
  dplyr::mutate(
    blast_joined,
    percent_query_coverage = alignment_length / query_length * 100,
    percent_identity = number_identical / query_length * 100,
    .after = percent_identity
  )
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
  fasta_string <- paste0(blast_sequences, query_accession, ".fasta")
  
  fasta_dir <- file.path(dirname(fasta_string), "Sequences")
  dir.create(path = fasta_dir, recursive = TRUE, showWarnings = FALSE)
  
  fasta_file <- stringr::str_replace_all(
    
    string = stringr::str_replace_all(
      string = basename(fasta_string), pattern = "\\W", replacement = "_"
    ),
    
    pattern = "_[fast]+$", replacement = ".fasta"
  )
  fasta_path <- file.path(fasta_dir, fasta_file)

  message("INFO: Writing fasta file - ", fasta_path)
  readr::write_lines(x = fasta, file = fasta_path, append = TRUE)
  
  message("Done")
}


merge_metadata <- function(blast_table, metadata){
  dplyr::left_join(blast_table, metadata, by = "sample")
}

extract_top_only <- function(blast){
  blast_by_sample <- dplyr::group_by(blast, sample)
  blast_filter_cov <- dplyr::filter(
    .data = blast_by_sample,
    percent_query_coverage == max(percent_query_coverage)
  )
  
  dplyr::filter(
    .data = blast_filter_cov,
    percent_identity == max(percent_identity)
  )
}

write_blast <- function(blast_merged, blast_results){
  message("INFO: Writing blast results table")
  readr::write_tsv(x = blast_merged, file = blast_results)
}


merge_blast <- function(blast, faidx, blast_columns, metadata_file, exclude_seqs, top_only, blast_results){
  
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
  blast_table <- recalculate_statistics(blast_joined)
  
  if(file.exists(metadata_file)){
    metadata <- read_metadata(metadata_file)
    blast_table <- merge_metadata(blast_table, metadata)
  }
  
  message("INFO: Cleaning output")
  blast_clean <- blast_table
  
  if (exclude_seqs){
    blast_fasta <- make_sequences(blast_table)
    blast_sequences <- stringr::str_remove(string = blast_results, pattern = "blast.tsv")
    blast_clean <- dplyr::select(blast_table, -c(subject_seq, query_seq))
    create_fasta(blast_fasta, blast_sequences)
  }
  if (top_only)
    blast_clean <- extract_top_only(blast = blast_clean)
  
  write_blast(blast_clean, blast_results)
  
  message("Done")
}
if (snakemake@params$debug | snakemake@params$exclude_seqs){
  tmp_dir <- file.path(snakemake@params[["outdir"]], "tmp")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  message("DEBUG: Saving image to ", tmp_file <- file.path(tmp_dir, "merge_blast.RData"))
  save.image(file = tmp_file)
}

message(paste("INFO: Merging", length(snakemake@input[["blast"]]), "samples"))

## Execute workflow
merge_blast(
  blast = snakemake@input$blast,
  faidx = snakemake@input$faidx,
  metadata_file = snakemake@params$metadata,
  exclude_seqs = snakemake@params$exclude_seqs,
  top_only = snakemake@params$top_only,
  blast_results = snakemake@output$blast_results
)

