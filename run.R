library(tidyverse)
library(stringr)
library(readr)

#read gene file
file_path <- "Homo_sapiens.gene_info.gz"

data <- read_tsv(gzfile(file_path))
cols_to_drop <- c("#tax_id", "LocusTag", "dbXrefs", "chromosome", "map_location", "description", 
                  "type_of_gene", "Symbol_from_nomenclature_authority", 
                  "Full_name_from_nomenclature_authority", "Nomenclature_status", 
                  "Other_designations", "Modification_date", "Feature_type")

filtered_data <- data[, !(names(data) %in% cols_to_drop)]

exploded_data <- filtered_data %>% 
  separate_rows(Synonyms, sep = "\\|")

head(exploded_data)

#read_gmt is a function to read a GMT file

read_gmt <- function(file_path) {
  # Check if the file exists
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  # Read the file line by line
  lines <- readLines(file_path)
  
  # Check if the file is empty
  if (length(lines) == 0) {
    stop("File is empty: ", file_path)
  }
  
  # Initialize the list to store the gene sets
  gene_sets <- list()
  
  # Process each line
  for(line in lines) {
    # Split the line into components
    split_line <- strsplit(line, "\t")[[1]]
    
    # Extract the gene set name, description, and genes
    gene_set_name <- split_line[1]
    description <- split_line[2]
    genes <- split_line[3:length(split_line)]
    
    # If there is no description, adjust accordingly
    if(description %in% genes) {
      genes <- genes[genes != description]
      description <- NA
    }
    
    # Add the gene set to the list
    gene_sets[[gene_set_name]] <- list(description = description, genes = genes)
  }
  
  return(gene_sets)
}

# Define the function to print all gene sets
print_all_gene_sets <- function(gene_sets) {
  for(gene_set in names(gene_sets)) {
    cat("Gene Set:", gene_set, "\n")
    cat("Description:", gene_sets[[gene_set]]$description, "\n")
    cat("Genes:", paste(gene_sets[[gene_set]]$genes, collapse = ", "), "\n\n")
  }
}

# Use the function
file_path <- "h.all.v2023.1.Hs.symbols.gmt"
gene_sets <- read_gmt(file_path)

# Check if there are any gene sets
if (length(gene_sets) == 0) {
  
  print("No gene sets found.")
} else {
  # Print all gene sets
  print_all_gene_sets(gene_sets)
}

symbol_to_id <- exploded_data %>%
  select(Symbol, GeneID) %>%
  deframe()

# Split synonyms and create a mapping
exploded_data <- exploded_data %>%
  mutate(Synonyms = str_split(Synonyms, pattern = "\\|")) %>%
  unnest(Synonyms)

synonyms_to_id <- exploded_data %>%
  select(Synonyms, GeneID) %>%
  deframe()

# Function to read the gmt file and replace the gene names with IDs
read_gmt <- function(file_path, symbol_to_id, synonyms_to_id) {
  lines <- readLines(file_path)
  
  gene_sets <- map(lines, function(line) {
    parts <- str_split(line, "\t")[[1]]
    
    gene_set_name <- parts[1]
    description <- parts[2]
    genes <- parts[3:length(parts)]
    
    # If description is present in genes, remove it
    if (description %in% genes) {
      genes <- setdiff(genes, description)
      description <- NULL
    }
    
    # Replace the gene names with their IDs
    genes <- map_chr(genes, function(gene) {
      if (gene %in% names(symbol_to_id)) {
        as.character(symbol_to_id[[gene]])
      } else if (gene %in% names(synonyms_to_id)) {
        as.character(synonyms_to_id[[gene]])
      } else {
        NA_character_
      }
    })
    
    list(gene_set_name = gene_set_name, description = description, genes = genes)
  })
  
  gene_sets
}

# Function to print all gene sets
print_all_gene_sets <- function(gene_sets) {
  for (gene_set in gene_sets) {
    cat(paste("Gene Set: ", gene_set$gene_set_name, "\n"))
    cat(paste("Description: ", gene_set$description, "\n"))
    cat(paste("Genes: ", paste(gene_set$genes, collapse = ", "), "\n\n"))
  }
}

# Function to write the gene sets to a .gmt file
write_gmt <- function(gene_sets, file_name) {
  file_conn <- file(file_name, "w")
  for (gene_set in gene_sets) {
    lines <- c(paste("Gene Set: ", gene_set$gene_set_name), 
               paste("Description: ", gene_set$description), 
               paste("Genes: ", paste(gene_set$genes, collapse = ", ")))
    writeLines(lines, file_conn)
    writeLines("", file_conn) # add a newline after each gene set
  }
  close(file_conn)
}

# Use the function
gene_sets <- read_gmt("h.all.v2023.1.Hs.symbols.gmt", symbol_to_id, synonyms_to_id)

# Print all gene sets
print_all_gene_sets(gene_sets)

