library(tidyverse)
library(stringr)
library(readr)

#read gene file
file_path <- "\home\abirami\Downloads\task\hiring_task\Homo_sapiens.gene_info.gz"

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
file_path <- "\home\abirami\Downloads\task\hiring_task\h.all.v2023.1.Hs.symbols.gmt"
gene_sets <- read_gmt(file_path)

# Check if there are any gene sets
if (length(gene_sets) == 0) {
  
  print("No gene sets found.")
} else {
  # Print all gene sets
  print_all_gene_sets(gene_sets)
}
