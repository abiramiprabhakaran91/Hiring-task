# Load necessary libraries
library(dplyr)
library(ggplot2)

# Read the tab-separated file
data <- read.delim("Homo_sapiens.gene_info", sep="\t")
column_name <- "chromosome"

item_counts <- data %>%
  group_by_at(column_name) %>%
  summarise(Count = n()) %>%
  ungroup()

# Rename the columns
item_counts <- item_counts %>%
  rename(chromosome = {{ column_name }}, gene_count = Count)

chromosome_order <- c(paste(1:22), "X", "Y", "MT", "Un")
item_counts$chromosome <- factor(item_counts$chromosome, levels = chromosome_order)

#print(item_counts$chromosome)

#drop na
cleaned_data <- na.omit(item_counts)
#print(cleaned_data)

plot <- ggplot(cleaned_data, aes(x = chromosome, y = gene_count, fill = Chromosome)) +
  geom_bar(stat = "identity", fill = "#737373")+
 theme_minimal() +
   labs(title = "Number of genes in each chromosome",
        x = "Chromosomes",
        y = "Gene count") +
   theme(plot.title = element_text(hjust = 0.5), axis.line.x.bottom=element_line(color="black"), axis.line.y.left=element_line(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
 
# Print the plot
#print(plot)
ggsave("task3.pdf", plot, width = 10, height = 6)
