## Load raw expression data for each developmental process, making sure gene ids are consistent
library(readr)
embryogenesis.phylo_expression_set <- read_delim("raw_data/Athaliana.csv", delim = "\t", 
                                               escape_double = FALSE, 
                                               col_types = cols(Phylostratum = col_integer()), 
                                               trim_ws = TRUE)

library(readxl)
germination.phylo_expression_set <- read_excel("raw_data/post-embryo data.xlsx", sheet = "Germination")
germination.phylo_expression_set$GeneID <- toupper(germination.phylo_expression_set$GeneID)

floral_transition.phylo_expression_set <- read_excel("raw_data/post-embryo data.xlsx", sheet = "Floral transition")
floral_transition.phylo_expression_set$GeneID <- toupper(floral_transition.phylo_expression_set$GeneID)

flower_development.phylo_expression_set <- read_excel("raw_data/post-embryo data.xlsx", sheet = "Flower Development")
flower_development.phylo_expression_set$GeneID <- toupper(flower_development.phylo_expression_set$GeneID)

# Apply sqrt transformation (temporary until GATAI can do it)
embryogenesis.phylo_expression_set[3:9] <- lapply(embryogenesis.phylo_expression_set[3:9], sqrt)
germination.phylo_expression_set[3:9] <- lapply(germination.phylo_expression_set[3:9], sqrt)
floral_transition.phylo_expression_set[3:9] <- lapply(floral_transition.phylo_expression_set[3:9], sqrt)
flower_development.phylo_expression_set[3:9] <- lapply(flower_development.phylo_expression_set[3:9], sqrt)

## Use updated phylomap of athaliana
library(phylomapr)
athaliana.phylo_map.uniprot <- phylomapr::Arabidopsis_thaliana.PhyloMap

# Need to convert gene ids from UNIPROT format into TAIR format
athaliana.phylo_map.tair.new <- phylomapr::convertID(phylomap = athaliana.phylo_map.uniprot, mart = "plants_mart", dataset = "athaliana_eg_gene")

# Update the gene ages
library(myTAI)
embryogenesis.phylo_expression_set.new <- myTAI::MatchMap(athaliana.phylo_map.tair.new, embryogenesis.phylo_expression_set[ , 2:9])
germination.phylo_expression_set.new <- myTAI::MatchMap(athaliana.phylo_map.tair.new, germination.phylo_expression_set[ , 2:9])
floral_transition.phylo_expression_set.new <- myTAI::MatchMap(athaliana.phylo_map.tair.new, floral_transition.phylo_expression_set[ , 2:9])
flower_development.phylo_expression_set.new <- myTAI::MatchMap(athaliana.phylo_map.tair.new, flower_development.phylo_expression_set[ , 2:9])

# Compare old phylomap with new phylomap
library(dplyr)
athaliana.phylo_map.tair.old <- unique(bind_rows(germination.phylo_expression_set[1:2], 
                                                 floral_transition.phylo_expression_set[1:2], 
                                                 flower_development.phylo_expression_set[1:2]))

library(ggplot2)
# Show distribution of gene ages
ggplot(athaliana.phylo_map.tair.old, aes(x = Phylostratum)) + geom_histogram()
ggplot(athaliana.phylo_map.tair.new, aes(x = Phylostratum)) + geom_histogram()

# Show how gene ages have changed between the old and new phylomap in a scatterplot
colnames(athaliana.phylo_map.tair.old)[1] <- "Phylostratum old"
athaliana.phylo_map.tair.comparison <- merge(athaliana.phylo_map.tair.old, athaliana.phylo_map.tair.new)
ggplot(athaliana.phylo_map.tair.comparison, aes(x = `Phylostratum old`,y = Phylostratum)) + geom_bin2d() + geom_abline()

# Save processed data to file
write.table(embryogenesis.phylo_expression_set.new, "processed_data/embryogenesis.csv", sep='\t', row.names=FALSE)
write.table(germination.phylo_expression_set.new, "processed_data/germination.csv", sep='\t', row.names=FALSE)
write.table(floral_transition.phylo_expression_set.new, "processed_data/floral_transition.csv", sep='\t', row.names=FALSE)
write.table(flower_development.phylo_expression_set.new, "processed_data/flower_development.csv", sep='\t', row.names=FALSE)
