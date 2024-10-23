## Load raw expression data for each developmental process
library(readr)
embryogenesis.phylo_expression_set <- read_delim("raw_data/Athaliana.csv", delim = "\t", 
                                               escape_double = FALSE, 
                                               col_types = cols(Phylostratum = col_integer()), 
                                               trim_ws = TRUE)
embryogenesis.expression_matrix <- embryogenesis.phylo_expression_set[ , 2:9]

library(readxl)
germination.phylo_expression_set <- read_excel("raw_data/post-embryo data.xlsx", sheet = "Germination")
germination.expression_matrix <- germination.phylo_expression_set[ , 2:9]

floral_transition.phylo_expression_set <- read_excel("raw_data/post-embryo data.xlsx", sheet = "Floral transition")
floral_transition.expression_matrix <- floral_transition.phylo_expression_set[ , 2:9]

flower_development.phylo_expression_set <- read_excel("raw_data/post-embryo data.xlsx", sheet = "Flower Development")
flower_development.expression_matrix <- flower_development.phylo_expression_set[ , 2:9]



## Use updated phylomap of athaliana
library(phylomapr)
athaliana.phylo_map.uniprot <- phylomapr::Arabidopsis_thaliana.PhyloMap

# Need to convert gene ids
athaliana.phylo_map.

athaliana.phylo_map.tair <- germination.phylo_expression_set[ , 1:2]

library(myTAI)
embryogenesis.phylo_expression_set <- myTAI::MatchMap(athaliana.phylo_map, embryogenesis.expression_matrix)
germination.phylo_expression_set <- myTAI::MatchMap(athaliana.phylo_map, germination.expression_matrix)
floral_transition.phylo_expression_set <- myTAI::MatchMap(athaliana.phylo_map, floral_transition.expression_matrix)
flower_development.phylo_expression_set <- myTAI::MatchMap(athaliana.phylo_map, flower_development.expression_matrix)
