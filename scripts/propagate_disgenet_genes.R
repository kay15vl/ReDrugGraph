#Created by: Karla Vela Lopez
#Created on: 4/1/24

#Input: Disease gene association data from DisGeNET
# Output: Propagated genes associated to of disease of interest (TB).

library(tidyverse)
library(ontologyIndex)
library(parallel)
library(data.table)

mondo_ont <- get_OBO(file = "data/raw/mondo_2024_03.obo",
                     propagate_relationships = "is_a",
                     extract_tags = "everything")
df_mondo_ont <- as.data.frame(mondo_ont)


#get disgenet data
disgenet_data <- read_tsv("data/bigdatasets/all_gene_disease_associations.tsv") %>% 
  filter(score >=.1) %>% 
  select(-DSI,-DPI,-diseaseType:-source) %>% 
  rename(disease_id=diseaseId, disease_name=diseaseName) %>% 
  mutate(disease_id = str_replace(disease_id, "^C", "UMLS:C"))  %>% setDT()


map_mondo_to_ontology <- function(dis_gene_data, obo_file) {
  
  mondo_to_disid <- as.data.frame(obo_file) %>%
    dplyr::select(
      -`format-version`:-ontology, -synonym, -comment:-union_of,
      -seeAlso:-is_class_level, -parents:-ancestors, -property_value:-def,
      -replaced_by:-is_a
    ) %>%
    mutate(xref = lapply(xref, function(x) gsub("; ", ";", x))) %>%
    separate_rows(xref, sep = ";") %>%
    filter(!is.na(obsolete) & obsolete != TRUE) %>%
    dplyr::rename(mondo_id = id, mondo_term = name, ontology_dis_id = xref) %>%
    filter(!is.na(ontology_dis_id) & ontology_dis_id != "") %>% 
    rename(disease_id = ontology_dis_id) %>% 
    mutate(disease_id = str_replace(disease_id, "Orphanet", "ORPHA"))
  
  # join mondo df to the gwas dataframe. Join used to map mondo to the dise
  convert_to_mondo <- dis_gene_data %>%
    left_join(mondo_to_disid, by = "disease_id") %>%
    relocate(c("mondo_id", "mondo_term"), .after = disease_id) %>%
    dplyr::select(-obsolete, -subset) %>%
    filter(!is.na(mondo_id)) %>%
    unique() 
  
  return(convert_to_mondo)
}

#funciton to get ancestor information
get_ances<- function(term_id, ontology_obo){
  
  # get parent terms
  ancestors<- get_ancestors(ontology_obo, terms= term_id )
  
  # init df
  parent_df <- tibble()
  
  # add data to df
  parent_df <- tibble(ID= term_id,
                      Parent= ancestors)
  
  # collapse ancestors to one row per go term id
  parent_df <- parent_df %>% filter(Parent != term_id) %>% 
    group_by(ID) %>% 
    summarise(parent_ct = n_distinct(Parent),
              Parent = paste(Parent, collapse = ","))
  return(parent_df)
  
}

# Function to propagate genes from child terms to parent terms
propagate_genes<- function(current_term, ance_df, disgenes_df){
  #get ancestor terms
  ancestor_terms <- ance_df %>% filter(ID == current_term) %>%select(Parent) %>%
    unlist() %>% strsplit( ",") %>% unlist() %>%trimws() %>%  as.vector()
  
  #find genes of term in dataframe
  current_term_disgenes <- disgenes_df %>% filter(mondo_id== current_term) 
  
  ancestor_disgenes <- disgenes_df %>% filter(mondo_id %in% ancestor_terms)
  
  # add to parent terms
  unique_mondo_ids <- unique(ancestor_disgenes$mondo_id)
  unique_mondo_term <- unique(ancestor_disgenes$mondo_term)
  
  # Repeat gene_id values for each unique mondo_id
  repeated_gene_ids <- rep(current_term_disgenes$gene_id, length(unique_mondo_ids))
  
  # Create a dataframe with repeated gene_id values for each unique mondo_id
  new_rows <- data.frame(
    gene_id = repeated_gene_ids,
    mondo_id = rep(unique_mondo_ids, each = nrow(current_term_disgenes)),
    mondo_term = rep(unique_mondo_term, each = nrow(current_term_disgenes))
  )
  
  propagated_df <-ancestor_disgenes %>%  bind_rows(new_rows, current_term_disgenes) %>% 
    unique()
  
  return(propagated_df)
}

####### DisGeNET Gene Collection #######
final_disgenet <- map_mondo_to_ontology(disgenet_data, df_mondo_ont)
# get rid of uneccesary col
final_disgenet <- final_disgenet %>% select(-geneSymbol, -disease_id,
                                            -disease_name) %>% 
  rename("gene_id"= "geneId") %>% unique()

write_tsv(final_disgenet, 
          file = "data/processed/not_propagated_disgenet_genes.tsv")
terms <- unique(final_disgenet$mondo_id)
ancestor_info <- mclapply(terms, get_ances,mondo_ont) %>% bind_rows() %>% 
  unique() %>% setDT()

propagated_genes <- mclapply(terms, propagate_genes, ancestor_info, final_disgenet) %>% 
  bind_rows() %>% unique() %>% setDT()

final_prop_df <- bind_rows(final_disgenet, propagated_genes) %>% 
  unique()

final_tb_data <- final_prop_df %>% filter(mondo_id == "MONDO:0018076")

write_tsv(final_tb_data, 
          file = "data/processed/tb_propagated_disgenet_genes.tsv")
save(final_prop_df, 
     file = "data/processed/propagated_disgenet_genes.RData")

print("Script successfully executed")
# diseases associated to tb phenos OMIM:300645, OMIM:619644, ORPHA:183675


