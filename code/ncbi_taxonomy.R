load("/Users/dlforrister/Downloads/the_data.Rdata")

#taxonomy functions come from this rmarkdown
#https://github.com/drostlab/myTAI/blob/master/vignettes/Taxonomy.Rmd


pfSpecies
nrow(pfSpecies)
trop_sp <- pfSpecies[sample(1:nrow(pfSpecies),243,replace = F),c("species_id","family","genus","accepted_name")]
temp_sp <- edSpecies[sample(1:nrow(edSpecies),243,replace = F),c("species_id","family","genus","accepted_name")]
species_list <- rbind(trop_sp,temp_sp)
nrow(species_list)

#install.packages('myTAI')
library(myTAI)

install.packages(c("taxize", "usethis"))
taxize::use_entrez()
# Create your key from your (brand-new) account's. 
# After generating your key set it as ENTREZ_KEY in .Renviron.
# ENTREZ_KEY='youractualkeynotthisstring'
# For that, use usethis::edit_r_environ()
usethis::edit_r_environ()



species_names <- species_list$accepted_name



ncbi_taxon_loop <- function(species_names){
  taxon_id <- NULL
  family_id <- NULL
  genus_id <- NULL
  taxon_tabl <- data.frame()
  no_data <- NULL
  for (i in c(species_names)) {
    tmp_id <- taxonomy(organism = i, db = "ncbi", output = "classification" )
    
    # NCBITaxon not found
    if (length(tmp_id) <= 1) {
      no_data <- c(no_data,i)
      # NCBITaxon found
    } else {
      taxon_id <- paste0(tmp_id$id[which(tmp_id$rank == "species")], "|", tmp_id$name[which(tmp_id$rank == "species")])
      family_id <- paste0(tmp_id$id[which(tmp_id$rank == "family")], "|", tmp_id$name[which(tmp_id$rank == "family")])
      genus_id <- paste0(tmp_id$id[which(tmp_id$rank == "genus")], "|", tmp_id$name[which(tmp_id$rank == "genus")])
      taxon_tabl <- rbind(taxon_tabl,data.frame("species_name" = i, "NCBIfamily" = family_id,"NCBIgenus" = genus_id,"NCBItaxon" = taxon_id ))
    }}
  
  
  return(taxon_tabl)
}

#species level
taxon_tabl <- ncbi_taxon_loop(species_list$accepted_name)

write.csv(taxon_tabl,"~/Documents_Mac/CODE_GIT_HUB_2017_Aug_31/smile_2021_plumnet/data/qemistree_test_1_evaluation/random_taxonomy.csv")
