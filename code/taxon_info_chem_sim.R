pfSpecies
nrow(pfSpecies)
trop_sp <- pfSpecies[sample(1:nrow(pfSpecies),243,replace = F),c("species_id","family","genus","accepted_name")]
temp_sp <- edSpecies[sample(1:nrow(edSpecies),243,replace = F),c("species_id","family","genus","accepted_name")]
species_list <- rbind(trop_sp,temp_sp)
nrow(species_list)

#install.packages('myTAI')
library(myTAI)

species_names <- species_list$accepted_name



ncbi_taxon_loop <- function(species_names){
taxon_id <- NULL
family_id <- NULL
genus_id <- NULL

for (i in c(species_names)) {
  tmp_id <- taxonomy(organism = i, db = "ncbi", output = "classification" )
  
  # NCBITaxon not found
  if (length(tmp_id) <= 1) {
    taxon_id <- c(taxon_id, 0)
    # NCBITaxon found
  } else {
    taxon_id <- c(taxon_id, paste0(tmp_id$id[which(tmp_id$rank == "species")], "|", tmp_id$name[which(tmp_id$rank == "species")]))
    family_id <- c(family_id, paste0(tmp_id$id[which(tmp_id$rank == "family")], "|", tmp_id$name[which(tmp_id$rank == "family")]))
    genus_id <-c(genus_id, paste0(tmp_id$id[which(tmp_id$rank == "genus")], "|", tmp_id$name[which(tmp_id$rank == "genus")]))
  }}

taxon_tabl <- data.frame("species_name" = species_names, "NCBIfamily" = family_id,"NCBIgenus" = genus_id,"NCBItaxon" = taxon_id )  
return(taxon_tabl)
}

#species level
taxon_tabl <- ncbi_taxon_loop(species_list$accepted_name)






#### SAMPLE DISTANCE
sample_dist <- read.csv("../../data/qemistree_test_1_evaluation/true_sample_distance.csv",stringsAsFactors = T)
row.names(sample_dist) <- sample_dist$filename
sample_dist <- t(sample_dist[,2:7])

species_dist <- as.dist(1-cor(sample_dist, method="pearson"))				# Distance matrix
species_single <- hclust(species_dist, method="single")					# Single linkage clustering
species_ward <- hclust(species_dist, method="ward.D")					# Ward clustering
species_complete <- hclust(species_dist, method="complete")				# Complete linkage clustering
species_centroid <- hclust(species_dist, method="centroid")				# Centroid clustering
species_median <- hclust(species_dist, method="median")				        # Median clustering
#Comparison between the distance matrix and binary matrices representing partitions 
coph1 <- cophenetic(species_single)							# Compute Patristic distances		
coph2 <- cophenetic(species_ward)
coph3 <- cophenetic(species_complete)
coph4 <- cophenetic(species_centroid)
coph5 <- cophenetic(species_median)
a <-cor(coph1, species_dist)								#Cophenetic correlations
b <-cor(coph2, species_dist)
c <-cor(coph3, species_dist)
d <-cor(coph4, species_dist)
e <-cor(coph5, species_dist)
method <- c("single","ward","complete", "centroid", "median")
mhc<- method[which.max(c(a,b,c,d,e))]

result_species <- pvclust(sample_dist, method.hclust=mhc, method.dist="correlation", use.cor="pairwise.complete.obs", nboot=1000,parallel=T,)
plot(result_species)
chemocode_tree<-as.phylo.pvclust(result_species)
write.tree(chemocode_tree,file=paste("../../data/qemistree_test_1_evaluation/true_sample_distance.csv.tre",sep=""))
