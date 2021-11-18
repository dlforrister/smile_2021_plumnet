library(muma)
library(pvclust)
library(vegan)
library(splitstackshape)
library(ape)


as.phylo.pvclust <- function (pvclust_object,nodelabels="AU",roundvalue=2,terminals=0.05,...) 
{
  x <- pvclust_object$hclust
  #number of groups in hclust object
  N <- dim(x$merge)[1]
  #getting edge matrix
  edge <- matrix(0L, 2 * N, 2)
  #getting edge.length vector
  edge.length <- numeric(2 * N)
  #getting vector of nodes
  node <- integer(N)
  node[N] <- N + 2L
  cur.nod <- N + 3L
  j <- 1L
  for (i in N:1) {
    edge[j:(j + 1), 1] <- node[i]
    for (l in 1:2) {
      k <- j + l - 1L
      y <- x$merge[i, l]
      
      if (y > 0) {
        edge[k, 2] <- node[y] <- cur.nod
        cur.nod <- cur.nod + 1L
        edge.length[k] <- x$height[i] - x$height[y]
      }
      else {
        edge[k, 2] <- -y
        if (terminals=="NA"){
          edge.length[k] <- x$height[i]
        }
        else {
          edge.length[k] <- terminals
        }
      }
    }
    j <- j + 2L
  }
  node.label <- numeric(N)
  if (nodelabels=="AU"){
    probs <- round(pvclust_object$edges[,1],digits=roundvalue)
    node.label[node-N-1] <- probs
    node.label[1] <- NA
  }
  else {
    if (nodelabels=="BP"){
      probs <- round(pvclust_object$edges[,2],digits=roundvalue)
      node.label[node-N-1] <- probs
      node.label[1] <- NA
    }
    else {
      node.label[] <- NA   
    }
  }
  if (is.null(x$labels)) 
    x$labels <- as.character(1:(N + 1))
  obj <- list(edge = edge, edge.length = edge.length/2, tip.label = x$labels, 
              Nnode = N,node.label=node.label)
  class(obj) <- "phylo"
  reorder(obj)
  
}




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

plot(species_centroid)

result_species <- pvclust(sample_dist, method.hclust=mhc, method.dist="correlation", use.cor="pairwise.complete.obs", nboot=1000,parallel=T,)

plot(result_species)
plot(result_species)


chemocode_tree<-as.phylo.pvclust(result_species)
write.tree(chemocode_tree,file="../../data/qemistree_test_1_evaluation/true_sample_distance.tre"))

mat <- dist(df, diag = TRUE, upper = FALSE)
species_dist <- as.dist(1-cor(sample_dist, method="pearson"),diag = TRUE, upper = FALSE)	
mat2 <- as.matrix(species_dist)
mat2[upper.tri(mat2, diag = FALSE)] <- ""


write.csv(mat2,"../../data/qemistree_test_1_evaluation/true_sample_distance_as_dist.csv")

