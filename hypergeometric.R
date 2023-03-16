# Install the following packages if you have not 
install.packages("phylolm")
install.packages("phytools")
install.packages("maps")
install.packages("ape")

# running hypergeometric distribution test  
# choose one of the four major phyla from plant-associated environment  
# Bacteroidetes, Firmicutes, Actinobacteria, Proteobacteria
rm(list = ls())
#setwd("~/Globus/bioinfo/BGC")

# for family-based BGC matrix
GENOME<- "genome_database.csv"
phylum = "Bacteroidetes"
BGC <- paste("BGCF_matrix_0.2_", phylum, ".csv", sep = "")

# PA versus NPA
suppressMessages(library(phytools))
print ("PA versus NPA")
print(BGC)
print(GENOME)

# Read in data
Metadata <- read.csv(GENOME, header=T)
bgc_count <- read.csv(BGC, header=T, row.names = 1)
#bgc_count<-t(bgc_count)
# first filter genomes without any BGC predicted
Metadata <- Metadata[match(rownames(bgc_count),Metadata$filename),]
# Choose one of the phyla (Here we choose Firmicutes)
Metadata <- Metadata[which(Metadata$taxonomy_NCBI_2 == phylum),]
Metadata <- droplevels(Metadata)
bgc_count <- bgc_count[match(Metadata$filename,rownames(bgc_count)),]
# Remove BGC classes that have no prediction in this phylum
bgc_count <- bgc_count[,which(colSums(bgc_count)!=0)]

# Read the phylogenetic tree
tree<-read.tree("Tree_of_All")

# Subset the phylogenetic tree
subtree <- drop.tip(phy=tree,tip=which(!(tree$tip.label%in%rownames(bgc_count))))

# functions for hypergeomteric distribution calculation
group_hyp <- function(tree,groups,metadata){
  # the total BGC number from PA/NPA bacteria of this phylum
  N.pa <- sum(groups[ metadata$Classification == "PA", ])
  N.npa <- sum(groups[ metadata$Classification == "NPA", ])
  # N.genomes <- nrow(Map)
  Res <- NULL
  for(group in 1:ncol(groups)){
    # gene <- 1
    print (group)
    group_name <- colnames(groups)[group]
    metadata$group <- groups[,group]
    
    # Binary version
    # q = # of PA genomes with this class of BGCs
    # m = total # of the BGCs of this class in the dataset
    # n = # of genomes without this class of BGCs
    # k = # of PA genomes in the dataset
    # lower.tail: calculating the probablity of getting more than q - 1 white balls 
    # the model is: test if it follows hypergeometric distribution, if not, then this gene 
    # is either enriched or depleted from 
    pval <- phyper(q = nrow(subset(metadata, Classification == "PA" & group > 0)) - 1,
                   m = sum(metadata$group > 0),
                   n = nrow(subset(metadata,group == 0)),
                   k = nrow(subset(metadata, Classification == "PA")),lower.tail=FALSE)
    score <- -log10(pval)
    # rawcounts
    pval2 <- phyper(q = sum(subset(metadata, Classification == "PA")$group) - 1, 
                    m = N.pa, n = N.npa, k = sum(metadata$group),lower.tail=FALSE)
    
    
    # Test for depletion binary version 
    pval_dep <- phyper(q = nrow(subset(metadata, Classification == "PA" & group > 0)),
                       m = sum(metadata$group > 0),
                       n = nrow(subset(metadata,group == 0)),
                       k = nrow(subset(metadata, Classification == "PA")),
                       lower.tail=TRUE)
    
    score_dep<--log10(pval_dep)
    
    # rawcounts
    pval2_dep <- phyper(q = sum(subset(metadata, Classification == "PA")$group),
                        m = N.pa, n = N.npa, k = sum(metadata$group),lower.tail=TRUE)
    score2_dep<--log10(pval2_dep)
    
    
    
    res <- data.frame(BGC_class = group_name, score_enriched_binary = score, p.value_enriched_binary = pval,
                      score_depletion_binary=score_dep,p.value_depletion_binary=pval_dep,
                      score_enriched_rawcounts = -log10(pval2), p.value_enriched_rawcounts = pval2,
                      score_depletion_rawcounts=score2_dep,p.value_depletion_rawcounts=pval2_dep)
    
    metadata$group <- NULL
    Res <- rbind(Res,res)
  }
  
  Res<-data.frame(BGC_class=Res$BGC_class,score_enriched_binary=Res$score_enriched_binary,
                  z.score_enriched_binary= (Res$score_enriched_binary - mean(Res$score_enriched_binary)) / sd(Res$score_enriched_binary),
                  p.value_enriched_binary=Res$p.value_enriched_binary,
                  score_depletion_binary=Res$score_depletion_binary,
                  z.score_depletion_binary=(Res$score_depletion_binary - mean(Res$score_depletion_binary)) / sd(Res$score_depletion_binary),
                  p.value_depletion_binary=Res$p.value_depletion_binary,
                  score_enriched_rawcounts=Res$score_enriched_rawcounts,
                  z.score_enriched_rawcounts=(Res$score_enriched_rawcounts - mean(Res$score_enriched_rawcounts)) / sd(Res$score_enriched_rawcounts),
                  p.value_enriched_rawcounts=Res$p.value_enriched_rawcounts,
                  score_depletion_rawcounts=Res$score_depletion_rawcounts,
                  z.score_depletion_rawcounts=(Res$score_depletion_rawcounts - mean(Res$score_depletion_rawcounts)) / sd(Res$score_depletion_rawcounts),
                  p.value_depletion_rawcounts=Res$p.value_depletion_rawcounts)
  return(Res)
}

finRes<-group_hyp(subtree,bgc_count,Metadata)
# BGC_class_based
# filename = paste("hypergeometric_", phylum ,"_pa_npa.csv", sep = "")
filename = paste("hypergeometric_", phylum ,"_BGCF.csv", sep = "")
print(filename)
write.csv(finRes,file = filename,
            row.names = FALSE, quote = FALSE)
