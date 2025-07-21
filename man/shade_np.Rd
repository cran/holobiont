\name{shade_np}
\alias{shade_np}
\title{
Non-phylogenetic Core Community Using the Shade & Stopnisek (2019) algorithm.
}
\description{
Called internally. Identifies all taxa and core taxa of the core community using a non-phylogenetic approach with a modified Shade and Stopnisek (2019) algorithm.
}
\usage{
shade_np(xv, mxtx,fi)
}
\arguments{
\item{xv}{
  (Required) The phyloseq object passed from main functions and containing microbial community data.
}

\item{mxtx}{
(Required) The maximum number of branches to add sequentially, as a percentage of the total branches.
}
\item{fi}{
(Required) The threshold for the percent increase in beta diversity used to identify the taxon at which point adding more taxa yields diminishing returns in explanatory power.
}

}
\details{
\code{shade_np} is used internally in the holobiont package to identify core edges of a microbial phylogeny for the non-phylogenetic approach using a modified Shade and Stopnisek (2019) algorithm. Briefly, taxa are ranked according to abundance first and then occupancy (i.e., occupancy is more important). Next, taxa are added sequentially, up to a maximum number of taxa. Each time a new taxon is added, the mean beta-diversity (Jaccard) between each pair of samples is calculated and is then averaged across all samples. Next, total beta-diversity (all taxa) between pairs samples is calculated, and again averaged across all samples. Finally, the percent increase in total beta-diversity is calculated for each new taxon added. A threshold is then selected based on the lowest ranked taxon that yields a minimum percent increase in beta-diversity.
}
\value{
This function returns a list of all taxa and core taxa.
}
\references{

Shade, Ashley, and Nejc Stopnisek. "Abundance-occupancy distributions to prioritize plant core microbiome membership." Current opinion in microbiology 49 (2019): 50-58.
}


\examples{
#Test with enterotype dataset
library(phyloseq)
library(ape)
library(phytools)
data(enterotype)

set.seed(1)

enterotype<-prune_taxa(taxa_names(enterotype)[2:21],enterotype)
enterotype<-prune_samples(sample_names(enterotype)[2:21],enterotype)

#Generate an example tree and label it with the names of the microbial taxa
enterotype_tree<-rtree(length(taxa_names(enterotype)))
enterotype_tree$tip.label<-taxa_names(enterotype)
enterotype_tree$node.label<-as.character(1:1:enterotype_tree$Nnode)

#Create a phyloseq object with a tree
example_phyloseq<-phyloseq(otu_table(enterotype),phy_tree(as.phylo(enterotype_tree)))
newtree<-bind.tip(phy_tree(example_phyloseq),tip.label='outgroup',edge.length=0.0001,position=0)

shade_np(example_phyloseq, NULL, 2)

}
