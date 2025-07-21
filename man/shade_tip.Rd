\name{shade_tip}
\alias{shade_tip}
\title{
Tip-based Core Community Phylogeny Using the Shade & Stopnisek (2019) algorithm.
}
\description{
Called internally. Identifies all edges and core edges of the core community phylogeny using a tip-based approach with a modified Shade and Stopnisek (2019) algorithm.
}
\usage{
shade_tip(xv, nt, mxtx,fi)
}
\arguments{
\item{xv}{
  (Required) The phyloseq object passed from main functions and containing microbial community data.
}
\item{nt}{
(Required) The microbial phylogeny passed from main functions.
}
\item{mxtx}{
(Required) The maximum number of branches to add sequentially, as a percentage of the total branches.
}
\item{fi}{
(Required) The threshold for the percent increase in beta diversity used to identify the taxon at which point adding more taxa yields diminishing returns in explanatory power.
}

}
\details{
\code{shade_tip} is used internally in the holobiont package to identify core edges of a microbial phylogeny for the tip-based approach using a modified Shade and Stopnisek (2019) algorithm. Briefly, taxa are ranked according to abundance first and then occupancy (i.e., occupancy is more important). Next, taxa are added sequentially, up to a maximum number of taxa. Each time a new taxon is added, the mean beta-diversity (UniFrac) between each pair of samples is calculated and is then averaged across all samples. Next, total beta-diversity (all taxa) between pairs samples is calculated, and again averaged across all samples. Finally, the percent increase in total beta-diversity is calculated for each new taxon added. A threshold is then selected based on the lowest ranked taxon that yields a minimum percent increase in beta-diversity.
}
\value{
This function returns a list of all edges and core edges, as well s all taxa and all core taxa.
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

shade_tip(example_phyloseq, newtree, NULL, 2)

}
