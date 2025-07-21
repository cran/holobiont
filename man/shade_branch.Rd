\name{shade_branch}
\alias{shade_branch}
\title{
Branch-based Core Community Phylogeny Using the Shade & Stopnisek (2019) algorithm.
}
\description{
Called internally. Identifies all edges and core edges of the core community phylogeny using a branch-based approach with a modified Shade and Stopnisek (2019) algorithm.
}
\usage{
shade_branch(xv, nt, mxtx,fi,st)
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
\item{st}{
(Required) The number of branches to include prior to testing for increases in beta diversity. The default is to use all branches that are present in every sample and to begin testing for branches that are missing from at least one sample.
}
}
\details{
\code{shade_branch} is used internally in the holobiont package to identify core edges of a microbial phylogeny for the branch-based approach using a modified Shade and Stopnisek (2019) algorithm.
}
\value{
This function returns a list of all edges and core edges.
}
\references{

Shade, Ashley, and Nejc Stopnisek. "Abundance-occupancy distributions to prioritize plant core microbiome membership." Current opinion in microbiology 49 (2019): 50-58. Briefly, branches are ranked according to abundance first and then occupancy (i.e., occupancy is more important). Next, starting from st branches with the highest occupancy/abundance, branches are added sequentially, up to a maximum number of branches. Each time a new branch is added, the mean beta-diversity (UniFrac) between each pair of samples is calculated and is then averaged across all samples. Next, total beta-diversity (full phylogeny) between pairs samples is calculated, and again averaged across all samples. Finally, the percent increase in total beta-diversity is calculated for each new bracnh added. A threshold is then selected based on the lowest ranked branch that yields a minimum percent increase in beta-diversity.
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

shade_branch(example_phyloseq, newtree, NULL, 2, NULL)

}
