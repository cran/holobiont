\name{coreUniFrac}
\alias{coreUniFrac}
\title{
UniFrac Distance Between Core Microbiomes from Two Habitats
}
\description{
Calculates the UniFrac distance between the core microbiomes from two different types of habitats based on either the tip-based or the branch-based core community phylogeny.
}
\usage{
coreUniFrac(x, grouping, core_fraction, mode='branch', rooted=TRUE)
}
\arguments{
  \item{x}{
  (Required) Microbial community data. This must be in the form of a phyloseq object and must contain, at a minimum, an OTU abundance table and a phylogeny.
}
\item{grouping}{
(Required) A vector specifying which samples belong to which habitat type.
}
\item{core_fraction}{
(Required) The fraction of samples that a microbial taxon must be found in to be considered part of the 'core' microbiome.
}
\item{mode}{
Whether to build a tip-based ('tip') or a branch-based ('branch') phylogeny. The default is 'branch'.
}
\item{rooted}{
Whether to include the root of the phylogeny. The default is TRUE, meaning that the root is necessarily included in all phylogenies. This requires that the input tree be rooted.
}

}
\details{
\code{coreUniFrac} calculates the UniFrac distance (Lozupone and Knight, 2005) between the core microbiomes from two different types of habitats based on either their tip-based or their branch-based core community phylogenies. In both cases, the UniFrac distance is calculated as the sum of the unique branch lengths in the core community phylogenies from each of the two habitats, divided by the total branch length of all branches in the core community phylogenies from the two habitats combined. For more details, see Bewick and Camper (2025).
}
\value{
This function returns the numeric value of the UniFrac distance between two core microbiomes.
}

\references{
Bewick, S.A. and Benjamin T. Camper. "Phylogenetic Measures of the Core Microbiome" <doi:TBD>

Lozupone, Catherine, and Rob Knight. "UniFrac: a new phylogenetic method for comparing microbial communities." Applied and environmental microbiology 71.12 (2005): 8228-8235.

McMurdie, Paul J., and Susan Holmes. "phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data." PloS one 8.4 (2013): e61217.

McMurdie, Paul J., and Susan Holmes. "Phyloseq: a bioconductor package for handling and analysis of high-throughput phylogenetic sequence data." Biocomputing 2012. 2012. 235-246.
}

\examples{
#Test with enterotype dataset
library(phyloseq)
library(ape)
library(phytools)
data(enterotype)

set.seed(1)

#Generate an example tree and label it with the names of the microbial taxa
enterotype_tree<-rtree(length(taxa_names(enterotype)))
enterotype_tree$tip.label<-taxa_names(enterotype)

#keep only those samples with gender identified
gendered<-which(!(is.na(sample_data(enterotype)$Gender)))
enterotypeMF<-prune_samples(sample_names(enterotype)[gendered],enterotype)

#Create a phyloseq object with a tree
example_phyloseq<-phyloseq(otu_table(enterotypeMF),phy_tree(as.phylo(enterotype_tree)))

coreUniFrac(example_phyloseq,grouping=sample_data(enterotypeMF)$Gender,0.5)

}
