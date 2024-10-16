\name{coreFaithsPD}
\alias{coreFaithsPD}
\title{
Faiths PD of the Core Microbiome
}
\description{
Calculates Faith's phylogenetic diversity (PD) of a core microbiome based on either the tip-based or the branch-based core community phylogeny.
}
\usage{
coreFaithsPD(x, core_fraction, mode='branch', rooted=TRUE)
}
\arguments{
\item{x}{
  (Required) Microbial community data. This must be in the form of a phyloseq object and must contain, at a minimum, an OTU abundance table and a phylogeny.
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
\code{coreFaithsPD} calculates Faith's PD (Faith, 1992) based on either the tip-based or the branch-based core community phylogeny. In both cases, Faith's PD is calculated as the sum of the branch lengths in the core community phylogeny. For more details, see Bewick and Camper (2025).
}
\value{
This function returns the numeric value of Faith's PD for the core microbiome.
}

\references{
Bewick, S.A. and Benjamin T. Camper. "Phylogenetic Measures of the Core Microbiome" <doi:TBD>

Faith, Daniel P. "Conservation evaluation and phylogenetic diversity." Biological conservation 61.1 (1992): 1-10.

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

#Create a phyloseq object with a tree
example_phyloseq<-phyloseq(otu_table(enterotype),phy_tree(as.phylo(enterotype_tree)))

coreFaithsPD(example_phyloseq,0.5)

}
