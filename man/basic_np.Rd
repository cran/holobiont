\name{basic_np}
\alias{basic_np}
\title{
Non-phylogenetic Core Community Phylogeny Using Thresholds
}
\description{
Called internally. Identifies all edges and core edges of the core community phylogeny using a non-phylogenetic approach with core thresholds.
}
\usage{
basic_np(xv, cf, abt1, abt2, abt3)
}
\arguments{
\item{xv}{
  (Required) The phyloseq object passed from main functions and containing microbial community data.
}
\item{cf}{
(Required) The fraction of samples that a microbial taxon must be found in to be considered part of the 'core' microbiome.
}
\item{abt1}{
(Required) The threshold for mean abundance across all samples.
}
\item{abt2}{
(Required) The threshold for maximum abundance in any sample.
}
\item{abt3}{
(Required) The threshold for the minimum abundance across all samples.
}
}
\details{
\code{basic_np} is used internally in the holobiont package to identify all taxa and all core taxa for the non-phylogenetic approach using thresholds.
}
\value{
This function returns a list of all taxa and all core taxa.
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

basic_np(example_phyloseq, 0.5, 0, 0, 0)

}
