\name{basic_tip}
\alias{basic_tip}
\title{
Tip-based Core Community Phylogeny Using Thresholds
}
\description{
Called internally. Identifies all edges and core edges of the core community phylogeny using a tip-based approach with core thresholds.
}
\usage{
basic_tip(xv, nt, cf, abt1, abt2, abt3,rt)
}
\arguments{
\item{xv}{
  (Required) The phyloseq object passed from main functions and containing microbial community data.
}
\item{nt}{
(Required) The microbial phylogeny passed from main functions.
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
\item{rt}{
Whether to include the root of the phylogeny. The default is TRUE, meaning that the root is necessarily included in all phylogenies. This requires that the input tree be rooted.
}
}
\details{
\code{basic_tip} is used internally in the holobiont package to identify core edges of a microbial phylogeny for the tip-based approach using thresholds.
}
\value{
This function returns a list of all edges and core edges, as well as all taxa and core taxa.
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

basic_tip(example_phyloseq, phy_tree(example_phyloseq), 0.5, 0, 0, 0,TRUE)


}
