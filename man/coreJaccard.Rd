\name{coreJaccard}
\alias{coreJaccard}
\title{
Jaccard Distance Between Core Microbiomes from Two Habitats
}
\description{
Calculates the Jaccard distance between the core microbiomes from two different types of habitats based on either the tip-based or the branch-based core community phylogeny.
}
\usage{
coreJaccard(x, grouping, core_fraction)
}
\arguments{
  \item{x}{
  (Required) Microbial community data. This must be in the form of a phyloseq object and must contain, at a minimum, an OTU abundance table.
}
\item{grouping}{
(Required) A vector specifying which samples belong to which habitat type.
}
\item{core_fraction}{
(Required) The fraction of samples that a microbial taxon must be found in to be considered part of the 'core' microbiome.
}

}
\details{
\code{coreJaccard} calculates the Jaccard distance (Jaccard, 1901 A and B) between the core microbiomes from two different types of habitats. Briefly, the Jaccard distance is calculated as the number of unique taxa in the core communities from each of the two habitats, divided by the total number of taxa in the two habitats combined. For more details, see Bewick and Camper (2025).
}
\value{
This function returns the numeric value of the Jaccard distance between two core microbiomes.
}

\references{
Bewick, S.A. and Benjamin T. Camper. "Phylogenetic Measures of the Core Microbiome" <doi:TBD>

Jaccard, P. (1901A) Étude comparative de la distribution florale dans une portion des Alpes et des Jura. Bulletin de la Société Vaudoise des Sciences Naturelles 37, 547-579.

Jaccard, P. (1901B) Distribution de la flore alpine dans le bassin des Dranses et dans quelques régions voisines. Bulletin de la Société Vaudoise des Sciences Naturelles 37, 241-272.

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
example_phyloseq<-phyloseq(otu_table(enterotypeMF))

coreJaccard(example_phyloseq,grouping=sample_data(enterotypeMF)$Gender,0.5)

}
