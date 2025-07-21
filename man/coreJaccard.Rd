\name{coreJaccard}
\alias{coreJaccard}
\title{
Jaccard Distance Between Core Microbiomes from Two Habitats
}
\description{
Calculates the Jaccard distance between the core microbiomes from two different types of habitats based on either the tip-based or the branch-based core community phylogeny.
}
\usage{
coreJaccard(x, grouping, core_fraction = 0.5, ab_threshold1 = 0,
ab_threshold2 = 0, ab_threshold3 = 0, selection='basic',
max_tax = NULL, increase_cutoff = 2)
}
\arguments{
  \item{x}{
  (Required) Microbial community data. This must be in the form of a phyloseq object and must contain, at a minimum, an OTU abundance table.
}
\item{grouping}{
(Required) A vector specifying which samples belong to which habitat type.
}
\item{core_fraction}{
The fraction of samples that a microbial taxon must be found in to be considered part of the 'core' microbiome. This variable is only used when selection = 'basic' and is ignored when selection = 'shade'. The default value is 0.5.
}
\item{ab_threshold1}{
The threshold for mean relative abundance across all samples. This variable is only used when selection = 'basic' and is ignored when selection = 'shade'. The default value is 0.
}
\item{ab_threshold2}{
The threshold for maximum relative abundance in any sample. This variable is only used when selection = 'basic' and is ignored when selection = 'shade'. The default value is 0.
}
\item{ab_threshold3}{
The threshold for minimum relative abundance across all samples. This variable is only used when selection = 'basic' and is ignored when selection = 'shade'. The default value is 0.
}
\item{selection}{
Whether to use thresholds ('basic') or the Shade and Stopnisek method ('shade') to define the core community. The default is 'basic'.
}
\item{max_tax}{
The maximum number of branches to add sequentially, as a percentage of the total branches when using the Shade and Stopnisek method. This variable is only used when selection = 'shade' and is ignored when selection = 'basic'.
}
\item{increase_cutoff}{
The threshold for the percent increase in beta diversity used to identify the taxon at which point adding more taxa yields diminishing returns in explanatory power. This variable is only used when selection = 'shade' and is ignored when selection = 'basic'.
}

}
\details{
\code{coreJaccard} calculates the Jaccard distance (Jaccard, 1901 A and B) between the core microbiomes from two different types of habitats using either basic thresholds or a modification of the Shade and Stopnisek (2019) algorithm. Briefly, the Jaccard distance is calculated as the number of unique taxa in the core communities from each of the two habitats, divided by the total number of taxa in the two habitats combined. For more details, see Bewick and Camper (2025).
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

Shade, Ashley, and Nejc Stopnisek. "Abundance-occupancy distributions to prioritize plant core microbiome membership." Current opinion in microbiology 49 (2019): 50-58.
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

coreJaccard(example_phyloseq,grouping=sample_data(enterotypeMF)$Gender)

}
