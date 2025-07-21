\name{coreUniFrac}
\alias{coreUniFrac}
\title{
UniFrac Distance Between Core Microbiomes from Two Habitats
}
\description{
Calculates the UniFrac distance between the core microbiomes from two different types of habitats based on either the tip-based or the branch-based core community phylogeny.
}
\usage{
coreUniFrac(x, grouping, core_fraction = 0.5, ab_threshold1 = 0,
ab_threshold2 = 0, ab_threshold3 = 0, mode='branch', selection='basic',
max_tax = NULL, increase_cutoff = 2, initial_branches = NULL, rooted=TRUE)
}
\arguments{
  \item{x}{
  (Required) Microbial community data. This must be in the form of a phyloseq object and must contain, at a minimum, an OTU abundance table and a phylogeny.
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
\item{mode}{
Whether to build a tip-based ('tip') or a branch-based ('branch') phylogeny. The default is 'branch'.
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
\item{initial_branches}{
The number of branches to include prior to testing for increases in beta diversity. The default is to use all branches that are present in every sample and to begin testing branches that are missing from at least one sample. This variable is only used when selection = 'shade' and is ignored when selection = 'basic'.
}
\item{rooted}{
Whether to include the root of the phylogeny. The default is TRUE, meaning that the root is necessarily included in all phylogenies. This requires that the input tree be rooted.
}

}
\details{
\code{coreUniFrac} calculates the UniFrac distance (Lozupone and Knight, 2005) between the core microbiomes from two different types of habitats based on either their tip-based or their branch-based core community phylogenies and using either basic thresholds or a modification of the Shade and Stopnisek (2019) algorithm. In both cases, the UniFrac distance is calculated as the sum of the unique branch lengths in the core community phylogenies from each of the two habitats, divided by the total branch length of all branches in the core community phylogenies from the two habitats combined. For more details, see Bewick and Camper (2025).
}
\value{
This function returns the numeric value of the UniFrac distance between two core microbiomes.
}

\references{
Bewick, S.A. and Benjamin T. Camper. "Phylogenetic Measures of the Core Microbiome" <doi:TBD>

Lozupone, Catherine, and Rob Knight. "UniFrac: a new phylogenetic method for comparing microbial communities." Applied and environmental microbiology 71.12 (2005): 8228-8235.

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
enterotype_tree$node.label<-as.character(1:1:enterotype_tree$Nnode)

#keep only those samples with gender identified
gendered<-which(!(is.na(sample_data(enterotype)$Gender)))
enterotypeMF<-prune_samples(sample_names(enterotype)[gendered],enterotype)

#Create a phyloseq object with a tree
example_phyloseq<-phyloseq(otu_table(enterotypeMF),phy_tree(as.phylo(enterotype_tree)))

coreUniFrac(example_phyloseq,grouping=sample_data(enterotypeMF)$Gender)

}
