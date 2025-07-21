\name{coreRichness}
\alias{coreRichness}
\title{
Richness of the Core Microbiome
}
\description{
Calculates richness of a core microbiome using a non-phylogenetic approach.
}
\usage{
coreRichness(x, core_fraction = 0.5, ab_threshold1 = 0, ab_threshold2 = 0,
ab_threshold3 = 0, selection='basic', max_tax = NULL, increase_cutoff = 2)
}
\arguments{
\item{x}{
  (Required) Microbial community data. This must be in the form of a phyloseq object and must contain, at a minimum, an OTU abundance table and a phylogeny.
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
\code{coreRichness} calculates richness (a count of microbial taxa) using either basic thresholds or a modification of the Shade and Stopnisek (2019) algorithm.
}
\value{
This function returns the numeric value of richness for the core microbiome.
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

#Generate an example tree and label it with the names of the microbial taxa
enterotype_tree<-rtree(length(taxa_names(enterotype)))
enterotype_tree$tip.label<-taxa_names(enterotype)

#Create a phyloseq object with a tree
example_phyloseq<-phyloseq(otu_table(enterotype),phy_tree(as.phylo(enterotype_tree)))

coreRichness(example_phyloseq)

}
