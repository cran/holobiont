\name{coreTree}
\alias{coreTree}
\title{
The Core Community Phylogeny
}
\description{
Identifies and plots the tip-based or the branch-based core community phylogeny based on the occurrence of abundance of different microbial lineages in a set of samples from a common habitat (e.g., type of host or environment).
}
\usage{
coreTree(x, core_fraction = 0.5, ab_threshold1 = 0,
ab_threshold2 = 0,ab_threshold3 = 0, mode='branch',
selection = 'basic',max_tax=NULL, increase_cutoff = 2,
initial_branches = NULL, NCcol = 'black', Ccol='red',rooted=TRUE,
branch.width=4,label.tips=FALSE, scaled = FALSE, remove_zeros=TRUE,
plot.chronogram=FALSE)
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
}\item{NCcol}{
The color to plot all branches of the phylogeny that are NOT part of the core community phylogeny. The default is black.
}
\item{Ccol}{
The color to plot all branches of the phylogeny that are ARE part of the core community phylogeny. The default is red.
}
\item{branch.width}{
The width to use when plotting the branches of the phylogeny. The default is 4.
}
\item{label.tips}{
Whether or not to label the tips of the phylogeny with the microbial taxon names. The default is FALSE.
}
\item{scaled}{
Whether or not to scale the branch lengths. The default is FALSE.
}
\item{remove_zeros}{
Whether or not to remove taxa that are missing from all samples prior to drawing the phylogeny. The default is TRUE.
}
\item{plot.chronogram}{
Whether to plot a phylogeny or a chronogram. The default is FALSE.
}

}
\details{
\code{coreTree} identifies either the tip-based or the branch-based core community phylogeny. For the tip-based core community phylogeny, individual microbial taxa are retained based on being present in a threshold number of samples or at a threshold abundance. Once core microbial taxa have been identified, they are used to reconstruct the core community phylogeny. For the branch-based core community phylogeny, the phylogenetic tree for the entire dataset is examined, branch-by-branch, to determine which branches should be retained based on being present in a threshold number of samples or at a threshold abundance. If \code{rooted = TRUE}, branches are counted based on individual sample phylogenies that include the root node. Likewise, the tip-based tree is forced to include the root. If \code{rooted = FALSE}, branches are counted based on individual sample phylogenies that span all taxa present in the sample. Similarly, the tip-based phylogeny is the tree that spans all core taxa, and may not include the root.  For more details, see Bewick and Camper (2025).
}
\value{
This function plots the phylogeny for the entire dataset in black and colors the branches that are part of the core community phylogeny in red. These colors can be altered using the \code{NCcol} and \code{Ccol} variables.
}

\references{
Bewick, S.A. and Benjamin T. Camper. "Phylogenetic Measures of the Core Microbiome" <doi:TBD>

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

coreTree(example_phyloseq,0.5)

}
