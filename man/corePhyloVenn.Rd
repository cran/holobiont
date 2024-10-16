\name{corePhyloVenn}
\alias{corePhyloVenn}
\title{
Phylogenetic Venn Diagram Comparing Core Microbiomes from Multiple Habitats
}
\description{
Generates a Venn diagram comparing and contrasting the shared and unique phylogenetic branch lengths of core microbiomes from two or more (up to seven) different types of habitats based on either the tip-based or the branch-based core community phylogeny.
}
\usage{
corePhyloVenn(x, grouping, core_fraction,
mode = 'branch',rooted=TRUE, ordered_groups=NULL,show_percentage=TRUE,decimal=2,
fill_color=c('red','orange','yellow','green','blue','purple','black'),fill_alpha=0.5,
stroke_color='black',stroke_alpha = 1, stroke_size = 1,stroke_linetype = "solid",
set_name_color = "black", set_name_size = 6,text_color = "black",text_size = 4)
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
\item{ordered_groups}{
When provided, specifies the order in which different habitats should be plotted. This order is matched to the color order specified in \code{fill_color} (see below).
}
\item{show_percentage}{
If \code{TRUE} the Venn diagram is shown as percentages of the total branch length in each compartment. If \code{FALSE} the Venn diagram is shown as absolute branch lengths.
}
\item{decimal}{
The number of decimal points to report in the Venn diagram compartments (either for percentages or absolute branch lengths).
}
\item{fill_color}{
A vector specifying the colors to use for each of the different habitats in the Venn diagram. The default is to use the six colors of the rainbow, followed by black.
}
\item{fill_alpha}{
The transparency of the fill colors in the Venn diagram.
}
\item{stroke_color}{
The color of the outlines in the Venn diagram.
}
\item{stroke_alpha}{
The transparency of the outlines in the Venn diagram.
}
\item{stroke_size}{
The width of the outlines in the Venn diagram.
}
\item{stroke_linetype}{
The style of the outlines in the Venn diagram.
}
\item{set_name_color}{
The color of the font used to label each habitat.
}
\item{set_name_size}{
The size of the font used to label each habitat
}
\item{text_color}{
The color of the font used to report the value in each compartment of the Venn diagram
}
\item{text_size}{
The size of the font used to report the value in each compartment of the Venn diagram
}

}
\details{
\code{corePhyloVenn} generates a Venn diagram showing either the percentage of the phylogenetic branch length or the absolute phylogenetic branch length shared between the core community phylogenies from all possible combinations of up to seven different habitat types. For more details, see Bewick and Camper (2025).
}
\value{
This function returns a plot of the phylogenetic Venn diagram.
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
library(ggplot2)
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

corePhyloVenn(example_phyloseq,grouping=sample_data(enterotypeMF)$Gender,0.5)

}
