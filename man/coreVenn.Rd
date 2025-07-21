\name{coreVenn}
\alias{coreVenn}
\title{
Venn Diagram Comparing (Non-Phylogenetic) Core Microbiomes from Multiple Habitats
}
\description{
Generates a Venn diagram comparing and contrasting the shared and unique taxa of core microbiomes from two or more (up to seven) different types of habitats using a non-phylogenetic approach.
}
\usage{
coreVenn(x, grouping, core_fraction = 0.5, ab_threshold1 = 0,
ab_threshold2 = 0, ab_threshold3 = 0, selection='basic', max_tax = NULL,
increase_cutoff = 2,ordered_groups=NULL,show_percentage=TRUE,decimal=2,
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
\code{coreVenn} generates a Venn diagram showing either the percentage of taxa or the taxon count shared between the core communities from all possible combinations of up to seven different habitat types. For more details, see Bewick and Camper (2025).
}
\value{
This function returns a plot of the Venn diagram.
}

\references{
Bewick, S.A. and Benjamin T. Camper. "Phylogenetic Measures of the Core Microbiome" <doi:TBD>

McMurdie, Paul J., and Susan Holmes. "phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data." PloS one 8.4 (2013): e61217.

McMurdie, Paul J., and Susan Holmes. "Phyloseq: a bioconductor package for handling and analysis of high-throughput phylogenetic sequence data." Biocomputing 2012. 2012. 235-246.

Shade, Ashley, and Nejc Stopnisek. "Abundance-occupancy distributions to prioritize plant core microbiome membership." Current opinion in microbiology 49 (2019): 50-58.

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
enterotype_tree$node.label<-as.character(1:1:enterotype_tree$Nnode)

#keep only those samples with gender identified
gendered<-which(!(is.na(sample_data(enterotype)$Gender)))
enterotypeMF<-prune_samples(sample_names(enterotype)[gendered],enterotype)

#Create a phyloseq object with a tree
example_phyloseq<-phyloseq(otu_table(enterotypeMF),phy_tree(as.phylo(enterotype_tree)))

coreVenn(example_phyloseq,grouping=sample_data(enterotypeMF)$Gender)

}
