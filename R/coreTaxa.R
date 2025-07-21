#' @importFrom phytools bind.tip
#' @importFrom phyloseq sample_names otu_table phy_tree prune_taxa taxa_names
#' @importFrom ape which.edge mrca nodepath plot.phylo chronos
#' @importFrom castor get_all_distances_to_root
#' @export
coreTaxa <- function(x,
                      core_fraction = 0.5,
                      ab_threshold1 = 0,
                      ab_threshold2 = 0,
                      ab_threshold3 = 0,
                      selection = 'basic',
                      max_tax=NULL,
                      increase_cutoff = 2,
                      remove_zeros=TRUE) {

  core<-core_fraction

  if (selection == 'basic'){

    #Call function to calculate all edges and core edges
    temp2<-basic_np(x,core,ab_threshold1,ab_threshold2,ab_threshold3)

  }else if (selection == 'shade'){

    #If no maximum taxon is specified, default to considering 1%
    if (is.null(max_tax)){max_tax<-1}
    temp2<-shade_np(x,max_tax,increase_cutoff)

  }else{

    stop('That method of core selection is not available. Please use basic or shade.')

  }

  return(temp2$core_taxa)

}
