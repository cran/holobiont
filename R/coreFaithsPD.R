#' @importFrom phytools bind.tip
#' @importFrom phyloseq sample_names otu_table phy_tree taxa_names
#' @importFrom ape which.edge mrca nodepath
#' @export
coreFaithsPD <- function(x,
                         core_fraction = 0.5,
                         ab_threshold1 = 0,
                         ab_threshold2 = 0,
                         ab_threshold3 = 0,
                         mode = 'branch',
                         selection = 'basic',
                         max_tax=NULL,
                         increase_cutoff = 2,
                         initial_branches=NULL,
                         rooted=TRUE) {

  core<-core_fraction

  #Using an unrooted tree in branch mode with abundance thresholds is not currently supported.  If an unrooted tree is provided with non-zero abundance thresholds and branch mode is requested, reroot it
  length_correct<-0.0001
  if (rooted == FALSE){
    length_correct<-0
    if (ab_threshold1+ab_threshold2+ab_threshold3>0){
      randoroot = sample(phy_tree(x)$tip.label, 1)
      warning("Use of an unrooted tree with non-zero abundance thresholds in branch mode is not supported.")
      warning("Randomly assigning root as -- ", randoroot, " -- in the phylogenetic tree provided.")
      phy_tree(x) <- root(phy=phy_tree(x), outgroup=randoroot, resolve.root=TRUE, interactive=FALSE)
      if( !is.rooted(phy_tree(x)) ){
        stop("Problem automatically rooting tree. Make sure your tree is rooted before attempting UniFrac calculation. See ?ape::root")
      }
      rooted = TRUE
      length_correct<-0.0001
    }
  }

  #add an outgroup so that there is the option of drawing all edges to the root rather than the minimal spanning tree
  newtree<-bind.tip(phy_tree(x),tip.label='outgroup',edge.length=0.0001,position=0)

  #if you are building a branch-based tree...
  if (mode=='branch'){

    #If using the basic selection mode...
    if (selection == 'basic'){

      #Call function to calculate all edges and core edges
      temp2<-basic_branch(x,newtree,core,ab_threshold1,ab_threshold2,ab_threshold3,rooted)

      #If using the shade selection mode...
    }else if (selection == 'shade'){

      #If no maximum taxon is specified, default to considering 1%
      if (is.null(max_tax)){max_tax<-1}

      #Call the function to calculate all edges and core edges
      temp2<-shade_branch(x,newtree,max_tax,increase_cutoff,initial_branches)

      #If a non-supported mode is selected print a warning and stop
    }else{

      stop('That method of core selection is not available. Please use basic or shade.')

    }

    #make a list of the core edges from the focal habitat
    core_edges<-temp2$core_edges


    #sum up the lengths of those edges, possibly making a corection for the outgroup length
    core_length<-sum(newtree$edge.length[core_edges])-length_correct

  #if you are building a tip-based tree...
  }else if (mode=='tip'){

    #If using the basic selection mode...
    if (selection == 'basic'){

      #Call function to calculate all edges and core edges
      temp2<-basic_tip(x,newtree,core,ab_threshold1,ab_threshold2,ab_threshold3,rooted)

      #If using the shade selection mode...
    }else if (selection == 'shade'){

      #If no maximum taxon is specified, default to considering 1%
      if (is.null(max_tax)){max_tax<-1}

      #Call the function to calculate all edges and core edges
      temp2<-shade_tip(x,newtree,max_tax,increase_cutoff)

      #If a non-supported mode is selected print a warning and stop
    }else{

      stop('That method of core selection is not available. Please use basic or shade.')

    }

    #make a list of the core edges from the focal habitat
    core_edges<-temp2$core_edges

    #sum up the lengths of those edges, possibly making a corection for the outgroup length
    core_length<-sum(newtree$edge.length[core_edges])-length_correct

    }

  #if a mode was entered that is not supported, print a warning
  else{warning('Warning: that mode is not supported')}

  #return the summed length of the edges
  return(core_length)
}
