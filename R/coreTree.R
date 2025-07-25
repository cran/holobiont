#' @importFrom phytools bind.tip
#' @importFrom phyloseq sample_names otu_table phy_tree prune_taxa taxa_names
#' @importFrom ape which.edge mrca nodepath plot.phylo chronos
#' @importFrom castor get_all_distances_to_root
#' @export
coreTree <- function(x, core_fraction = 0.5, ab_threshold1 = 0, ab_threshold2 = 0, ab_threshold3 = 0, mode = 'branch', selection = 'basic', max_tax=NULL, increase_cutoff = 2, initial_branches=NULL,NCcol = 'black',Ccol='red',rooted=TRUE,branch.width=4,label.tips=FALSE, scaled = FALSE, remove_zeros=TRUE,plot.chronogram=FALSE) {

  core<-core_fraction

  #remove taxa that are not present in any sample
  if (remove_zeros==TRUE){
    x<-prune_taxa(taxa_names(x)[which(rowSums(sign(otu_table(x)))>0)],x)
  }

  #Using an unrooted tree in branch mode with abundance thresholds is not currently supported.  If an unrooted tree is provided with non-zero abundance thresholds and branch mode is requested, reroot it
  if (rooted == FALSE){
    if (ab_threshold1+ab_threshold2+ab_threshold3>0){
      randoroot = sample(phy_tree(x)$tip.label, 1)
      warning("Use of an unrooted tree with non-zero abundance thresholds in branch mode is not supported.")
      warning("Randomly assigning root as -- ", randoroot, " -- in the phylogenetic tree provided.")
      phy_tree(x) <- root(phy=phy_tree(x), outgroup=randoroot, resolve.root=TRUE, interactive=FALSE)
      if( !is.rooted(phy_tree(x)) ){
        stop("Problem automatically rooting tree. Make sure your tree is rooted before attempting UniFrac calculation. See ?ape::root")
      }
      rooted = TRUE
    }
  }

  #add an outgroup so that there is the option of drawing all edges to the root rather than the minimal spanning tree
  newtree<-bind.tip(phy_tree(x),tip.label='outgroup',edge.length=0.0001,position=0)

  #if you are building a branch-based tree...
  if (mode=='branch'){

    if (selection == 'basic'){

    #Call function to calculate all edges and core edges
      temp<-basic_branch(x,newtree,core,ab_threshold1,ab_threshold2,ab_threshold3,rooted)

    }else if (selection == 'shade'){

      #If no maximum taxon is specified, default to considering 1%
      if (is.null(max_tax)){max_tax<-round(length(newtree$tip.label)*0.01)}
      temp<-shade_branch(x,newtree,max_tax,increase_cutoff,initial_branches)

    }else{

      warning('That method of core selection is not available. Please use basic or shade.')

    }
    edges<-temp$edges
    core_edges<-temp$core_edges

    #if you are building a tip-based tree
  }else if (mode=='tip'){

    if (selection == 'basic'){

      #Call function to calculate all edges and core edges
      temp<-basic_tip(x,newtree,core,ab_threshold1,ab_threshold2,ab_threshold3,rooted)

    }else if (selection == 'shade'){

      #If no maximum taxon is specified, default to considering 1%
      if (is.null(max_tax)){max_tax<-round(length(newtree$tip.label)*0.01)}
      temp<-shade_tip(x,newtree,max_tax,increase_cutoff)
    }else{

      warning('That method of core selection is not available. Please use basic or shade.')

    }

    edges<-temp$edges
    core_edges<-temp$core_edges

  }

  #if the mode doesn't match one of the options, print a warning
  else{warning('Warning, that mode is not supported')}

  #make a state vector to designate which edges are core and which are not; label all edges with 1 (not-core) initially
  states<-rep(1,length(phy_tree(x)$edge.length))
  #label all of the core edges with 2
  states[core_edges]<-2

  #color code the states 1 and 2
  branch_col<-c(NCcol,Ccol)

  #plot the tree (omitting the outgroup... since it was added at the root, it is the last edge and can be ignored)
  if (plot.chronogram==FALSE){
    #plot the tree
    if (scaled == FALSE){
      tt<-plot.phylo(phy_tree(x),edge.color=branch_col[states],edge.width = branch.width,show.tip.label=label.tips,cex=0.5)
    }else{
      temptree<-phy_tree(x)
      temptree$edge.length <- sqrt(temptree$edge.length+1e-2)
      tt<-plot.phylo(temptree,edge.color=branch_col[states],edge.width = branch.width,show.tip.label=label.tips,cex=0.5)
    }
  }else{
    tt<-plot.phylo(chronos(phy_tree(x)),edge.color=branch_col[states],edge.width=branch.width,show.tip.label=label.tips,cex = 0.5 )
  }
}
