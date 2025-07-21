#' @importFrom phytools bind.tip
#' @importFrom phyloseq sample_names otu_table phy_tree taxa_names
#' @importFrom ape which.edge mrca nodepath
#' @export
coreUniFrac <- function(x,
                        grouping,
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

  #find the number of different habitat types (e.g. hosts or environments) that are being compared
  group_count<-length(unique(grouping))
  group_id<-unique(grouping)

  #make a list of the samples from the first habitat
  group1<-sample_names(x)[which(grouping==group_id[1])]
  #make a list of the samples from the second habitat
  group2<-sample_names(x)[which(grouping==group_id[2])]

  x1<-prune_samples(group1,x)
  x2<-prune_samples(group2,x)

  #add an outgroup so that there is the option of drawing all edges to the root rather than the minimal spanning tree
  newtree<-bind.tip(phy_tree(x),tip.label='outgroup',edge.length=0.0001,position=0)

  #if you are building a branch-based tree...
  if (mode=='branch'){

    #If using the basic selection mode...
    if (selection == 'basic'){

      temp1<-basic_branch(x1,newtree,core,ab_threshold1,ab_threshold2,ab_threshold3,rooted)
      temp2<-basic_branch(x2,newtree,core,ab_threshold1,ab_threshold2,ab_threshold3,rooted)

    }else if (selection == 'shade'){

      #If no maximum taxon is specified, default to considering 1%
      if (is.null(max_tax)){max_tax<-1}

      #Call the function to calculate all edges and core edges
      temp1<-shade_branch(x1,newtree,max_tax,increase_cutoff,initial_branches)
      temp2<-shade_branch(x2,newtree,max_tax,increase_cutoff,initial_branches)

       #If a non-supported mode is selected print a warning and stop
    }else{

      stop('That method of core selection is not available. Please use basic or shade.')

    }


  }else if (mode=='tip'){

    #If using the basic selection mode...
    if (selection == 'basic'){

      #Call function to calculate all edges and core edges
      temp1<-basic_tip(x1,newtree,core,ab_threshold1,ab_threshold2,ab_threshold3,rooted)
      temp2<-basic_tip(x2,newtree,core,ab_threshold1,ab_threshold2,ab_threshold3,rooted)

      #If using the shade selection mode...
    }else if (selection == 'shade'){

      #If no maximum taxon is specified, default to considering 1%
      if (is.null(max_tax)){max_tax<-1}

      #Call the function to calculate all edges and core edges
      temp1<-shade_tip(x1,newtree,max_tax,increase_cutoff)
      temp2<-shade_tip(x2,newtree,max_tax,increase_cutoff)

      #If a non-supported mode is selected print a warning and stop
    }else{

      stop('That method of core selection is not available. Please use basic or shade.')

    }
  }

  core_edges1<-temp1$core_edges
  core_edges2<-temp2$core_edges



    #find the edges that are core to the first habitat but not the second
    in1only<-which(!(core_edges1 %in% core_edges2))
    #find the edges that are core to the second habitat but not the first
    in2only<-which(!(core_edges2 %in% core_edges1))

    #divide the length of the unique core edges by the total length or the core edges
    if (length(in1only)>0){
      A<-sum(newtree$edge.length[core_edges1[in1only]])
    }else{A<-0}
    if (length(in2only)>0){
      B<-sum(newtree$edge.length[core_edges2[in2only]])
    }else{B<-0}

    unifrac<-(A+B)/(sum(newtree$edge.length[unique(c(core_edges1,core_edges2))])-length_correct)



  #return the unifrac distance
  return(unifrac)

}
