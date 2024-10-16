#' @importFrom phytools bind.tip
#' @importFrom phyloseq sample_names otu_table phy_tree prune_taxa taxa_names
#' @importFrom ape which.edge mrca nodepath plot.phylo chronos
#' @export
coreTree <- function(x, core_fraction, mode = 'branch',NCcol = 'black',Ccol='red',rooted=TRUE,branch.width=4,label.tips=FALSE, remove_zeros=TRUE,plot.chronogram=FALSE) {

  core<-core_fraction

  #remove taxa that are not present in any sample
  if (remove_zeros==TRUE){
    x<-prune_taxa(taxa_names(x)[which(rowSums(sign(otu_table(x)))>0)],x)
  }

  #add an outgroup so that there is the option of drawing all edges to the root rather than the minimal spanning tree
  newtree<-bind.tip(phy_tree(x),tip.label='outgroup',edge.length=0.0001,position=0)

  #if you are building a branch-based tree...
  if (mode=='branch'){

    #initialize a list of edges that are present in the samples
    edges<-c()

    #for each sample...
    for (i in 1:length(sample_names(x))){

      if (rooted == TRUE){
      #find the taxa in the sample (including the outgroup so that you draw the tree back to the root)
        nz<-c('outgroup',taxa_names(x)[which(otu_table(x)[,i]>0)])
      }
      else{
        nz<-taxa_names(x)[which(otu_table(x)[,i]>0)]
      }
      #find the edges associated with those taxa (drawn back to the tree root)
      edges<-c(edges,which.edge(newtree,nz))
    }


    #find counts of the number of times each edge appeared across all the samples
    branch_counts<-table(edges)


    #pull out the edges that were present in at least a core threshold number of samples
    core_branch<-which(branch_counts>=core*length(sample_names(x)))

    #make a list of the core edges
    core_edges<-as.integer(names(core_branch))

    #if the root is not included
    if (rooted==FALSE){

      #find the nodes associated with core edges
      nodes<-unique(c(newtree$edge[,1][core_edges],newtree$edge[,2][core_edges]))
      #find the mrca of each node in the tree
      cc<-mrca(newtree,full=TRUE)
      #find the mrca of each node associated with a core edge
      mrca_matrix<-cc[nodes,nodes]
      #find the unique mrcas for the core edge nodes
      mrca_list<-unique(as.vector(mrca_matrix))
      #find the unique mrcas plus core edge nodes
      mrca_list<-unique(mrca_list,nodes)
      #identify mrcas missing from the list of nodes associated with core edges
      missing<-mrca_list[which(!(mrca_list %in% nodes))]
      if (length(missing)>0){
      for (i in 1:length(missing)){
        for (j in 1:length(nodes)){
          #find the nodes connecting the missing mrcas to the nodes associated with core edges
          mrca_list<-c(mrca_list,nodepath(newtree,from=missing[i],to=nodes[j]))
        }
      }
      }
      #find the edges associated with all the nodes (core and mrcas)
      all_core_edges<-intersect(which(newtree$edge[,1] %in% mrca_list),which(newtree$edge[,2] %in% mrca_list))
      #add the missing edges to the core edges
      core_edges<-unique(c(core_edges,all_core_edges))
      }


  #if you are building a tip-based tree
  }else if (mode=='tip'){

    #find all of the edges in the tree (not generated automatically as it was for branch-based tree)
    alltaxa<-taxa_names(x)
    edges<-which.edge(newtree,alltaxa)

    #if you want to include the root in the tip-based tree
    if (rooted==TRUE){
      #find all of the taxa present in a threshold number of samples; add the outgroup so that edges are drawn back to the root
      coretaxa<-c('outgroup',taxa_names(x)[which(rowSums(sign(otu_table(x)))>=core*length(sample_names(x)))])
    }else{
      #find all of the taxa present in a threshold number of samples; do not add the outgroup so that edges are drawn on the minimal spanning tree
      coretaxa<-taxa_names(x)[which(rowSums(sign(otu_table(x)))>=core*length(sample_names(x)))]
    }

    #find the edges associated with all of the core taxa identified above
    core_edges<-which.edge(newtree,coretaxa)
  }

  #if the mode doesn't match one of the options, print a warning
  else{warning('Warning, that mode is not supported')}

  #make a state vector to designate which edges are core and which are not; label all edges with 1 (not-core) initially
  states<-rep(1,length(edges))
  #label all of the core edges with 2
  states[core_edges]<-2

  #color code the states 1 and 2
  branch_col<-c(NCcol,Ccol)

  #plot the tree (omitting the outgroup... since it was added at the root, it is the last edge and can be ignored)
  if (plot.chronogram==FALSE){
    tt<-plot.phylo(phy_tree(x),edge.color=branch_col[states],edge.width=branch.width,show.tip.label=label.tips )
  }else{
    tt<-plot.phylo(chronos(phy_tree(x)),edge.color=branch_col[states],edge.width=branch.width,show.tip.label=label.tips )
  }
}
