#' @importFrom phytools bind.tip
#' @importFrom phyloseq sample_names otu_table phy_tree taxa_names
#' @importFrom ape which.edge mrca nodepath
#' @export
coreFaithsPD <- function(x, core_fraction, mode = 'branch',rooted=TRUE) {

  core<-core_fraction

  #add an outgroup so that there is the option of drawing all edges to the root rather than the minimal spanning tree
  newtree<-bind.tip(phy_tree(x),tip.label='outgroup',edge.length=0.0001,position=0)

  #if you are building a branch-based tree...
  if (mode=='branch'){

    #initialize a list of all the edges associated with taxa present in the system
    edges<-c()

    #for each sample...
    for (i in 1:length(sample_names(x))){

      #if you are including the root...
      if (rooted==TRUE){

        #find the taxa present in the sample; include the outgroup so that you find a tree spanning the root; account for core = 0 by including taxa not present in any samples
        if (core>0){
          nz<-c('outgroup',taxa_names(x)[which(otu_table(x)[,i]>0)])
        }else{
          nz<-c('outgroup',taxa_names(x)[which(otu_table(x)[,i]>=0)])
        }
        #if you include the outgroup, put in a correction for its length
        length_correct<-0.0001

      #if you are not including the root
      }else{

        #find the taxa present in the sample; do NOT include the outgroup because you are not including the root; account for core = 0 by including taxa not present in any samples
        if (core>0){
          nz<-taxa_names(x)[which(otu_table(x)[,i]>0)]
        }else{
          nz<-taxa_names(x)[which(otu_table(x)[,i]>=0)]
        }
        #if you did not include the outgroup, there is no correction for its length
        length_correct<-0
      }

      #find the edges associated with those taxa
      edges<-c(edges,which.edge(newtree,nz))
    }

    #find counts of the number of times each edge appeared across all samples
    branch_counts<-table(edges)

    #pull out the core edges that were present in at least a core number of samples
    core_branch<-which(branch_counts>=core*length(sample_names(x)))

    #make a list of the core edges for inornatus
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


    #sum up the lengths of those edges, possibly making a corection for the outgroup length
    core_length<-sum(newtree$edge.length[core_edges])-length_correct

  #if you are building a tip-based tree...
  }else if (mode=='tip'){

    #if you are including the root...
    if (rooted==TRUE){

      #make a list of all the taxa present in a threshold number of samples; include the outgroup so that you draw lines back to the root
      coretaxa<-c('outgroup',taxa_names(x)[which(rowSums(sign(otu_table(x)))>=core*length(sample_names(x)))])

      #find the edges associated with those taxa
      core_edges<-which.edge(newtree,coretaxa)

      #sum up the lengths of those edges, correcting for the edge of the added outgroup
      core_length<-sum(newtree$edge.length[core_edges])-0.0001

    #if you are not including the root...
    }else{

      #make a list of all the taxa present in a threshold number of samples; do NOT include the outgroup so that you do not draw lines back to the root
      coretaxa<-taxa_names(x)[which(rowSums(sign(otu_table(x)))>=core*length(sample_names(x)))]

      #find the edges associated with those taxa
      core_edges<-which.edge(phy_tree(x),coretaxa)

      #sum up the lengths of those edges, there is no correction necessary
      core_length<-sum(phy_tree(x)$edge.length[core_edges])

    }

  #if a mode was entered that is not supported, print a warning
  }else{warning('Warning: that mode is not supported')}

  #return the summed length of the edges
  return(core_length)
}
