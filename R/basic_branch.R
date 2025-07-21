#' @importFrom castor get_all_distances_to_root
#' @importFrom ape seq_root2tip
#' @export

basic_branch <-function(xv, nt, cf, abt1, abt2, abt3,rt){


  #initialize a list of edges that are present in the samples
  edges<-c()

  #for each sample...
  for (i in 1:length(sample_names(xv))){

    if (rt == TRUE){

      #find the taxa in the sample (including the outgroup so that you draw the tree back to the root)
      nz<-c('outgroup',taxa_names(xv)[which(otu_table(xv)[,i]>0)])

    }
    else{

      #find the taxa in the sample (not including an outgroup so that you don't draw the tree back to the root)
      nz<-taxa_names(xv)[which(otu_table(xv)[,i]>0)]

    }

    #find the edges associated with those taxa
    edges<-c(edges,which.edge(nt,nz))
  }

  #find counts of the number of times each edge appeared across all the samples
  branch_counts<-table(edges)

  #pull out the edges that were present in at least a core threshold number of samples
  core_branch<-which(branch_counts>cf*length(sample_names(xv)))

  #make a list of the core edges based on incidence threshold
  core_edges<-as.integer(names(core_branch))

  #If there are non-zero abundance thresholds
  if (abt1+abt2+abt3>0){

    #Ensure that the outgroup is the last tip, if it isn't, stop and print a warning
    outno<-which(nt$tip.label=='outgroup')
    tipno<-length(nt$tip.label)
    if (outno!=tipno){
      warning('WARNING: outgroup is not at the end of the tree')
      stop('The outgroup is not the final tip in the tree. Please ensure that the outgroup is the final tip in the tree')}

    #Calculate the total nodes = internal + tips in the tree
    node_no<-length(nt$tip.label)+length(nt$node.label)
    #Initialize a matrix where columns = nodes and rows = tree tips
    nodecount<-matrix(0,ncol=node_no,nrow=length(nt$tip.label))
    #Find the nodes in the path between the tip and the root
    toroot<-idk2(nt$edge, length(nt$tip.label), nt$Nnode)
    #Each node x tip element equals 0 if the node is not part of the path from the tip to the root and 1 if it is
    for (k in 1:tipno){
      nodecount[k,toroot[[k]]]<-1
    }

    #Make sure that the otu table taxa are in the same order as the tree taxa
    match_order<-match(rownames(otu_table(xv)),nt$tip.label)
    oldotu<-otu_table(xv)
    newotu<-oldotu[order(match_order),]

    #Normalize the otu table
    newotu<-decostand(newotu,method='total',2)

    #Add the outgroup as the final row
    outgroup<-rep(0,length(newotu[2,]))
    newotu<-rbind(newotu,outgroup)

    #Initialize the abundance vector (vector of abundances associated with each edge)
    node_abs1<-c()  #mean abundance across samples for each node
    node_abs2<-c()  #maximum abundance across samples for each node
    node_abs3<-c()  #minimum abundance across samples for each node

    #For each node in the tree
    for (k in 1:node_no){

      #If it leads to more than 1 tip

      if (length(which(nodecount[,k]==1))>1){

        #Sum up all of the relative abundance from each tip contributing to the node
        node_abs<-colSums(newotu[which(nodecount[,k]==1),])
        #Find the mean value across samples
        node_abs1<-c(node_abs1,mean(node_abs))
        #Find the maximum value across samples
        node_abs2<-c(node_abs2,max(node_abs))
        #Find the minimum value across samples
        node_abs3<-c(node_abs3,min(node_abs))

      }else{

        #Find the mean value across samples
        node_abs1<-c(node_abs1,mean(newotu[which(nodecount[,k]==1),]))
        #Find the mean value across samples
        node_abs2<-c(node_abs2,max(newotu[which(nodecount[,k]==1),]))
        #Find the mean value across samples
        node_abs3<-c(node_abs3,min(newotu[which(nodecount[,k]==1),]))

      }
    }

    #Reorder the node abundance list to reflect the associated edge by using the number of the associated child node for each edge
    #Multiply the abundances by 100 to reflect percentages
    edge_abs1<-100*node_abs1[nt$edge[,2]]
    edge_abs2<-100*node_abs2[nt$edge[,2]]
    edge_abs3<-100*node_abs3[nt$edge[,2]]

    #Find the edges above the abundance thresholds
    abs1<-which(edge_abs1>=abt1)
    abs2<-which(edge_abs2>=abt2)
    abs3<-which(edge_abs3>=abt3)

    #Ensure that the core edges not only meet the occupancy criterion, but also the abundance criterion
    core_edges<-intersect(intersect(intersect(core_edges,abs1),abs2),abs3)


  }

  #if the root is not included (only implemented for abundance thresholds of zero!)
  if (rt==FALSE){

    #find the nodes associated with core edges
    nodes<-unique(c(nt$edge[,1][core_edges],nt$edge[,2][core_edges]))

    #find the mrca of each node in the tree
    cc<-mrca(nt,full=TRUE)

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
          mrca_list<-c(mrca_list,nodepath(nt,from=missing[i],to=nodes[j]))
        }
      }
    }

    #find the edges associated with all the nodes (core and mrcas)
    all_core_edges<-intersect(which(nt$edge[,1] %in% mrca_list),which(nt$edge[,2] %in% mrca_list))

    #add the missing edges to the core edges
    core_edges<-unique(c(core_edges,all_core_edges))

  }

  output<-list(edges=edges,core_edges=core_edges)

  return(output)

}


idk2 <- function (edge, nbtip, nbnode)
  (.Call)(seq_root2tip, edge, nbtip, nbnode)

