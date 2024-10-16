#' @importFrom phytools bind.tip
#' @importFrom phyloseq sample_names otu_table phy_tree taxa_names
#' @importFrom ape which.edge mrca nodepath
#' @export
coreUniFrac <- function(x, grouping, core_fraction, mode = 'branch',rooted=TRUE) {

  core<-core_fraction
  #find the number of different habitat types (e.g. hosts or environments) that are being compared
  group_count<-length(unique(grouping))
  group_id<-unique(grouping)

  #make a list of the samples from the first habitat
  group1<-sample_names(x)[which(grouping==group_id[1])]
  #make a list of the samples from the second habitat
  group2<-sample_names(x)[which(grouping==group_id[2])]

  #if you are building a branch-based tree...
  if (mode=='branch'){

    #add an outgroup so that there is the option of drawing all edges to the root rather than the minimal spanning tree
    newtree<-bind.tip(phy_tree(x),tip.label='outgroup',edge.length=0.0001,position=0)

    #initialize lists of edges that are present in each of the two habitats
    edges1<-c()
    edges2<-c()

    #for each sample...
    for (i in 1:length(sample_names(x))){

      #if you are including the root...
      if (rooted == TRUE){

        #find the taxa present in the sample; include the outgroup so that lines are drawn back to the root; account for the case where core = 0 by including taxa that aren't present in any samples
        if (core>0){
          nz<-c('outgroup',taxa_names(x)[which(otu_table(x)[,i]>0)])
        }else{
          nz<-c('outgroup',taxa_names(x)[which(otu_table(x)[,i]>=0)])
        }
        lengther<-0.0001

        #if you aren't including the root...
      }else{

        #find the taxa present in the sample; do NOT include the outgroup so that lines are not drawn back to the root; account for the case where core = 0 by including taxa that aren't present in any samples
        if (core>0){
          nz<-taxa_names(x)[which(otu_table(x)[,i]>0)]
        }else{
          nz<-taxa_names(x)[which(otu_table(x)[,i]>=0)]
        }
        lengther<-0

      }

      #find the edges associated with those taxa
      if (sample_names(x)[i] %in% group1){
        edges1<-c(edges1,which.edge(newtree,nz))
      }else{
        edges2<-c(edges2,which.edge(newtree,nz))
      }

    }

    #find counts of the number of times each edge appeared across each habitat
    branch_counts1<-table(edges1)
    branch_counts2<-table(edges2)

    #pull out the core edges that were present in at least a core number of samples from each habitat
    core_branch1<-which(branch_counts1>=core*length(group1))
    core_branch2<-which(branch_counts2>=core*length(group2))

    #make a list of the core edges for each habitat
    core_edges1<-as.integer(names(core_branch1))
    core_edges2<-as.integer(names(core_branch2))

    #if the root is not included
    if (rooted==FALSE){

      #find the nodes associated with core edges
      nodes1<-unique(c(newtree$edge[,1][core_edges1],newtree$edge[,2][core_edges1]))
      nodes2<-unique(c(newtree$edge[,1][core_edges2],newtree$edge[,2][core_edges2]))
      #find the mrca of each node in the tree
      cc<-mrca(newtree,full=TRUE)
      #find the mrca of each node associated with a core edge
      mrca_matrix1<-cc[nodes1,nodes1]
      mrca_matrix2<-cc[nodes2,nodes2]
      #find the unique mrcas for the core edge nodes
      mrca_list1<-unique(as.vector(mrca_matrix1))
      mrca_list2<-unique(as.vector(mrca_matrix2))
      #find the unique mrcas plus core edge nodes
      mrca_list1<-unique(mrca_list1,nodes1)
      mrca_list2<-unique(mrca_list2,nodes2)
      #identify mrcas missing from the list of nodes associated with core edges
      missing1<-mrca_list1[which(!(mrca_list1 %in% nodes1))]
      missing2<-mrca_list2[which(!(mrca_list2 %in% nodes2))]
      if (length(missing1)>0){
      for (i in 1:length(missing1)){
        for (j in 1:length(nodes1)){
          #find the nodes connecting the missing mrcas to the nodes associated with core edges
          mrca_list1<-c(mrca_list1,nodepath(newtree,from=missing1[i],to=nodes1[j]))
        }
      }
      }
      if (length(missing2)>0){
      for (i in 1:length(missing2)){
        for (j in 1:length(nodes2)){
          #find the nodes connecting the missing mrcas to the nodes associated with core edges
          mrca_list2<-c(mrca_list2,nodepath(newtree,from=missing2[i],to=nodes2[j]))
        }
      }
      }
      #find the edges associated with all the nodes (core and mrcas)
      all_core_edges1<-intersect(which(newtree$edge[,1] %in% mrca_list1),which(newtree$edge[,2] %in% mrca_list1))
      all_core_edges2<-intersect(which(newtree$edge[,1] %in% mrca_list2),which(newtree$edge[,2] %in% mrca_list2))

      #add the missing edges to the core edges
      core_edges1<-unique(c(core_edges1,all_core_edges1))
      core_edges2<-unique(c(core_edges2,all_core_edges2))
    }

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

    unifrac<-(A+B)/(sum(newtree$edge.length[unique(c(core_edges1,core_edges2))])-lengther)

  #if you are building a tip-based tree...
  }else if (mode == 'tip'){

    #find the taxa that are core to the first habitat
    coretaxa1<-taxa_names(x)[which(rowSums(sign(otu_table(x)[,group1]))>=core*length(group1))]
    #find the taxa that are core to the second habitat
    coretaxa2<-taxa_names(x)[which(rowSums(sign(otu_table(x)[,group2]))>=core*length(group2))]

    #add an outgroup so that there is the option of drawing all edges to the root rather than the minimal spanning tree
    newtree<-bind.tip(phy_tree(x),tip.label='outgroup',edge.length=0.0001,position=0)
    lengther<-0.0001

    #find the core edges, including any root edges, for the first habitat
    core_edges1<-which.edge(newtree,c('outgroup',coretaxa1))
    #find the core edges, including any root edges, for the second habitat
    core_edges2<-which.edge(newtree,c('outgroup',coretaxa2))
    #find the minimal spanning tree for both habitats (this does not include the root edges)
    spanlist<-which.edge(phy_tree(x),unique(c(coretaxa1,coretaxa2)))
    #find the minimal spanning trees for each habitat independently
    habitatspanlist1<-which.edge(phy_tree(x),coretaxa1)
    habitatspanlist2<-which.edge(phy_tree(x),coretaxa2)

    #if you do not want to include the root, remove the edges not associated with the particular habitat (these are not in the minimal spanning tree for that habitat)
    if (rooted==FALSE){
      core_edges1<-core_edges1[which(core_edges1 %in% habitatspanlist1)]
      core_edges2<-core_edges2[which(core_edges2 %in% habitatspanlist2)]
      lengther<-0
    }

    #find the edges that are core in the first habitat but not the second
    in1only<-which(!(core_edges1 %in% core_edges2))
    #find the edges that are core in the second habitat but not the first
    in2only<-which(!(core_edges2 %in% core_edges1))
    #divide the length of the unique core edges by the total length of the core edges
    if (length(in1only)>0){
      A<-sum(newtree$edge.length[core_edges1[in1only]])
    }else{A<-0}
    if (length(in2only)>0){
      B<-sum(newtree$edge.length[core_edges2[in2only]])
    }else{B<-0}

    unifrac<-(A+B)/(sum(newtree$edge.length[unique(c(core_edges1,core_edges2))])-lengther)
    }

  #return the unifrac distance
  return(unifrac)

}
