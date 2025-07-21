#' @importFrom castor get_all_distances_to_root
#' @importFrom dplyr left_join
#' @importFrom data.table last
#' @importFrom ape seq_root2tip
#' @export

shade_branch <-function(xv, nt,mxtx,fi,st){

  #This is only available for a rooted tree...

  #Pull out otu tables for abundance and presence/absence
  otu<-otu_table(xv)
  otu_PA<-sign(otu_table(xv))

  #initialize a list of edges that are present in the samples
  edges<-c()

  #for each sample...
  for (i in 1:length(sample_names(xv))){

    #find the taxa in the sample (including the outgroup so that you draw the tree back to the root)
    nz<-c('outgroup',taxa_names(xv)[which(otu[,i]>0)])

    #find the edges associated with those taxa (drawn back to the tree root)
    edges<-c(edges,which.edge(nt,nz))

  }

  #find counts of the number of times each edge appeared across all the samples
  branch_counts<-table(edges)

  #Ensure that the outgroup is the last tip, if it isn't, stop and print a warning
  outno<-which(nt$tip.label=='outgroup')
  tipno<-length(nt$tip.label)
  if (outno!=tipno || length(outno)==0){
    print('WARNING: outgroup is not at the end of the tree or outgroup is not labelled correctly')
    stop('The outgroup is not the final tip in the tree. Please ensure that the outgroup is the final tip in the tree')}

  #Calculate the total nodes = internal + tips in the trree
  node_no<-length(nt$tip.label)+length(nt$node.label)

  #Initialize a matrix where columns = nodes and rows = tree tips
  nodecount<-matrix(0,ncol=node_no,nrow=length(nt$tip.label))

  #Find the nodes in the path between the tip and the root
  toroot<-idk(nt$edge, length(nt$tip.label), nt$Nnode)

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

  #Find the occupancy of each branch that is present
  temper <- branch_counts/ncol(otu)


  #Find mean abundance of each taxon (present or otherwise)
  otu_ab<-apply(decostand(newotu, method="total", MARGIN=2),1, mean)     # mean relative abundance

  #For each node in the tree
  node_ab<-c()
  for (k in 1:node_no){
    #sum up the abundance associated with that node based on abundances of its children tips
    node_ab<-c(node_ab,sum(otu_ab[which(nodecount[,k]==1)]))
  }

  #assign node abundances to their respective edges
  branch_ab<-node_ab[nt$edge[,2]]

  #Initialize a vector of edge occupancies
  branch_occ<-rep(0,length(branch_ab))

  #For each branch (present or otherwise) put the occupancy in the vector based on the branch number
  for (k in 1:length(temper)){
    hitter<-as.integer(rownames(temper)[k])
    branch_occ[hitter] <- temper[k]
  }

  #Make a dataframe of occupancy and abundance
  occ_abun <- data.frame(branch_occ=branch_occ, branch_ab=branch_ab, branch_length=nt$edge.length)

  #Sort first by abundance and then by occupancy
  occ_abun_sorted<-occ_abun[with(occ_abun, order(-branch_occ, -branch_ab)), ]

  #Make a list of the names of the branches (actually numbers) in ranked order
  branch_ranked<-rownames(occ_abun_sorted)

  #Make a matrix with rows as taxon names and columns as nodes (1 = that node is on path to that taxon, 0 = it is not)
  taxa_branches<-as.matrix(nodecount[1:length(rownames(newotu)),])

  #Multiply each sample relative abundance by branches associated with the corresponding taxa
  #Gives a matrix with the relative abundance associated with each node in the phylogenetic tree for each sample
  rbyl<-apply(newotu,2,function(x) colSums(as.vector(t(x))*taxa_branches))

  #Reorder so it is the order of the branches, not nodes
  rbyl<-rbyl[nt$edge[,2],]

  #Reorder so it is in the order of branches with highest to lowest occupancy
  rbyl<-rbyl[as.integer(branch_ranked),]
  rownames(rbyl)<-branch_ranked

  #Initialize list of weighted UniFrac contributions
  BCaddition<-NULL

  #Begin with all the branches present on every sample
  if (is.null(st)){
    start_tree<-length(which(occ_abun_sorted$branch_occ==1))
  }else{start_tree<-st}
  start_tree<-max(2,start_tree)


  #Consider only the branches present on every sample
  rbyl_pruned<-rbyl[1:start_tree,]

  #Make a vector of the these branch lengths
  bl_pruned<-occ_abun_sorted[1:start_tree,]$branch_length

  #Make a vector of all branch lengths
  bl_full<-occ_abun_sorted$branch_length

  #Calculate the contribution to weighted Unifrac distance from the first two branches
  bc<-apply(combn(length(colnames(rbyl)), 2),2,function(y) sum(as.vector(abs(rbyl_pruned[,y[1]]-rbyl_pruned[,y[2]]))*as.vector(bl_pruned))/sum(as.vector(abs(rbyl[,y[1]]+rbyl[,y[2]]))*as.vector(bl_full)))
  bc_names<-apply(combn(length(colnames(rbyl)), 2), 2, function(y) paste(colnames(rbyl)[y], collapse=' - '))

  #Store the Bray-Curtis contributions in the list
  df_s <- data.frame(bc_names,bc)
  names(df_s)[2] <- length(which(occ_abun_sorted$branch_occ==1))
  BCaddition <- rbind(BCaddition,df_s)

  max_branches<-round(mxtx*length(branch_ranked)/100)
  max_branches<-max(max_branches,start_tree+10)

  #For the remaining taxa up to the highest taxon you want to consider...
  for(i in 2:max_branches){

    new_tree<-start_tree+i

    #Consider only the next branches
    rbyl_pruned<-rbyl[1:new_tree,]
    #Make a vector of branch lengths
    bl_pruned<-occ_abun_sorted[1:new_tree,]$branch_length

    #Calculate the contribution to weighted Unifrac distance from the first two branches
    bc<-apply(combn(length(colnames(rbyl)), 2),2,function(y) sum(as.vector(abs(rbyl_pruned[,y[1]]-rbyl_pruned[,y[2]]))*as.vector(bl_pruned))/sum(as.vector(abs(rbyl[,y[1]]+rbyl[,y[2]]))*as.vector(bl_full)))
    bc_names<-apply(combn(length(colnames(rbyl)), 2), 2, function(y) paste(colnames(rbyl)[y], collapse=' - '))

    #Store the UniFrac contributions in the list
    df_s <- data.frame(bc_names,bc)
    names(df_s)[2] <- new_tree
    BCaddition <- left_join(BCaddition,df_s, by='bc_names')

  }

  #If you have not considered the weighted UniFrac contributions including all taxa, do a final calculation for the full microbiota
  if (max_branches<length(nt$edge[,2])){

    #Find the total Bray-Curtis distances
    bc<-apply(combn(length(colnames(rbyl)), 2),2,function(y) sum(as.vector(abs(rbyl[,y[1]]-rbyl[,y[2]]))*as.vector(bl_full))/sum(as.vector(abs(rbyl[,y[1]]+rbyl[,y[2]]))*as.vector(bl_full)))
    bc_names<-apply(combn(length(colnames(rbyl)), 2), 2, function(y) paste(colnames(rbyl)[y], collapse=' - '))

    #Store the total Bray-Curtis distances in the list and rename it as BCfull
    df_full <- data.frame(bc_names,bc)
    names(df_full)[2] <- length(nt$edge[,2])
    BCfull <- left_join(BCaddition,df_full, by='bc_names')

  }else{

    #If you have already considered all the branches...
    BCfull<-BCaddition
  }


  #Name the rows of the BCfull list based on the pairs of samples each row considers
  rownames(BCfull) <- BCfull$bc_names
  #Store the BCfull matrix in second, temporary location
  temp_BC <- BCfull
  #Remove the bc_names column from the temporary location so you can turn it into a matrix
  temp_BC$bc_names <- NULL
  #Turn it into a matrix
  temp_BC_matrix <- as.matrix(temp_BC)

  #Find the mean weighted UniFrac contributions for the first taxon, the first two taxa, the first three taxa... etc.
  MeanBC<-colMeans(temp_BC_matrix,na.rm=TRUE)
  #Express the mean weighted UniFrac contributions as a fraction of the total weighted UniFrac index for each pair of samples
  proportionBC<-MeanBC/max(MeanBC)
  #Find the increase in weighted UniFrac contribution (factor) for each added taxon
  Increase=c(0,MeanBC[-1]/MeanBC[-length(MeanBC)])

  #Make a dataframe of mean BC, proportion BC and Increase in BC contributions
  Rank<-1:1:length(Increase)
  BC_ranked<-data.frame(Rank,MeanBC,proportionBC,Increase)
  #Remove the final row (the full microbiota)
  BC_ranked <- BC_ranked[-nrow(BC_ranked),]

  #Find the taxon
  criterion<-1+fi/100
  lastCall <- last(as.numeric(as.character(BC_ranked$Rank[(BC_ranked$Increase>=criterion)])))

  if (length(which(BC_ranked$Increase>=criterion))==0){
    lastCall<-0
    warning('Including only initial branches, consider lowering initial branches')
  }
  #find the edges associated with all of the core taxa identified above
  endme<-length(which(occ_abun_sorted$branch_occ==1))+lastCall
  #print(occ_abun_sorted)
  core_edges<-as.integer(branch_ranked[1:endme])

  output<-list(edges=edges,core_edges=core_edges)

  return(output)

}

idk <- function (edge, nbtip, nbnode)
  (.Call)(seq_root2tip, edge, nbtip, nbnode)

