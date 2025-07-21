#' @export
basic_tip <-function(xv, nt,cf, abt1, abt2, abt3,rt){

  #Find the otu table from the phyloseq object
  otu<-otu_table(xv)

  #find all of the edges in the tree
  alltaxa<-taxa_names(xv)
  edges<-which.edge(nt,alltaxa)


  #if you want to include the root in the tip-based tree
  if (rt==TRUE){

    #find all of the taxa present in a threshold number of samples; add the outgroup so that edges are drawn back to the root
    coretaxa<-c('outgroup',taxa_names(xv)[which(rowSums(sign(otu))>cf*length(sample_names(xv)))])

    #normalize the otu table
    notu<-decostand(otu,method='total',2)

    #Find all of the taxa present at a threshold mean abundance
    coreab1<-c('outgroup',taxa_names(xv)[which(rowMeans(notu)>=abt1/100)])

    #Find all of the taxa present at a threshold max abundance
    coreab2<-c('outgroup',taxa_names(xv)[which(apply(notu, 1, FUN = max)>=abt2/100)])

    #Find all of the taxa present at a threshold min abundance
    coreab3<-c('outgroup',taxa_names(xv)[which(apply(notu, 1, FUN = min)>=abt3/100)])

    #Find the intersection of all the taxa
    coretaxa<-intersect(intersect(intersect(coretaxa,coreab1),coreab2),coreab3)

  }else{

    #find all of the taxa present in a threshold number of samples; do not add the outgroup so that edges are drawn on the minimal spanning tree
    coretaxa<-taxa_names(xv)[which(rowSums(sign(otu))>cf*length(sample_names(xv)))]

    #normalize the otu table
    notu<-decostand(otu,method='total',2)

    #Find all of the taxa present at a threshold mean abundance
    coreab1<-c(taxa_names(x)[which(rowMeans(notu)>=abt1/100)])

    #Find all of the taxa present at a threshold max abundance
    coreab2<-c(taxa_names(x)[which(apply(notu, 1, FUN = max)>=abt2/100)])

    #Find all of the taxa present at a threshold min abundance
    coreab3<-c(taxa_names(x)[which(apply(notu, 1, FUN = min)>=abt3/100)])

    #Find the intersection of all the taxa
    coretaxa<-intersect(intersect(intersect(coretaxa,coreab1),coreab2),coreab3)

  }

  #find the edges associated with all of the core taxa identified above
  core_edges<-which.edge(nt,coretaxa)

  #output the list of edges and core edges
  output<-list(taxa=alltaxa,edges=edges,core_taxa=coretaxa,core_edges=core_edges)

  return(output)

}



