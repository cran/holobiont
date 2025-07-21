#' @importFrom phytools bind.tip
#' @importFrom phyloseq sample_names otu_table phy_tree taxa_names
#' @importFrom ape which.edge mrca nodepath
#' @importFrom vegan decostand
#' @export
basic_np <-function(xv,cf, abt1, abt2, abt3){

  #Find the otu table from the phyloseq object
  otu<-otu_table(xv)

  #find all of the edges in the tree
  alltaxa<-taxa_names(xv)


  #find all of the taxa present in a threshold number of samples; do not add the outgroup so that edges are drawn on the minimal spanning tree
  if (cf>0){
  coretaxa<-taxa_names(xv)[which(rowSums(sign(otu))>cf*length(sample_names(xv)))]
  }else{
    coretaxa<-taxa_names(xv)[which(rowSums(sign(otu))>=1E-12)]
  }
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

  #output the list of edges and core edges
  output<-list(taxa=alltaxa,core_taxa=coretaxa)

  return(output)

}



