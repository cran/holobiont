#' @importFrom phyloseq sample_names otu_table phy_tree taxa_names
#' @export
coreJaccard <- function(x, grouping, core_fraction) {

  #find the number of different habitat types (e.g. hosts or environments) that are being compared
  group_count<-length(unique(grouping))
  group_id<-unique(grouping)

  #make a list of the samples from the first habitat
  group1<-sample_names(x)[which(grouping==group_id[1])]
  #make a list of the samples from the second habitat
  group2<-sample_names(x)[which(grouping==group_id[2])]

  if (core_fraction > 0){
  group1core<-taxa_names(x)[which(rowSums(sign(otu_table(x)[,which(grouping==group_id[1])]))>core_fraction*length(group1))]
  group2core<-taxa_names(x)[which(rowSums(sign(otu_table(x)[,which(grouping==group_id[2])]))>core_fraction*length(group2))]
  }else{
    group1core<-taxa_names(x)[which(rowSums(sign(otu_table(x)[,which(grouping==group_id[1])]))>1E-12)]
    group2core<-taxa_names(x)[which(rowSums(sign(otu_table(x)[,which(grouping==group_id[2])]))>1E-12)]
  }

  Jaccard<-1-length(intersect(group1core,group2core))/length(union(group1core,group2core))
  return(Jaccard)
}
