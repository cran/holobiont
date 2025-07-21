#' @importFrom phyloseq sample_names otu_table phy_tree taxa_names
#' @export
coreJaccard <- function(x,
                        grouping,
                        core_fraction = 0.5,
                        ab_threshold1 = 0,
                        ab_threshold2 = 0,
                        ab_threshold3 = 0,
                        selection = 'basic',
                        max_tax=NULL,
                        increase_cutoff = 2) {


  core<-core_fraction

  #find the number of different habitat types (e.g. hosts or environments) that are being compared
  group_count<-length(unique(grouping))
  group_id<-unique(grouping)

  #make a list of the samples from the first habitat
  group1<-sample_names(x)[which(grouping==group_id[1])]
  #make a list of the samples from the second habitat
  group2<-sample_names(x)[which(grouping==group_id[2])]

  x1<-prune_samples(group1,x)
  x2<-prune_samples(group2,x)


  if (selection == 'basic'){

    #Call function to calculate all edges and core edges
    temp1<-basic_np(x1,core,ab_threshold1,ab_threshold2,ab_threshold3)
    temp2<-basic_np(x2,core,ab_threshold1,ab_threshold2,ab_threshold3)

  }else if (selection == 'shade'){

    #If no maximum taxon is specified, default to considering 1%
    if (is.null(max_tax)){max_tax<-1}

    temp1<-shade_np(x1,max_tax,increase_cutoff)
    temp2<-shade_np(x2,max_tax,increase_cutoff)

  }else{

    stop('That method of core selection is not available. Please use basic or shade.')

  }

  group1core<-temp1$core_taxa
  group2core<-temp2$core_taxa


  Jaccard<-1-length(intersect(group1core,group2core))/length(union(group1core,group2core))
  return(Jaccard)

}
