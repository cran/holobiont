#' @export
shade_np <-function(xv,mxtx,fi){

  #find all of the edges in the tree
  alltaxa<-taxa_names(xv)

  #Pull out otu tables from the phyloseq object
  otu<-otu_table(xv)
  #Turn it into a presence/absence table as well
  otu_PA<-sign(otu_table(xv))

  #Find occupancy of each taxon
  otu_occ <- rowSums(otu_PA)/ncol(otu_PA)

  #Find mean relative abundance of each taxon
  otu_ab<-apply(decostand(otu, method="total", MARGIN=2),1, mean)     # mean relative abundance

  #Make a dataframe of occupancy and abundance
  occ_abun <- data.frame(otu_occ=otu_occ, otu_ab=otu_ab)

  #Sort first by abundance and then by occupancy
  occ_abun_sorted<-occ_abun[with(occ_abun, order(-otu_occ, -otu_ab)), ]

  #Find the ranked list
  otu_ranked<-rownames(occ_abun_sorted)

  #Initialize list of Bray-Curtis contributions
  BCaddition<-NULL

  #Pick the highest ranked taxon
  otu_start=which(taxa_names(otu)==otu_ranked[1])

  #Make an otu table matrix
  start_matrix <- as.matrix(otu[otu_start,])

  #Find the contribution to Bray-Curtis distance made by that taxon
  bc <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(sum(colSums(otu)[c(x[1],x[2])])))
  bc_names<-apply(combn(ncol(start_matrix), 2), 2, function(y) paste(colnames(start_matrix)[y], collapse=' - '))

  #Store the Bray-Curtis contributions in the list
  df_s <- data.frame(bc_names,bc)
  names(df_s)[2] <- 1
  BCaddition <- rbind(BCaddition,df_s)

  max_tips<-round(mxtx*length(otu[,1])/100)
  max_tips<-max(10,max_tips)


  #For the remaining taxa up to the highest taxon you want to consider...
  for(i in 2:max_tips){

    #Pick the next highest ranked taxon
    otu_add<-which(taxa_names(otu)==otu_ranked[i])

    #Add it to the otu table
    add_matrix <- as.matrix(otu[otu_add,])
    start_matrix <- rbind(start_matrix, add_matrix)

    #Find the contributions to Bray-Curtis distance made by the taxa up to the current higest ranked taxa
    bc <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(sum(colSums(otu)[c(x[1],x[2])])))
    bc_names<-apply(combn(ncol(start_matrix), 2), 2, function(y) paste(colnames(start_matrix)[y], collapse=' - '))

    #Store the Bray-Curtis contributions in the list
    df_s <- data.frame(bc_names,bc)
    names(df_s)[2] <- i
    BCaddition <- left_join(BCaddition,df_s, by='bc_names')

  }

  #If you have not considered the Bray-Curtis contributions including all taxa, do a final calculation for the full microbiota
  if (max_tips<length(taxa_names(otu))){

    #Find the total Bray-Curtis distances
    bc <- apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]- otu[,x[2]]))/(sum(colSums(otu)[c(x[1],x[2])])))
    bc_names<-apply(combn(ncol(otu), 2), 2, function(y) paste(colnames(otu)[y], collapse=' - '))

    #Store the total Bray-Curtis distances in the list and rename it as BCfull
    df_full <- data.frame(bc_names,bc)
    names(df_full)[2] <- length(rownames(otu))
    BCfull <- left_join(BCaddition,df_full, by='bc_names')

  }else{
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

  #Find the mean Bray-Curtis contributions for the first taxon, the first two taxa, the first three taxa... etc.
  MeanBC<-colMeans(temp_BC_matrix)
  #Express the mean Bray-Curtis contributions as a fraction of the total Bray-Curtis index for each pair of samples
  proportionBC<-MeanBC/max(MeanBC)
  #Find the increase in Bray-Curtis contribution (factor) for each added taxon
  Increase=c(0,MeanBC[-1]/MeanBC[-length(MeanBC)])

  #Make a dataframe of mean BC, proportion BC and Increase in BC contributions
  Rank<-1:1:length(Increase)
  BC_ranked<-data.frame(Rank,MeanBC,proportionBC,Increase)
  #Remove the final row (the full microbiota)
  BC_ranked <- BC_ranked[-nrow(BC_ranked),]

  #Find the taxon
  criterion<-1+fi/100
  lastCall <- last(as.numeric(as.character(BC_ranked$Rank[(BC_ranked$Increase>=criterion)])))

  coretaxa<-otu_ranked[1:lastCall]


  output<-list(core_taxa=coretaxa,taxa=alltaxa)

  return(output)

}




