#' @importFrom ggplot2 coord_cartesian
#' @importFrom phytools bind.tip
#' @importFrom phyloseq sample_names otu_table phy_tree taxa_names prune_samples
#' @importFrom ape which.edge mrca nodepath
#' @importFrom utils combn
#' @export
#'
corePhyloVenn <- function(x,
                          grouping,
                          core_fraction = 0.5,
                          ab_threshold1 = 0,
                          ab_threshold2 = 0,
                          ab_threshold3 = 0,
                          mode = 'branch',
                          selection = 'basic',
                          max_tax=NULL,
                          increase_cutoff = 2,
                          initial_branches=NULL,
                          rooted=TRUE,
                          ordered_groups=NULL,
                          show_percentage=TRUE,
                          decimal=2,
                          fill_color=c('red','orange','yellow','green','blue','purple','black'),
                          fill_alpha=0.5,
                          stroke_color='black',
                          stroke_alpha = 1,
                          stroke_size = 1,
                          stroke_linetype = "solid",
                          set_name_color = "black",
                          set_name_size = 6,
                          text_color = "black",
                          text_size = 4) {

  core<-core_fraction

  #Using an unrooted tree in branch mode with abundance thresholds is not currently supported.  If an unrooted tree is provided with non-zero abundance thresholds and branch mode is requested, reroot it
  if (rooted == FALSE){
    if (ab_threshold1+ab_threshold2+ab_threshold3>0){
      randoroot = sample(phy_tree(x)$tip.label, 1)
      warning("Use of an unrooted tree with non-zero abundance thresholds in branch mode is not supported.")
      warning("Randomly assigning root as -- ", randoroot, " -- in the phylogenetic tree provided.")
      phy_tree(x) <- root(phy=phy_tree(x), outgroup=randoroot, resolve.root=TRUE, interactive=FALSE)
      if( !is.rooted(phy_tree(x)) ){
        stop("Problem automatically rooting tree. Make sure your tree is rooted before attempting UniFrac calculation. See ?ape::root")
      }
      rooted = TRUE
    }
  }

  #find the number of different habitat types (e.g. hosts or environments) that are being compared
  group_count<-length(unique(grouping))

  #if no group order is specified, pick arbitrarily
  if (is.null(ordered_groups)){
    group_id<-unique(grouping)
    #otherwise use the specified group order
  }  else{group_id<-ordered_groups}

  #Venn diagrams can only be drawn for 7 or less habitats
  #If more than 7 habitats are entered, print out a warning
  if (group_count>7){
    warning('Warning: Too many habitat types!')
  }

  #otherwise proceed...
  else{

    #add an outgroup so that there is the option of drawing all edges to the root rather than the minimal spanning tree
    newtree<-bind.tip(phy_tree(x),tip.label='outgroup',edge.length=0.0001,position=0)

    #if you are building a branch-based tree...
    if (mode=='branch'){

      #initialize a list of lists where the main list is of the habitats and the sublists are the samples associated with each habitat
      grouplist<-list()

      #initialize a list of lists where the main list is the habitats and the sublists are the edges associated with core microbes from each habitat
      edgelist<-list()

      #for each habitat...
      for (i in 1:group_count){
        #find the samples from that habitat
        temp<-list(sample_names(x)[which(grouping==group_id[i])])

        #put them into the list of samples in each habitat
        grouplist<-append(grouplist,temp)

        #Make a phyloseq object containing only samples from the focal habitat
        xg<-prune_samples(sample_names(x)[which(grouping==group_id[i])],x)

        #If using the basic selection mode...
        if (selection == 'basic'){

          #Call function to calculate all edges and core edges
          temp2<-basic_branch(xg,newtree,core,ab_threshold1,ab_threshold2,ab_threshold3,rooted)
        #If using the shade selection mode...
        }else if (selection == 'shade'){

          #If no maximum taxon is specified, default to considering 1%
          if (is.null(max_tax)){max_tax<-1}

          #Call the function to calculate all edges and core edges
          temp2<-shade_branch(xg,newtree,max_tax,increase_cutoff,initial_branches)

        #If a non-supported mode is selected print a warning and stop
        }else{

          stop('That method of core selection is not available. Please use basic or shade.')

        }

        #make a list of the core edges from the focal habitat
        core_edges<-temp2$core_edges

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

        #append the edges from the focal habitat to the list of lists where the main list is of the habitats and the sublists are of the core edges for each habitat
        edgelist<-append(edgelist,list(core_edges))

      }
    }

    #if you are building a tip-based tree...
    else if (mode=='tip'){


      #initialize a list of lists where the main list is of the habitats and the sublists are the samples associated with each habitat
      grouplist<-list()

      #initialize a list of lists where the main list is the habitats and the sublists are the names of the core taxa associated with each habitat
      corelist<-list()

      #initialize a vector of all the names of the core taxa across any/all habitats
      allcorelist<-c()

      #for each habitat...
      for (i in 1:group_count){

        #find the samples from that habitat
        temp<-list(sample_names(x)[which(grouping==group_id[i])])

        #put them into the list of samples in the habitat list
        grouplist<-append(grouplist,temp)

        #make a phyloseq object for habitat
        xg<-prune_samples(sample_names(x)[which(grouping==group_id[i])],x)

        if (selection == 'basic'){

          #Call function to calculate all edges and core edges
          temp2<-basic_tip(xg,newtree,core,ab_threshold1,ab_threshold2,ab_threshold3,rooted)

        }else if (selection == 'shade'){

          #If no maximum taxon is specified, default to considering 1%
          if (is.null(max_tax)){max_tax<-1}
          temp2<-shade_tip(xg,newtree,max_tax,increase_cutoff)

        }else{

          warning('That method of core selection is not available. Please use basic or shade.')

        }

        #put them into the list of taxa in the habitat list
        corelist<-append(corelist,list(temp2$core_taxa))

        #put them into the vector of core taxa from any/all habitats
        allcorelist<-unique(c(temp$core_taxa,allcorelist))

      }

      #find the edges associated with core taxa from any/all habitats (this is the minimal spanning tree)
      spanlist<-which.edge(newtree,allcorelist)
      habitatspanlist<-c()
      for (i in 1:group_count){
        habitatspanlist<-c(habitatspanlist,which.edge(newtree,corelist[[i]]))
      }

      #make a list of lists, with the main list being habitats and the sublists being the edges associated with core taxa in each habitat
      edgelist<-list()
      for (i in 1:group_count){
        #initially, include the edges associated with the root
        edgelist<-append(edgelist,list(which.edge(newtree,c('outgroup',corelist[[i]]))))
      }


      #if you are including the root...
      if (rooted ==TRUE){
        #if you are not including the root...
      }else{
        #for each habitat...
        for (i in 1:group_count){
          #remove all edges that are not part of the minimal spanning tree for core microbes from each habitat
          edgelist[[i]]<-edgelist[[i]][which(edgelist[[i]] %in% habitatspanlist[[i]])]
        }
      }

      #if a mode that is not supported is entered, print a warning
    }else{stop('Warning: that mode is not supported')}

    #initialize a list of possible habitat combinations
    combos<-list()

    #initialize a list of lists, with the main lists being habitat combinations and the sublists being the branch lengths shared by the habitat combinations
    intersections<-list()

    #initialize a list of total branch length shared by each habitat combination
    lengths<-c()

    #for anywhere from 1 to n combinations (i.e., {1,2} is a length 2 combination, {1,4,5} is a length 3 combination), where n is the total number of habitats...
    for (i in 1:group_count){


      #list all of the possible combinations of that size (e.g., {1,2},{1,3},{2,3})
      ff<-combn(group_count,i)

      #count the possible combinations of that size
      no_combinations<-length(ff[1,])

      #find the length of the combination
      combination_length<-length(ff[,1])

      #for each possible combination of the focal size...
      for (j in 1:no_combinations){

        #add the particular combination of habitats to your list of possible habitat combinations
        combos<-append(combos,list(ff[,j]))

        #find the habitats that are outside the particular combination
        outside<-which(!(1:1:group_count %in% ff[,j]))

        #find the edges in the first habitat of the combination
        shared_temp<-edgelist[[ff[1,j]]]


        #if there are more habitats in the combination...
        if (combination_length>1){

          #find the edges that are shared by all habitats in the combination
          for (k in 2:combination_length){
            shared_temp<-intersect(shared_temp,edgelist[[ff[k,j]]])
          }

        }

        #if there are any habitats outside the combination...
        if (length(outside)>0){

          #for each habitat outside...
          for (g in 1:length(outside)){
            #remove the edges that it shares with the focal combination of habitats (shared edges must be exclusive to the focal habitat combination)
            shared_temp<-setdiff(shared_temp,edgelist[[outside[g]]])
          }

        }

        #add the edges that are exclusively shared with the focal combination of habitats to the list of shared edges for each habitat combination
        intersections<-append(intersections,list(shared_temp))


        #add the lengths of the branches shared with the focal combination of habitats to the list of lengths for each habitat combination
        if (rooted==TRUE){

          #if the root was included, then you need to remove the branch length for the added outgroup that was shared with everyone
          if (length(shared_temp)>0){
          lengths<-c(lengths,sum(newtree$edge.length[shared_temp])-0.0001)
          }else{lengths<-c(lengths,0)}

        }else{

          #if the root was not included, then there is no correction to branch lengths
          if (length(shared_temp)>0){
            lengths<-c(lengths,sum(newtree$edge.length[shared_temp]))
          }else{lengths<-c(lengths,0)}

        }
      }
    }

    if (sum(lengths)==0){stop('There are no core taxa in any group')}

    #if the Venn diagram is shown as a percentage of the branch lengths...
    if (show_percentage){
      #round to the decimal desired, multiple by 100 and divide by the total length
      fractions<-round(10^(decimal+2)*lengths/sum(lengths))
    }else{
      #otherwise just round to the decimal desired
      fractions<-round(10^(decimal)*lengths)
    }
    #note that the above are in whole numbers... so if you want 25.56% you have 2556 in that slot...
    #this is for making the dummy dataframe that can be read by the Venn diagram package
    #find the cumulative percentages or branch lengths across the list of possible habitat combinations
    cumfractions<-cumsum(fractions)

    #initialize a dataframe for the Venn diagram package with FALSE for every entry
    #This dataframe has as many rows as the total value of the cumulative percentages or branch lengths in your cumfractions vector
    #So, for instance, if you want percentages with 0 decimals, you'd have 100 rows of FALSE entries... if you want percentages with 2 decimals you'd have 10000 rows of FALSE entries
    final<-data.frame(rep(FALSE,cumfractions[length(cumfractions)]))
    #There should be as many columns as there are habitats
    for (i in 2:group_count){
      final<-cbind(final,rep(FALSE,cumfractions[length(cumfractions)]))
    }
    #for each row in the dataframe
    counter<-1
    for (i in 1:length(cumfractions)){
      if (counter<cumfractions[i]){
      #find the habitat combinations associated with those percentages and set them to TRUE
      final[counter:cumfractions[i],combos[i][[1]]]<-TRUE
      counter<-cumfractions[i]+1
      }
    }

    #name the dataframe columns based on the habitats
    colnames(final)<-group_id
    #Plot the Venn diagram
    ggvenn2(final,columns=NULL,show_elements=FALSE,show_percentage,digits=decimal,fill_color,fill_alpha,stroke_color,stroke_alpha,stroke_size,stroke_linetype,set_name_color,set_name_size,text_color,text_size)+coord_cartesian(clip="off")
  }
}
