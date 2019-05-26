
getAncestors<-phytools:::getAncestors


plotBestPathTree<-function(tree,best_path_counts,branch_counts_df,path_scores_df, position_in_branch, best_node,...){
	plot(tree, cex=0.15, edge.col=ifelse(branch_counts_df$Edge %in% unique(path_scores_df$stopped_edges), yes="red", no='lightgrey'), show.tip.label = T, edge.width = 1, tip.color = 0)
	par(new=T)
	plot(tree, cex=0.15, edge.col=ifelse(branch_counts_df$Edge %in% best_path_counts$Edge, yes=3, no=0), show.tip.label = T, edge.width = 1, tip.color = "grey3")
	edgelabels(edge=branch_counts_df$Edge[branch_counts_df$conflict>0],pch=20, col=alpha("red", 0.5), cex=log((branch_counts_df$conflict[branch_counts_df$conflict>0])+1)/2)
	edgelabels(edge=branch_counts_df$Edge[branch_counts_df$support>0],pch=20, col=alpha("darkgreen", 0.7),cex=log((branch_counts_df$support[branch_counts_df$support>0])+1)/2)
	edgelabels(col="darkred", frame="none",edge=branch_counts_df$Edge[branch_counts_df$conflict>0], text=branch_counts_df$conflict[branch_counts_df$conflict>0], cex=0.3)
	edgelabels(col="black", frame="none",edge=branch_counts_df$Edge[branch_counts_df$support>0], text=branch_counts_df$support[branch_counts_df$support>0], cex=0.3)
	
	edgeLen<-tree$edge.length[best_path_counts$Edge[best_path_counts$Node2==best_node]]
	estimated_loc_at_branch=edgeLen*position_in_branch


	if(position_in_branch==0){
		estimated_loc_at_branch=edgeLen*0

	    try(x <- getphylo_x(tree, getAncestors(tree,best_node)[1]))
	    try(y <- getphylo_y(tree, best_node))
   		points(x+edgeLen,y, col="black", bg=alpha("yellow", 0.3), pch=21, cex=1)
	} else if(position_in_branch>0){
		estimated_loc_at_branch=edgeLen*position_in_branch

	    # try(x <- getphylo_x(tree, getAncestors(tree,best_node)[1]))
   	 	try(x <- getphylo_x(tree, best_node))
   		try(y <- getphylo_y(tree, best_node))
   		points(x-(edgeLen-estimated_loc_at_branch),y, col="black", bg=alpha("yellow", 0.3), pch=21, cex=1)
	}
}


estimatePlotDimensions<-function(tree){
	height=dim(tree$edge)[1]/75
	width=height*(2/3)

	#in case size too small, give a minimum size to prevent issues with plotting
	if (height<5 | width<5){
		height<-5
		width<-5
	}

	sizes<-list(height, width)
	return(sizes)
}


#make full report of counts at every path
makeCountsEveryPath<-function(paths_list, branch_counts_df){
	count_all_paths_df<-matrix(ncol=8,nrow=0)
	colnames(count_all_paths_df) <-c("path",colnames(branch_counts_df))
	count_all_paths_df<-data.frame(count_all_paths_df)
	for (path_num in 1:length(paths)){
		tmp_df<-getCountsforPath(paths_list[[path_num]], branch_counts_df, "nodes")
		tmp_df$path<-path_num
		count_all_paths_df<-rbind(count_all_paths_df,tmp_df)
	}
	return(count_all_paths_df)
}

#decide best path - the one containing higher number of derived alleles.
chooseBestPath<-function(path_scores,branch_counts_df){
	best_path_number<-path_scores$path[which(path_scores$total_derived==max(path_scores$total_derived))]
	if (length(best_path_number)==1){
		best_path<-paths[[best_path_number]]
		best_path_counts<-getCountsforPath(best_path, branch_counts_df, "nodes")
		last_node_w_support<-best_path_counts$Node2[which(best_path_counts$support>0)[length(which(best_path_counts$support>0))]]
		updated_best_path<-c(best_path_counts$Node1[1],best_path_counts$Node2[1:which(best_path_counts$Node2==last_node_w_support)])
	} else if (length(best_path_number)>1){
		best_path_numbers<-best_path_number
		best_path<-as.numeric(names(which(table(unlist(paths[c(best_path_numbers)]))==length(best_path_numbers))))
		best_path_counts<-getCountsforPath(best_path, branch_counts_df, "nodes")
		last_node_w_support<-best_path_counts$Node2[which(best_path_counts$support>0)[length(which(best_path_counts$support>0))]]
		updated_best_path<-c(best_path_counts$Node1[1],best_path_counts$Node2[1:which(best_path_counts$Node2==last_node_w_support)])
	}
	return(updated_best_path)
}




# plotting
getphylo_x <- function(tree, node) {
    if(is.character(node)) {
        node <- which(c(tree$tip.label, tree$node.label)==node)
    }
    pi <- tree$edge[tree$edge[,2]==node, 1]
    if (length(pi)) {
        ei<-which(tree$edge[,1]==pi & tree$edge[,2]==node)
        tree$edge.length[ei] + Recall(tree, pi)
    } else {
        if(!is.null(tree$root.edge)) {
            tree$root.edge
        } else {
            0
        }
    }
}

getphylo_y <- function(tree, node) {
    if(is.character(node)) {
        node <- which(c(tree$tip.label, tree$node.label)==node)
    }
    ci <- tree$edge[tree$edge[,1]==node, 2]
    if (length(ci)==2) {
        mean(c(Recall(tree, ci[1]), Recall(tree, ci[2])))
    } else if (length(ci)==0) {
        Ntip <- length(tree$tip.label)
        which(tree$edge[tree$edge[, 2] <= Ntip, 2] == node)
    } else {
        stop(paste("error", length(ci)))
    }
}




estimatePositionInBranch<-function(best_path){
	best_path_counts<-getCountsforPath(best_path, branch_counts, "nodes")
	lastEdgeCounts<-best_path_counts[which(best_path_counts$support>0)[length(which(best_path_counts$support>0))],]
	# lastEdgeCounts<-best_path_counts[dim(best_path_counts)[1],]
	if (lastEdgeCounts$conflict>0 & lastEdgeCounts$support>0){
		total_counts_obs<-(lastEdgeCounts$support+lastEdgeCounts$conflict)
		position_in_branch<-lastEdgeCounts$support/total_counts_obs
	} else if (lastEdgeCounts$conflict==0 & lastEdgeCounts$support>0){
		position_in_branch<-0
	}
	return(position_in_branch)
}


makePaths<-function(tree){
	paths<-list()
	for (tip_num in 1:length(tree$tip.label)){
		paths[[tip_num]]<-c(rev(getAncestors(tree, which(tree$tip.label==tree$tip.label[tip_num]))),tip_num)
	}
	return(paths)
}


getCountsforPath<-function(path_vect, branch_counts_df, mode){
	if (mode=="nodes"){
		nodes_in_path<-data.frame(Node1=path_vect[1:length(path_vect)-1], Node2=path_vect[2:length(path_vect)])
		rel_branch_count<-merge(branch_counts_df, nodes_in_path)
		rel_branch_count<-rel_branch_count[order(rel_branch_count$Edge),]
	} else if (mode=="edges"){
		edges<-path_vect
		rel_branch_count<-branch_counts_df[match(edges,branch_counts_df$Edge),]
		rel_branch_count<-rel_branch_count[order(rel_branch_count$Edge),]
	}
	return(rel_branch_count)
}



createPathScoresTab<-function(){
	path_scores<-matrix(ncol=4,nrow=0)
	colnames(path_scores)<-c("path","total_derived","total_ancestral", "stopped_edges")
	path_scores<-data.frame(path_scores)
}





traversePaths<-function(list_of_paths,branch_counts_df, maximumTolerance){


	path_number=0
	edges_walked<-list()
	stop_edges<-NULL
	path_scores<-createPathScoresTab()

	for (path in list_of_paths){
		path_number=path_number+1

		rel_branch_count<-getCountsforPath(path, branch_counts_df, "nodes")
		rel_branch_count[order(rel_branch_count$Edge),]
		edges<-rel_branch_count$Edge
		supporting<-rel_branch_count$support
		conflicting<-rel_branch_count$conflict

		max_tol<-maximumTolerance

		derived_sum=0
		ancestral_sum=0
		edges_walked[path_number] <- edges_walked[path_number]
		stop_edge<-NULL

		for (edge_num in 1:length(edges)){
			edge<-edges[edge_num]
			# edges_walked[[path_number]]<-c(edges_walked[[path_number]],edge)
			support_count<-supporting[edge_num]
			conflict_count<-conflicting[edge_num]		
			if (edge %in% unique(stop_edges)){
				break
			} else if (conflict_count>as.numeric(max_tol)){
				stop_edge<-edge
				stop_edges<-unique(c(stop_edges, stop_edge))
				break
			} else {
				derived_sum=derived_sum+support_count
				ancestral_sum=ancestral_sum+conflict_count
				# edges_walked[[path_number]]<-c(edges_walked[[path_number]],edge)
			}
		}
		if (is.null(stop_edge)){
			stop_edge<-NA
		}
		currentPathScore<-data.frame(path=path_number, total_derived=derived_sum,total_ancestral=ancestral_sum, stopped_edges=stop_edge)
		path_scores<-rbind(path_scores, currentPathScore)
	}
	return(path_scores)
}



makeSNPStatusOutput<-function(counts_df){
	counts_df$derived_allele<-NA
	counts_df$ancestral_allele<-NA
	
	counts_df$derived_allele[counts_df$status=='-']<-as.character(counts_df$REF[counts_df$status=='-'])
	counts_df$ancestral_allele[counts_df$status=='-']<-as.character(counts_df$ALT[counts_df$status=='-'])
	counts_df$derived_allele[counts_df$status=='+']<-as.character(counts_df$ALT[counts_df$status=='+'])
	counts_df$ancestral_allele[counts_df$status=='+']<-as.character(counts_df$REF[counts_df$status=='+'])
	
	counts_df$derived_allele_count<-NA
	counts_df$ancestral_allele_count<-NA
	counts_df$derived_allele_count[counts_df$status=='-']<-counts_df$REFreads[counts_df$status=='-']
	counts_df$ancestral_allele_count[counts_df$status=='-']<-counts_df$ALTreads[counts_df$status=='-']
	counts_df$derived_allele_count[counts_df$status=='+']<-counts_df$ALTreads[counts_df$status=='+']
	counts_df$ancestral_allele_count[counts_df$status=='+']<-counts_df$REFreads[counts_df$status=='+']

	counts_df<-counts_df[c('chr','pos','marker','hg','branch','derived_allele_count','ancestral_allele','allele_status')]
	return(counts_df)
}



assignAncientCallsToBranch<-function(input_calls, sites_info){
	#merge pileup_calls with tree site info
	merged_dfs<-merge(sites_info, input_calls)
	
	#exclude missing
	merged_dfs<-merged_dfs[merged_dfs$geno!=-9,]
	
	# make tables of derived and ancestral SNP count
	# the derived SNPs are those which match the alleles which define a given branch
	# the ancestral SNPs are those which are in conflict with the alleles which define a given branch
	derived<-rbind(merged_dfs[merged_dfs$geno==0 & merged_dfs$status=='-',], merged_dfs[merged_dfs$geno==1 & merged_dfs$status=='+',])
	ancestral<-rbind(merged_dfs[merged_dfs$geno==0 & merged_dfs$status=='+',],merged_dfs[merged_dfs$geno==1 & merged_dfs$status=='-',])


	if (dim(derived)[1]>0 & dim(ancestral)[1]>0){
		derived$allele_status<-'Der'
		ancestral$allele_status<-'Anc'
	} else if (dim(derived)[1]==0 & dim(ancestral)[1]>0){
		ancestral$allele_status<-'Anc'
	} else if (dim(derived)[1]>0 & dim(ancestral)[1]==0){
		derived$allele_status<-'Der'
	}
	calls_on_branches<-list(der=derived, anc=ancestral)
	return(calls_on_branches)
}



makeBranchStatusTable<-function(counts_table,edge_table){
	table_w_derived<-counts_table[counts_table$allele_status=="Der",]
	table_w_ancestral<-counts_table[counts_table$allele_status=="Anc",]
	
	#summarize counts of derived and ancestral counts
	der_table<-table(table_w_derived$branch)
	anc_table<-table(table_w_ancestral$branch)
	
	#convert to df
	der_table<-data.frame(der_table)
	if(dim(der_table)[1]==0){
		der_table<-data.frame(branch=0, support=0)
	}
	colnames(der_table)<-c('branch', 'support')

	anc_table<-data.frame(anc_table)
	if(dim(anc_table)[1]==0){
		anc_table<-data.frame(branch=0, support=0)
	}
	colnames(anc_table)<-c('branch', 'conflict')
	
	#merge, including those branches for which there is only either ancestral or derived data
	counts<-merge(der_table,anc_table, by='branch', all=T)
	counts$branch<-as.numeric(as.character(counts$branch))
	counts<-counts[order(counts$branch),]

	colnames(counts)[1]<-"Edge"

	#get relevant info for making the branch_count table
	branch_info<-edge_table[c('Edge','Node1','Node2','hg','snp_count')]
	
	#merge, also creating a row for the edges for which no data was obs in the input sample 
	branch_counts<-(merge(branch_info,counts, all=T, by="Edge"))
	branch_counts<-branch_counts[c("Edge","Node1","Node2","support","conflict","snp_count","hg")]
	
	#replace NAs with 0. No data observed at these branches.
	branch_counts$conflict[is.na(branch_counts$conflict)]<-0
	branch_counts$support[is.na(branch_counts$support)]<-0
	branch_counts<-branch_counts[branch_counts$Edge!=0,]

	return(branch_counts)
}




makeHaplogroupStatusOutput<-function(counts_df){
	derived_hgs<-data.frame(table(counts_df$hg[!is.na(counts_df$hg) & counts_df$allele_status=="Der" ]))
	colnames(derived_hgs)<-c('hg', 'derived_count')
	ancestral_hgs<-data.frame(table(counts_df$hg[!is.na(counts_df$hg) & counts_df$allele_status=="Anc" ]))
	colnames(ancestral_hgs)<-c('hg', 'ancestral_count')
	
	merged_hgs<-merge(derived_hgs, ancestral_hgs, by='hg', all=T)
	merged_hgs$hg<-as.character(merged_hgs$hg)
	merged_hgs<-merged_hgs[order(merged_hgs$hg),]
	merged_hgs$derived_count[is.na(merged_hgs$derived_count)]<-0
	merged_hgs$ancestral_count[is.na(merged_hgs$ancestral_count)]<-0
	return(merged_hgs)
}




plotAncDerSNPTree<-function(branch_counts_table){
	support_branch_counts<-branch_counts_table[branch_counts_table$support>0,]
	conflict_branch_counts<-branch_counts_table[branch_counts_table$conflict>0,]

  	plot(tree, cex=0.1, edge.col="grey", show.tip.label=F)

  	if(dim(support_branch_counts)[1]>0){
 		edgelabels(edge=support_branch_counts$Edge,pch=20, col=alpha("darkgreen", 0.7),cex=log((support_branch_counts$support)+1)/2)
  	}
  	if(dim(conflict_branch_counts)[1]>0){
	    edgelabels(edge=conflict_branch_counts$Edge,pch=20, col=alpha("red", 0.5), cex=log((conflict_branch_counts$conflict)+1)/2)
	}

  	if(dim(support_branch_counts)[1]>0){
   		edgelabels(edge=support_branch_counts$Edge, text=support_branch_counts$support,frame = "none", cex=0.3)
  	}
  	if(dim(conflict_branch_counts)[1]>0){
   		edgelabels(edge=conflict_branch_counts$Edge, text=conflict_branch_counts$conflict,frame = "none", cex=0.3)
	}

}

get_vcf<-function(vcf_name){
	all_content = readLines(vcf_name)
	skip = all_content[-c(grep("CHROM",all_content))]
	vcf <- read.table(textConnection(skip), stringsAsFactors=F)
	header<-unlist(strsplit(all_content[grep("CHROM", all_content)[length(grep("CHROM", all_content))]], '\t'))
    colnames(vcf)<-make.names(header)
    #if T alleles read as TRUE, convert to character T.
    vcf$REF[vcf$REF==TRUE]<-"T"
    vcf$ALT[vcf$ALT==TRUE]<-"T"
	all_content<-NULL
	return(vcf)
}
