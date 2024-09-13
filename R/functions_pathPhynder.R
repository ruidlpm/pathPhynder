
tmparg <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", tmparg[grep("^--file=", tmparg)])))
# tmpstr<-system('bash -l',input=c("shopt -s expand_aliases","type pathPhynder"), intern=T)
packpwd<-paste0(gsub('functions_pathPhynder.R','',gsub('\'','',gsub('.*.Rscript ','',scriptPath))),'')

getAncestors<-phytools:::getAncestors

get_score<-function(x){
	above<-x[1:dim(x)[1]-1,]
	onbranch<-x[dim(x)[1],]
	score_string<-paste0("[",sum(above$support),"-",sum(above$conflict), ";",sum(onbranch$support),"-",sum(onbranch$conflict), "]")
	return(score_string)
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



plotBestPathTree<-function(tree,best_path_counts,branch_counts_df,path_scores_df, position_in_branch, best_node,...){
	plot(tree, cex=0.15, edge.col=ifelse(branch_counts_df$Edge %in% unique(path_scores_df$stopped_edges), yes="red", no='lightgrey'), show.tip.label = T, edge.width = 1, tip.color = 0)
	par(new=T)
	plot(tree, cex=0.15, edge.col=ifelse(branch_counts_df$Edge %in% best_path_counts$Edge, yes=3, no=0), show.tip.label = T, edge.width = 1, tip.color = "grey3")
	edgelabels(edge=branch_counts_df$Edge[branch_counts_df$conflict>0],pch=20, col=alpha("red", 0.5), cex=log((branch_counts_df$conflict[branch_counts_df$conflict>0])+1)/2)
	edgelabels(edge=branch_counts_df$Edge[branch_counts_df$support>0],pch=20, col=alpha("darkgreen", 0.7),cex=log((branch_counts_df$support[branch_counts_df$support>0])+1)/2)
	edgelabels(col="darkred", frame="none",edge=branch_counts_df$Edge[branch_counts_df$conflict>0], text=branch_counts_df$conflict[branch_counts_df$conflict>0], cex=0.3)
	edgelabels(col="black", frame="none",edge=branch_counts_df$Edge[branch_counts_df$support>0], text=branch_counts_df$support[branch_counts_df$support>0], cex=0.3)

	if(length(branch_counts_df$Edge[branch_counts_df$support>0 & branch_counts_df$conflict>0])>0){
	edgelabels(col=1, pos=3, frame="none",edge=branch_counts_df$Edge[branch_counts_df$support>0 & branch_counts_df$conflict>0], 
		text=paste0('+',branch_counts_df$support[branch_counts_df$support>0 & branch_counts_df$conflict>0], '/-',branch_counts_df$conflict[branch_counts_df$support>0 & branch_counts_df$conflict>0]), cex=0.3)

	}

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



# plotting (thanks https://stackoverflow.com/questions/25624986/draw-on-a-phylogeny-edge-with-r)
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

		if (dim(rel_branch_count)[1]==0){
			print("no SNP counts at branches")
		} else {
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


get_nested_hg<-function(in_hg){
		if (in_hg=="A0-T"){
			hg<-toupper(letters)[!toupper(letters) %in% c("F", "K",'P')]
		} else if (in_hg=="BT"){
			hg<-toupper(letters)[!toupper(letters) %in% c("A","F", "K",'P')]
		} else if (in_hg=="CT"){
			hg<-toupper(letters)[!toupper(letters) %in% c("A","B","F", "K",'P')]
		} else if (in_hg=="F"){
			hg<-c(toupper(letters[7:11]),toupper(c("m","n","o","q","r","s")))
		} else if (in_hg=="DE"){
			hg<-c('D','E')
		} else if (in_hg=="CF"){
			hg<-c('C',c(toupper(letters[7:11]),toupper(c("m","n","o","q","r","s"))))
		} else if (in_hg=="GHIJK"){
			hg<-c(toupper(letters[7:11]),toupper(c("m","n","o","q","r","s")))
		} else if (in_hg=="HIJK"){
			hg<-c(toupper(letters[8:11]),toupper(c("m","n","o","q","r","s")))
		} else if (in_hg=="IJK"){
			hg<-c(toupper(letters[9:11]),toupper(c("m","n","o","q","r","s")))
		} else if (in_hg=="IJ"){
			hg<-c('I','J')
		} else if (in_hg=="LT"){
			hg<-c('L','T')
		} else if (in_hg=="K"){
			hg<-toupper(c("m","n","o","q","r","s"))
		} else if (in_hg=="NO"){
			hg<-c('N','O')
		} else if (in_hg=="P1~orK2b2a~"){
			hg<-c('R','S')
		} else if (in_hg=="P1orK2b2a"){
			hg<-c('R','S')
		} else if (in_hg=="PorK2b2"){
			hg<-c('R','S')
		} else if (in_hg=="P"){
			hg<-c('R','S')
		} else if (in_hg=="P1"){
			hg<-c('R','S')
		} else {
			return(in_hg)
		}
		hg<-hg[!hg %in% c("U", "V", "W", "X", "Y" ,"Z")]
	return(hg)
}


trim_hg<-function(hg){
	trimmed_hg<-paste(unlist(strsplit(hg,''))[1:nchar(hg)-1], collapse="")
	return(trimmed_hg)
}



determineHG<-function(isogg_hgs_out, best_path_report){

	hgs_table <-read.table(isogg_hgs_out, h=T, sep='\t')

# new
# hgs_table<-hgs_table[!hgs_table$haplogroup %in% c('A0-T','BT','CT','F','DE','CF','GHIJK','HIJK','IJK','IJ','LT','K','NO','P','P1'),]

	#remove branches without known isogg hgs
	j<-best_path_report[best_path_report$hg!='',]

	if (dim(j)[1]==0){
		print('Not enough data for hg assignment')
	} else {
		hgs<-NULL
		hg<-NULL
		rel_snps<-NULL

		#get last line of best path (last known isogg hg)
		hgs<-unlist(strsplit(j$hg[length(j$hg)], split= ';'))

		#many branches have multiple haplogroups
		#the fact that a sample was placed in a branch
		#does not imply derived state at all SNPs which define it
		#therefore, conservatively, we choose smallest hg
		#which should also be the most ancestral

		# hg_unfilt<-hgs[which(nchar(hgs)== min(nchar(hgs)))]

# opted instead for checking tree hgs with evidence of derived markers at isogg snps
hg_unfilt<-NULL
for (i in hgs){
	tmp<-hgs_table[hgs_table$haplogroup==i,]
	if (dim(tmp)[1]==0){
		next
	} else if (sum(tmp$status=="Der")>0){
		print (i)
		print(sum(tmp$status=="Der"))
		hg_unfilt<-c(hg_unfilt,i)
	}
	
}

if (is.null(hg_unfilt)){
		hg_unfilt<-hgs[which(nchar(hgs)== min(nchar(hgs)))]
	} else {
		hg_unfilt<-hg_unfilt[order(nchar(hg_unfilt), hg_unfilt)]

	}
	# if nested macrohaplogroup CT, for example, then give letters C to T
	if (length(hg_unfilt)==1){
		hg<-get_nested_hg(hg_unfilt)
	} else if (length(hg_unfilt)>1){
		hg<-NULL
		for (indiv_hg in 1:length(hg_unfilt)){
			hg[i]<-get_nested_hg(hg_unfilt[indiv_hg])
		}
	}


	ders<-data.frame(matrix(ncol=3, nrow=0))
	colnames(ders)<-c('hg','Der','Anc')
	relevant_table<-data.frame(matrix(ncol=3, nrow=0))
	colnames(relevant_table)<-c('hg','Der','Anc')

	#loop through haplogroups to account for the possibility of more than one hg defining snp at a branch
	for (possible_hg in hg){
		
		#get lines in isogg hg table which contain the hg string
		rel_snps<-hgs_table[grep(gsub('~','', possible_hg), hgs_table$haplogroup),]

		#get only Anc and Der base count columns, excluding ambiguous sites
		rel_snps<-rel_snps[rel_snps$status %in% c('Anc', 'Der'),]

		# if Isogg SNPs found
		if (dim(rel_snps)[1]>0){
			tmp<-as.matrix(table(rel_snps$haplogroup, rel_snps$status))
			# if there isn't a Der field, then there are no derived SNPs
			# trim 1 letter of hg
			if(! 'Der' %in% rel_snps$status){
				possible_hg_trimmed<-possible_hg
				i=0
				while(! 'Der' %in% rel_snps$status){
					i=i+1
					if( possible_hg_trimmed!=''){
						possible_hg_trimmed<-trim_hg(possible_hg_trimmed)
						rel_snps<-hgs_table[grep(gsub('~','', possible_hg_trimmed),hgs_table$haplogroup),]
						rel_snps<-rel_snps[rel_snps$status %in% c('Anc', 'Der'),]
						tmp<-as.matrix(table(rel_snps$haplogroup, rel_snps$status))
					} else {
						print('not possible to go beyond hgs_along_path.txt')
						break
					}
				} 
			}
	
			if (! 'Anc' %in% rel_snps$status) {
				rel_table<-data.frame(hg=rownames(tmp),Der=tmp[,1], Anc=0 )
			} else {
				rel_table<-data.frame(hg=rownames(tmp),Der=tmp[,2], Anc= tmp[,1] )
			}
			rownames(rel_table)<-NULL
		} else {
			try(rel_table<-data.frame(hg=j[dim(j)[1],]$hg,Der=j[dim(j)[1],]$support, Anc= j[dim(j)[1],]$conflict))
			if (dim(rel_table)[1]>0){
				paste0("Haplogroup (present only in tree, not in ISOGG): ",rel_table$hg)
				df_out<-(rel_table)
			} else{
				print('Unable to determine haplogroup')
			}
		}
		relevant_table<-unique(rbind(relevant_table,rel_table))
		ders<-unique(rbind(ders,rel_table[rel_table$Der>0,]))
	}

	rel_table<-relevant_table

	queue<-rev(ders$hg)
	queue<-queue[!queue %in% c("BT","CF","CT","DE","GHIJK","HIJK","IJ","IJK","LT[K1]","LT[K1]~","LTorK1","LorK1a","NO1[K1a1]","NO1[K1a1]~","NO[K2a]","NO[K2a]~","NorK2a1a","P1orK2b2a","P1~orK2b2a~","PorK2b2","TorK1b")]

	unlikely<-NULL
	unlikely_table<-data.frame(matrix(ncol=3, nrow=0))
	colnames(unlikely_table)<-colnames(ders)

	separator<-data.frame(matrix(ncol=3, nrow=1, c('Inconsistent assignments:','', '')))
	colnames(separator)<-colnames(ders)

	miniseparator<-data.frame(matrix(ncol=3, nrow=1, c('---------------','', '')))
	colnames(miniseparator)<-colnames(ders)
	
	empty<-data.frame(matrix(ncol=3, nrow=1, 'NA'))
	colnames(empty)<-colnames(ders)

	# print(queue)

	searched<-''
	for (current_hg in queue){
		hg_name<-current_hg
		while(current_hg!=''){
			if (current_hg %in% searched){
				current_hg<-trim_hg(current_hg)
			}
			if(dim(rel_table[gsub('~','',rel_table$hg)==current_hg,] )[1]==0){
				searched<-c(searched, current_hg)
				next
			} else {
				if (sum(rel_table[gsub('~','',rel_table$hg)==current_hg,]$Anc)==0 & sum(rel_table[gsub('~','',rel_table$hg)==current_hg,]$Der)>0){
					searched<-c(searched, current_hg)
					next
				} else if (sum(rel_table[gsub('~','',rel_table$hg)==current_hg,]$Anc)>0 & sum(rel_table[gsub('~','',rel_table$hg)==current_hg,]$Der)>0){
					searched<-c(searched, current_hg)
					next
				} else if (sum(rel_table[gsub('~','',rel_table$hg)==current_hg,]$Anc)>0){
					unlikely_table<-rbind(unlikely_table,miniseparator)
					unlikely<-c(unlikely,hg_name)
					unlikely_table<-rbind(unlikely_table,rel_table[gsub('~','',rel_table$hg)==current_hg,])
					unlikely_table<-rbind(unlikely_table,rel_table[rel_table$hg==hg_name,])
					searched<-c(searched, current_hg)
					current_hg<-''
					# print(rel_table[gsub('~','',rel_table$hg)==current_hg,])
					
				}
			}
		}
		searched<-c(searched, current_hg)
	}

	##### in path ##########

	hgs_in_path<-unique(unlist(strsplit(j$hg[1:length(j$hg)], split= ';')))
	rel_snps_in_path<-hgs_table[hgs_table$haplogroup %in%  hgs_in_path,]
	rel_snps_in_path<-rel_snps_in_path[rel_snps_in_path$status %in% c('Anc', 'Der'),]
	tmptab<-as.matrix(table(rel_snps_in_path$haplogroup, rel_snps_in_path$status))
	
	if (dim(tmptab)[1]>0){
	
		if (! 'Anc' %in% colnames(tmptab)){
			rel_table_in_path<-data.frame(hg=rownames(tmptab),Der=tmptab[,1], Anc= 0 )
		} else if (! 'Der' %in% colnames(tmptab)){
			rel_table_in_path<-data.frame(hg=rownames(tmptab),Der=0, Anc= tmptab[,1] )
		} else {
			rel_table_in_path<-data.frame(hg=rownames(tmptab),Der=tmptab[,2], Anc= tmptab[,1] )
		}
	} else {
		rel_table_in_path<-data.frame(hg='',Der='', Anc= '' )
	}
	
	rownames(rel_table_in_path) <- NULL
	
	rel_table_in_path<-rel_table_in_path[match(hgs_in_path, rel_table_in_path$hg),]
	rel_table_in_path<-rel_table_in_path[!is.na(rel_table_in_path$hg),]
	rel_table_in_path<-rel_table_in_path[rel_table_in_path$Der>0,]
	rel_table_in_path<-rel_table_in_path[!(rel_table_in_path$hg %in% unlikely_table$hg),]

	likely_assignment_table<-ders[!ders$hg %in% unlikely_table$hg,]

	ordering_file=paste0(gsub("R$","",packpwd), "data/hap_order.txt")

	ordering<-read.table(ordering_file)
	likely_assignment_table<-likely_assignment_table[order(match(likely_assignment_table$hg, ordering$V1)),]
		
	if (dim(likely_assignment_table)[1]>0){
		assignment<-data.frame(matrix(ncol=3, nrow=1, c('Haplogroup:',likely_assignment_table$hg[length(likely_assignment_table$hg)], '')))
		colnames(assignment)<-colnames(ders)
	} else {
		assignment<-data.frame(matrix(ncol=3, nrow=1, c('Haplogroup:',rel_table_in_path$hg[length(rel_table_in_path$hg)], '')))
		colnames(assignment)<-colnames(ders)
	}


	# make table

	rel_table_in_path<-rel_table_in_path[!rel_table_in_path$hg %in% likely_assignment_table$hg,]

		
	likely_assignment_table$hg[which(!likely_assignment_table$hg %in% hgs_in_path)] <- paste0('+ ', likely_assignment_table$hg[which(!likely_assignment_table$hg %in% hgs_in_path)])
		
	rel_table_in_path<-cbind(rel_table_in_path,data.frame(type=rep("ISOGG & tree",length(rel_table_in_path$hg))))
	likely_assignment_table<-cbind(likely_assignment_table,data.frame(type=rep("likely",length(likely_assignment_table$hg))))

	assignment$type<-"ISOGG result"
	miniseparator$type<-""
	separator$type<-""
		if (dim(unlikely_table)[1]>0){
			unlikely_table
			unlikely_table$type<-''

			df<-rbind(rel_table_in_path,likely_assignment_table)
			# df<-rel_table_in_path
			z<-NULL
			z<-(j[!j$hg %in% gsub('\\+ ','',df$hg),])
			try(z$type<-"tree")
			z2<-data.frame(hg=z$hg, Der=z$support, Anc=z$conflict, type=z$type)
			df<-rbind(df, z2)
			df_out<-df[order(gsub('\\+ ','',df$hg)),]
			df_out<-rbind(df_out,miniseparator, assignment,miniseparator,miniseparator, separator, unlikely_table)
		} else {
			# df_out<-(rbind(unique(rel_table_in_path),unique(likely_assignment_table), miniseparator, assignment))
			df<-rbind(rel_table_in_path,likely_assignment_table)
			# df<-rel_table_in_path
			z<-NULL
			z<-(j[!j$hg %in% gsub('\\+ ','',df$hg),])
			try(z$type<-"tree")
			z2<-data.frame(hg=z$hg, Der=z$support, Anc=z$conflict, type=z$type)
			df<-rbind(df, z2)
			df_out<-df[order(gsub('\\+ ','',df$hg)),]
			df_out<-rbind(df_out, miniseparator, assignment)
		}	
	}
		# df_out<-df_out[grep(';',df_out$hg, invert=T),]
		# df_out<-df_out[order(match(gsub('\\+ ','',df_out$hg), ordering$V1)),]
	
	df_out<-df_out[df_out$type!='tree',]

	df_out<-df_out[!(df_out$Der==0 & df_out$Anc==0),]
	df_out<-df_out[df_out$type!="tree",]
	df_out$type<-NULL
	print(df_out)

	return(df_out)
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




# makeHaplogroupStatusOutput<-function(counts_df){

# 	derived_hgs<-data.frame(table(counts_df$hg[!is.na(counts_df$hg) & counts_df$allele_status=="Der" ]))
# 	colnames(derived_hgs)<-c('hg', 'derived_count')
# 	ancestral_hgs<-data.frame(table(counts_df$hg[!is.na(counts_df$hg) & counts_df$allele_status=="Anc" ]))
# 	colnames(ancestral_hgs)<-c('hg', 'ancestral_count')
	
# 	merged_hgs<-merge(derived_hgs, ancestral_hgs, by='hg', all=T)
# 	merged_hgs$hg<-as.character(merged_hgs$hg)
# 	merged_hgs<-merged_hgs[order(merged_hgs$hg),]
# 	merged_hgs$derived_count[is.na(merged_hgs$derived_count)]<-0
# 	merged_hgs$ancestral_count[is.na(merged_hgs$ancestral_count)]<-0
# 	return(merged_hgs)
# }


# makeHaplogroupStatusOutput<-function(counts_df){
# 	mat<-data.frame(matrix(ncol=2, nrow=0))
# 	df_mat<-data.frame(mat)

# 	derived_hgs<-df_mat
	
# 	ancestral_hgs<-df_mat
	
# 	if (dim(data.frame(table(counts_df$hg[!is.na(counts_df$hg) & counts_df$allele_status=="Der" ])))[1]>0){
# 		derived_hgs<-data.frame(table(counts_df$hg[!is.na(counts_df$hg) & counts_df$allele_status=="Der" ]))
# 	}
# 	if (dim(data.frame(table(counts_df$hg[!is.na(counts_df$hg) & counts_df$allele_status=="Anc" ])))[1]>0){
# 		ancestral_hgs<-data.frame(table(counts_df$hg[!is.na(counts_df$hg) & counts_df$allele_status=="Anc" ]))
# 	}
# 	colnames(derived_hgs)<-c('hg', 'derived_count')
# 	colnames(ancestral_hgs)<-c('hg', 'ancestral_count')

# 	merged_hgs<-merge(derived_hgs, ancestral_hgs, by='hg', all=T)
# 	merged_hgs$hg<-as.character(merged_hgs$hg)
# 	merged_hgs<-merged_hgs[order(merged_hgs$hg),]
# 	merged_hgs$derived_count[is.na(merged_hgs$derived_count)]<-0
# 	merged_hgs$ancestral_count[is.na(merged_hgs$ancestral_count)]<-0
# 	return(merged_hgs)
# }


makeHaplogroupStatusOutput<-function(counts_df){
	mat<-data.frame(matrix(ncol=2, nrow=0))
	df_mat<-data.frame(mat)

	derived_hgs<-df_mat
	
	ancestral_hgs<-df_mat
	counts_df<-counts_df[!is.na(counts_df$hg),]
	counts_df$hg<-as.character(counts_df$hg)
	if (dim(data.frame(table(counts_df$hg[as.character(counts_df$hg_strand)==as.character(counts_df$status) & counts_df$allele_status=="Der" ])))[1]>0){
		derived_hgs<-data.frame(table(counts_df$hg[as.character(counts_df$hg_strand)==as.character(counts_df$status) & counts_df$allele_status=="Der" ]))
	}
	if (dim(data.frame(table(counts_df$hg[as.character(counts_df$hg_strand)==as.character(counts_df$status) & counts_df$allele_status=="Anc" ])))[1]>0){
		ancestral_hgs<-data.frame(table(counts_df$hg[as.character(counts_df$hg_strand)==as.character(counts_df$status) & counts_df$allele_status=="Anc" ]))
	}

	colnames(derived_hgs)<-c('hg', 'derived_count')
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

  	plot(tree, cex=0.1, edge.col="grey", show.tip.label=T)

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
