
tmparg <- commandArgs(trailingOnly = F)
scriptPath <- normalizePath(dirname(sub("^--file=", "", tmparg[grep("^--file=", tmparg)])))
packpwd<-scriptPath

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
	for (path_num in 1:length(paths_list)){
		tmp_df<-getCountsforPath(paths_list[[path_num]], branch_counts_df, "nodes")
		tmp_df$path<-path_num
		count_all_paths_df<-rbind(count_all_paths_df,tmp_df)
	}
	return(count_all_paths_df)
}

#decide best path - the one containing higher number of derived alleles.
chooseBestPath<-function(path_scores,branch_counts_df,paths){
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


estimatePositionInBranch<-function(best_path,branch_counts_df){
	best_path_counts<-getCountsforPath(best_path, branch_counts_df, "nodes")
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
			stop_edge<-NULL

			for (edge_num in 1:length(edges)){
				edge<-edges[edge_num]
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


#produces a sortable key that orders haplogroup names in genuine ancestor-first
#(pre-order) tree order, for display purposes only.
hg_display_sort_key<-function(hg_vector){
	#pre-order backbone position for every named macrohaplogroup/paragroup label and
	#single root letter, matching the phylogeny already encoded in get_nested_hg()
	#(cross-checked against yfull.com/tree/ during development)
	backbone<-c(
		"A0-T"=1, "A"=2, "BT"=3, "B"=4, "CT"=5, "DE"=6, "D"=7, "E"=8,
		"CF"=9, "C"=10, "F"=11, "GHIJK"=11, "G"=12, "HIJK"=13, "H"=14,
		"IJK"=15, "IJ"=16, "I"=17, "J"=18, "K"=19,
		"LT"=20, "LTorK1"=20, "L"=21, "LorK1a"=21, "T"=22, "TorK1b"=22,
		"M"=23, "NO"=24, "N"=25, "NorK2a1a"=25, "O"=26,
		"P"=27, "P1"=27, "PorK2b2"=27, "P1orK2b2a"=27, "P1~orK2b2a~"=27,
		"Q"=28, "R"=29, "S"=30
	)

	get_key<-function(hg){
		#exact match against a known macro/paragroup label (checked on the raw string,
		#since some of these labels have their own literal tildes, e.g. "P1~orK2b2a~")
		if (hg %in% names(backbone)){
			root_pos<-unname(backbone[hg])
			tokens<-character(0)
		} else {
			tokens<-regmatches(hg, gregexpr("[0-9]+|[A-Za-z]+|[^0-9A-Za-z]+", hg))[[1]]
			if (length(tokens)==0) return("999")
			root_letter<-toupper(tokens[1])
			root_pos<-unname(backbone[root_letter])
			if (is.na(root_pos)) root_pos<-999
			tokens<-tokens[-1]
		}
		#zero-pad numeric tokens so e.g. "10" sorts after "9", not between "1" and "2"
		padded<-vapply(tokens, function(t){
			if (grepl("^[0-9]+$", t)) sprintf("n%05d", as.numeric(t)) else t
		}, character(1))
		paste0(sprintf("%03d", root_pos), "_", paste(padded, collapse="_"))
	}

	vapply(hg_vector, get_key, character(1), USE.NAMES = FALSE)
}


#Determines the sample's haplogroup by walking the tree's own resolved best
#path (best_path_report), from the tip upward, looking for the deepest edge
#backed by unambiguous derived evidence - i.e. der count >= min_der and
#der/(der+anc) >= conflict_ratio, the same ratio-based tolerance philosophy as
#maximumTolerance elsewhere in this pipeline. This replaces the previous
#hap_order.txt string-ranking approach, which had no notion of the tree's real
#topology and could rank a poorly-supported candidate above a better-supported
#one whenever the file happened to list it first.
#
#Once the deepest tree-resolved edge is found, it tries to sharpen the answer
#by searching the full ISOGG marker table for names that extend it further
#(not every real branch point is resolved as its own edge in the input tree).
#A deeper name is only adopted if it is the unique, unambiguous winner: when
#two or more deeper candidates diverge into different, mutually-exclusive
#sub-lineages (neither's name is a prefix of the other's) rather than forming
#one consistent chain, adopting either would overclaim a precision the data
#doesn't actually support, so the search stops at the deepest point they all
#agree on.
determineHG<-function(isogg_hgs_out, best_path_report, min_der=1, conflict_ratio=0.7, min_der_absolute=5, edge_df=NULL){

	hgs_table<-read.table(isogg_hgs_out, h=T, sep='\t')
	hgs_table<-hgs_table[hgs_table$status %in% c('Anc', 'Der'),]

	#remove branches without known isogg hgs
	j<-best_path_report[best_path_report$hg!='',]

	empty_result<-function(msg){
		print(msg)
		list(summary=data.frame(Haplogroup='Haplogroup:', Der_SNP_count='undetermined', Anc_SNP_count=''),
		     marker_detail=data.frame(haplogroup=character(0), POS=integer(0), marker=character(0), mutation=character(0), ANC_count=integer(0), DER_count=integer(0), status=character(0)))
	}

	if (dim(j)[1]==0){
		return(empty_result('Not enough data for hg assignment'))
	}

	#aggregate Der/Anc read counts per exact ISOGG haplogroup label, across
	#every marker tested for that label in this sample
	tmp<-as.matrix(table(hgs_table$haplogroup, hgs_table$status))
	der_counts<-setNames(rep(0, nrow(tmp)), rownames(tmp))
	anc_counts<-setNames(rep(0, nrow(tmp)), rownames(tmp))
	if ('Der' %in% colnames(tmp)) der_counts<-setNames(tmp[,'Der'], rownames(tmp))
	if ('Anc' %in% colnames(tmp)) anc_counts<-setNames(tmp[,'Anc'], rownames(tmp))
	all_hg_names<-rownames(tmp)
	all_hg_norm<-gsub('~$', '', all_hg_names)
	all_hg_nchar<-nchar(all_hg_norm)

	get_der_anc<-function(lab){
		d<-if (lab %in% names(der_counts)) unname(der_counts[lab]) else 0
		a<-if (lab %in% names(anc_counts)) unname(anc_counts[lab]) else 0
		list(der=d, anc=a, ratio=if ((d+a)>0) d/(d+a) else 0)
	}

	#a branch's defining SNPs didn't all arise in one instant - the branch
	#length represents real elapsed time, and its markers accumulated one by
	#one across it. An ancient individual can genuinely sit partway along a
	#long branch: derived at the earlier-arising markers, still ancestral at
	#the later ones (which hadn't happened yet in their lifetime, even though
	#every modern descendant - who by definition survived to the branch's
	#endpoint - now carries all of them). That's a real biological pattern,
	#not noise, so a candidate also qualifies on a large absolute count of
	#independently-agreeing derived markers, even with a poor ratio, since
	#e.g. 13 independent markers agreeing is strong evidence on its own
	#regardless of how many later-arising ones read ancestral
	qualifies<-function(d, a, r) (d>=min_der & r>=conflict_ratio) | (d>=min_der_absolute)

	norm_hg<-function(hg) gsub('~$', '', hg)

	#which entries of `candidates` (vectorized) are genuine one-or-more-level-
	#deeper ISOGG extensions of `base` - not just strings that happen to share
	#`base` as a literal prefix? ISOGG names strictly alternate letter/digit
	#tokens, so a real extension's first added character must switch type
	#(letter->digit or digit->letter) from base's own last character. This
	#matters most when base is a bare macrohaplogroup letter (e.g. "I"): plain
	#startsWith() can't distinguish a true child ("I1", letter->digit, valid)
	#from an unrelated name that merely starts with the same letter by
	#coincidence ("IJK", "IJ" - both real ANCESTORS of "I" via a
	#macrohaplogroup letter-jump, letter->letter, not a valid extension at all)
	real_extensions<-function(base, candidates){
		nb<-norm_hg(base); nc<-gsub('~$', '', candidates)
		is_prefix<-startsWith(nc, nb) & nchar(nc)>nchar(nb)
		base_last_is_digit<-grepl('[0-9]', substr(nb, nchar(nb), nchar(nb)))
		next_is_digit<-grepl('[0-9]', substr(nc, nchar(nb)+1, nchar(nb)+1))
		is_prefix & (base_last_is_digit!=next_is_digit)
	}

	#two haplogroup names are "compatible" if one is a nomenclature prefix of
	#the other (i.e. literally ancestor/descendant on the same named lineage),
	#as opposed to diverging into different, mutually-exclusive sub-lineages
	compatible<-function(a, b){
		na<-norm_hg(a); nb<-norm_hg(b)
		startsWith(na, nb) || startsWith(nb, na)
	}

	#depth is computed on the normalized (tilde-stripped) name - a trailing "~"
	#marks the same nominal tree position as its bare counterpart, not a level
	#of its own. Counting it as a token made e.g. "I2a1b1a1a1~" (11 tokens)
	#tie with the genuinely one-level-deeper "I2a1b1a1a1b" (also 11 tokens,
	#since "b" replaces "~"), and which.max()'s first-match tie-break then
	#kept the shallower tilde name over a real, uncontradicted, deeper answer
	hg_depth<-function(hg) length(regmatches(norm_hg(hg), gregexpr("[0-9]+|[A-Za-z]+|[^0-9A-Za-z]+", norm_hg(hg)))[[1]])

	#among a set of qualifying candidate labels, return the deepest one that is
	#compatible with every other qualifying candidate - i.e. the deepest point
	#not contradicted by another, equally-supported but divergent candidate.
	#Returns NULL if every candidate conflicts with at least one other.
	#cand_der (optional, parallel to cand_hg) lets a divergent competitor be
	#ignored when it's clearly weak (a single read) and the candidate under
	#test is comfortably stronger (>=3 reads) - a lone read on an unrelated
	#branch is far more likely a miscall than a real competing lineage, and
	#shouldn't be enough to veto an otherwise well-supported candidate. This
	#is intentionally asymmetric: a weak candidate is still correctly blocked
	#by a strong one, it just can't block a strong one in return.
	safest_deepest<-function(cand_hg, cand_der=NULL){
		if (length(cand_hg)==0) return(NULL)
		if (is.null(cand_der)) cand_der<-rep(NA_real_, length(cand_hg))
		is_safe<-vapply(seq_along(cand_hg), function(i){
			others<-setdiff(seq_along(cand_hg), i)
			all(vapply(others, function(j){
				if (compatible(cand_hg[i], cand_hg[j])) return(TRUE)
				!is.na(cand_der[j]) && !is.na(cand_der[i]) && cand_der[j]==1 && cand_der[i]>=3
			}, logical(1)))
		}, logical(1))
		safe<-cand_hg[is_safe]
		if (length(safe)==0) return(NULL)
		safe[which.max(vapply(safe, hg_depth, integer(1)))]
	}

	#---- walk the resolved best-path from the tip upward, edge by edge,
	#looking for the deepest edge with a safe, unambiguous qualifying label ----
	report<-data.frame(hg=character(0), Der=numeric(0), Anc=numeric(0))
	tree_hg<-NULL
	resolved_row_i<-NULL
	for (i in nrow(j):1){
		labels<-unlist(strsplit(j$hg[i], split=';'))
		da<-lapply(labels, get_der_anc)
		der<-vapply(da, function(x) x$der, numeric(1))
		anc<-vapply(da, function(x) x$anc, numeric(1))
		ratio<-vapply(da, function(x) x$ratio, numeric(1))
		report<-rbind(report, data.frame(hg=labels, Der=der, Anc=anc))
		qual_mask<-qualifies(der, anc, ratio)
		qual_labels<-labels[qual_mask]
		if (length(qual_labels)>0 && is.null(tree_hg)){
			winner<-safest_deepest(qual_labels, der[qual_mask])
			if (!is.null(winner)){
				tree_hg<-winner
				resolved_row_i<-i
				break
			}
			#the deepest edge's own qualifying labels diverge into mutually-
			#exclusive lineages with no safe resolution here - keep walking up
			#towards the root for an unambiguous qualifying edge
		}
	}

	if (is.null(tree_hg)){
		return(empty_result('Unable to determine haplogroup: no tree edge along the best path has unambiguous derived support'))
	}

	#---- try to sharpen tree_hg by searching the full ISOGG table for names
	#that extend it further than the input tree itself resolves, adopting a
	#deeper name only if it is the unique, unambiguous winner ----
	base_norm<-norm_hg(tree_hg)
	deeper_candidates<-all_hg_names[real_extensions(tree_hg, all_hg_names) & all_hg_names!=tree_hg]
	downstream_conflict<-data.frame(hg=character(0), Der=numeric(0), Anc=numeric(0))
	rejected_contradicted<-data.frame(hg=character(0), Der=numeric(0), Anc=numeric(0), reason=character(0))
	if (length(deeper_candidates)>0){
		deeper_da<-lapply(deeper_candidates, get_der_anc)
		deeper_der<-vapply(deeper_da, function(x) x$der, numeric(1))
		deeper_anc<-vapply(deeper_da, function(x) x$anc, numeric(1))
		deeper_ratio<-vapply(deeper_da, function(x) x$ratio, numeric(1))
		report<-rbind(report, data.frame(hg=deeper_candidates, Der=deeper_der, Anc=deeper_anc))
		qualifying_deeper<-deeper_candidates[qualifies(deeper_der, deeper_anc, deeper_ratio) & deeper_der>0]
		#reject a candidate outright if any position strictly between the
		#tree-grounded base and the candidate itself reads confirmed ancestral
		#(Der==0, Anc>0). On a haploid, non-recombining Y chromosome a sample
		#can't be derived at a deep SNP while ancestral at one that has to
		#precede it in the same lineage; a single derived read past a real
		#ancestral chain like that is far more likely a miscalled/recurrent
		#SNP than genuine evidence, so this isn't treated as a borderline call.
		#Only NON-tilde (ISOGG-confirmed-position) markers count as
		#contradiction evidence here - a "~" marker is ISOGG's own uncertain-
		#position designation, and it disagreeing with a confirmed marker
		#elsewhere at the "same" name is expected, not an error. A single
		#confirmed non-tilde ancestral read is sufficient to reject a
		#candidate.
		is_tilde<-grepl('~$', deeper_candidates)
		#for each qualifying candidate, the specific contradicting position(s)
		#(if any) - NA when not contradicted. Captured (not just a TRUE/FALSE)
		#so the rejection reason can be surfaced in the report instead of
		#silently disappearing: the contradicting marker is very often at the
		#SAME normalized position as the rejected candidate (e.g. a confirmed
		#"X" reading ancestral vs its own "X~" reading derived), which is
		#exactly the case the display's neighbour-injection step (below)
		#deliberately excludes - by design, for a different reason (stopping
		#the *accepted* answer's own synonym from appearing as a fake new
		#conflict) - so without capturing it here, this reason is invisible
		#anywhere in the output even though it's the actual cause of rejection
		contradiction_reason<-vapply(qualifying_deeper, function(cand){
			cn<-norm_hg(cand)
			mask<-!is_tilde & deeper_candidates!=cand & startsWith(cn, deeper_candidates) & nchar(deeper_candidates)>nchar(base_norm)
			bad<-mask & deeper_der==0 & deeper_anc>0
			if (!any(bad)) return(NA_character_)
			paste(sprintf("%s (%dder/%danc)", deeper_candidates[bad], deeper_der[bad], deeper_anc[bad]), collapse="; ")
		}, character(1))
		is_contradicted<-!is.na(contradiction_reason)
		rejected_contradicted<-data.frame(
			hg=qualifying_deeper[is_contradicted],
			Der=deeper_der[match(qualifying_deeper[is_contradicted], deeper_candidates)],
			Anc=deeper_anc[match(qualifying_deeper[is_contradicted], deeper_candidates)],
			reason=contradiction_reason[is_contradicted]
		)
		qualifying_deeper<-qualifying_deeper[!is_contradicted]
		winner<-safest_deepest(qualifying_deeper, deeper_der[match(qualifying_deeper, deeper_candidates)])
		if (!is.null(winner)) tree_hg<-winner
		#downstream candidates with real derived signal (Der>0) that were NOT
		#adopted as the final answer - either because they didn't clear the
		#conflict_ratio bar on their own, or because they conflicted with
		#another, equally-supported but divergent candidate (safest_deepest
		#returned NULL/a different winner). Surfaced explicitly here so a
		#reader can see e.g. "derived at one defining SNP, but ancestral at N
		#others on the same branch" without having to scan the full table -
		#see isogg_marker_detail.txt for the individual SNPs behind each count.
		#Excludes ancestors-or-self of the final tree_hg (e.g. a bare vs
		#tilde-suffixed variant of the same node) - those are already part of
		#the accepted chain shown above, not a real conflict.
		is_ancestor_or_self<-vapply(deeper_candidates, function(c) startsWith(norm_hg(tree_hg), norm_hg(c)), logical(1))
		conflict_mask<-deeper_der>0 & deeper_candidates!=tree_hg & !is_ancestor_or_self
		if (any(conflict_mask)){
			downstream_conflict<-data.frame(hg=deeper_candidates[conflict_mask], Der=deeper_der[conflict_mask], Anc=deeper_anc[conflict_mask], role='candidate')
			#for context, also surface each conflicting candidate's own direct
			#parent (even though it has no Der of its own) whenever that parent
			#has any recorded data at all - e.g. a candidate showing der=1/anc=1
			#is easier to read once you can see its parent already read
			#der=0/anc=2, i.e. the conflict was building before this exact node -
			#and likewise its own direct children, if any have recorded data,
			#since a conflicting candidate's local neighbourhood (not just the
			#accepted backbone's) is what a reader needs to judge it
			find_parent<-function(hg){
				nc<-norm_hg(hg)
				is_ancestor<-all_hg_norm!=nc & startsWith(nc, all_hg_norm)
				if (!any(is_ancestor)) return(NULL)
				all_hg_names[is_ancestor][which.max(all_hg_nchar[is_ancestor])]
			}
			find_children<-function(hg){
				nc<-norm_hg(hg); nchar_c<-nchar(nc)
				is_descendant<-all_hg_norm!=nc & startsWith(all_hg_norm, nc)
				if (!any(is_descendant)) return(character(0))
				cand<-all_hg_names[is_descendant]; cand_norm<-all_hg_norm[is_descendant]; cand_nchar<-all_hg_nchar[is_descendant]
				is_direct<-vapply(seq_along(cand), function(i) !any(all_hg_norm!=cand_norm[i] & all_hg_norm!=nc &
					all_hg_nchar>nchar_c & all_hg_nchar<cand_nchar[i] &
					startsWith(cand_norm[i], all_hg_norm) & startsWith(all_hg_norm, nc)), logical(1))
				cand[is_direct]
			}
			#tag each injected neighbour with its relationship to the actual
			#candidate it's context for - "above the branch" (ancestor) or
			#"below the branch" (descendant) - so a reader doesn't have to
			#infer it purely from tree-order position in the printed table
			neighbour_rows<-lapply(unique(downstream_conflict$hg), function(h){
				p<-find_parent(h)
				kids<-find_children(h)
				near<-c(if (!is.null(p)) p else character(0), kids)
				near_role<-c(if (!is.null(p)) 'ancestor' else character(0), rep('descendant', length(kids)))
				#exclude anything at the SAME normalized position as the
				#accepted answer (e.g. a bare vs tilde-suffixed synonym) - not
				#just an exact string match - that position is already part
				#of the accepted chain shown above, not a new conflict
				keep<-norm_hg(near)!=norm_hg(tree_hg) & !(near %in% downstream_conflict$hg)
				near<-near[keep]; near_role<-near_role[keep]
				if (length(near)==0) return(NULL)
				rows<-lapply(seq_along(near), function(i){ nda<-get_der_anc(near[i]); data.frame(hg=near[i], Der=nda$der, Anc=nda$anc, role=near_role[i]) })
				rows<-do.call(rbind, rows)
				rows[rows$Der>0 | rows$Anc>0,]
			})
			neighbour_rows<-neighbour_rows[!vapply(neighbour_rows, is.null, logical(1))]
			if (length(neighbour_rows)>0) downstream_conflict<-unique(rbind(downstream_conflict, do.call(rbind, neighbour_rows)))
			downstream_conflict<-downstream_conflict[order(hg_display_sort_key(downstream_conflict$hg)),]
		}
	}

	#---- topology-aware note: the reach-sharpening search above is a
	#string-prefix search, so it's blind to a macrohaplogroup letter-jump
	#(e.g. "F" -> "I" - a real parent/child pair in the tree, but "I" isn't
	#a string extension of "F"). If tree_hg is stuck at a real, named tree
	#node (edge_df$hg literally bundles it), its TRUE children are known
	#regardless of naming scheme. Any of them showing real derived evidence
	#is a genuine sign a deeper placement likely exists - just not
	#reachable without a full tree-walk past this point (a separate, known
	#limitation - see the chooseBestPath tie-collapse issue). Purely
	#informational: this never changes tree_hg itself, since confirming a
	#letter-jump candidate properly requires the real walk, not a hint
	viable_deeper_path<-NULL
	if (!is.null(edge_df)){
		own_edge<-edge_df[vapply(strsplit(edge_df$hg, ';'), function(labs) tree_hg %in% labs, logical(1)) & edge_df$hg!='', ]
		if (dim(own_edge)[1]>0){
			tree_node<-own_edge$Node2[1]
			#explore the subtree below tree_node, but PRUNE any branch the
			#moment it hits confirmed-ancestral evidence (same rule
			#is_contradicted uses above: any bundled label reading der==0 &
			#anc>0) - a real macrohaplogroup split is very often preceded by
			#a short run of near-zero-length, unresolved edges (e.g. F ->
			#"F;GHIJK" -> HIJK -> IJK -> IJ -> I -> ...), and most tree edges
			#for any one low-coverage ancient sample simply have zero reads
			#either way, so following only a single still-viable child stalls
			#immediately - it can't tell "the true path, just uncovered yet"
			#from "a genuine sibling, also uncovered yet". Exploring every
			#not-yet-disproven branch instead, and letting real ancestral
			#evidence prune away the wrong ones as it's encountered, avoids
			#pulling in noise from unrelated macrohaplogroups: an actually-
			#unrelated branch will eventually hit its own disconfirming
			#marker and get dropped, while capped total work below bounds
			#the worst case (a sample with almost no informative reads
			#anywhere)
			collected<-data.frame(hg=character(0), Der=numeric(0), Anc=numeric(0))
			frontier<-tree_node; visited<-tree_node; budget<-500
			while (length(frontier)>0 && budget>0){
				kids<-edge_df[edge_df$Node1 %in% frontier & edge_df$hg!='' & !(edge_df$Node2 %in% visited), ]
				if (dim(kids)[1]==0) break
				budget<-budget-nrow(kids)
				next_frontier<-c()
				for (i in seq_len(nrow(kids))){
					labs<-strsplit(kids$hg[i], ';')[[1]]
					labs<-labs[labs!='']
					if (length(labs)==0) next
					da<-lapply(labs, get_der_anc)
					der<-vapply(da, function(x) x$der, numeric(1))
					anc<-vapply(da, function(x) x$anc, numeric(1))
					if (any(der==0 & anc>0)) next #pruned - confirmed ancestral somewhere on this edge
					keep<-der>0 & norm_hg(labs)!=norm_hg(tree_hg)
					if (any(keep)) collected<-rbind(collected, data.frame(hg=labs[keep], Der=der[keep], Anc=anc[keep]))
					next_frontier<-c(next_frontier, kids$Node2[i])
				}
				visited<-c(visited, next_frontier)
				frontier<-next_frontier
			}
			if (dim(collected)[1]>0) viable_deeper_path<-collected[order(-collected$Der),]
		}
	}

	#---- root-to-tip backbone: re-derive the same per-edge winning label from
	#the root down to the edge actually resolved above (e.g. CT, GHIJK, G,
	#...), so the full accepted chain is visible even when the main walk
	#above jumped straight to a deep, multi-label-bundled edge without
	#needing to separately visit every shallower edge on the way there ----
	backbone<-data.frame(hg=character(0), Der=numeric(0), Anc=numeric(0))
	if (!is.null(resolved_row_i) && resolved_row_i>1){
		for (i in 1:(resolved_row_i-1)){
			labels<-unlist(strsplit(j$hg[i], split=';'))
			if (length(labels)==0) next
			da<-lapply(labels, get_der_anc)
			der<-vapply(da, function(x) x$der, numeric(1))
			anc<-vapply(da, function(x) x$anc, numeric(1))
			ratio<-vapply(da, function(x) x$ratio, numeric(1))
			winner<-safest_deepest(labels[qualifies(der, anc, ratio)])
			#a polytomy can bundle the same already-resolved label onto more
			#than one consecutive edge (e.g. edge "F" then edge "F;GHIJK") -
			#skip re-adding it if it's identical to the row just added
			if (!is.null(winner)){
				if (dim(backbone)[1]==0 || backbone$hg[dim(backbone)[1]]!=winner){
					wda<-get_der_anc(winner)
					backbone<-rbind(backbone, data.frame(hg=winner, Der=wda$der, Anc=wda$anc))
				}
			} else {
				last_hg<-if (dim(backbone)[1]>0) backbone$hg[dim(backbone)[1]] else NA
				keep<-labels!=last_hg | is.na(last_hg)
				if (any(keep)) backbone<-rbind(backbone, data.frame(hg=labels[keep], Der=der[keep], Anc=anc[keep]))
			}
		}
	}

	#---- assemble the report: the accepted root-to-tip backbone only - every
	#haplogroup that is an ancestor of, or equal to, the final answer, with any
	#recorded data at all (Der>0 or Anc>0). Once the walk resolves to tree_hg,
	#the display stops there too: candidates beyond it (whether zero-signal
	#siblings/children or real conflicts) don't belong in this table - a real
	#conflict is already shown in its own section below, and a zero-signal
	#branch adds nothing once we've already stopped short of it. Ancestors are
	#kept even when Der==0 (e.g. an intermediate level the reach-sharpening
	#step passed through that only had ancestral reads of its own, right
	#before the specific deeper SNP that tipped the answer) - that's still
	#useful "it read ancestral right up to this exact point" context, unlike a
	#rejected sibling/descendant, which is noise once rejected.
	#Note: the ancestor-of-tree_hg check below is a plain string-prefix test,
	#which only holds within one nomenclature family (e.g. R1b1a1b1a1a2c1a is
	#a real string-extension of R1b1a1b1a1a2) - it does NOT hold across a
	#macrohaplogroup letter change (e.g. "I" is a true ancestor of "IJ" in the
	#real tree, but isn't a string-prefix match). `backbone` rows don't need
	#this check at all - they were already derived from a real positional walk
	#of the actual best-path edges, so they're genuine ancestors by
	#construction regardless of naming scheme; only `report` (built from
	#tree_hg's own single-family deeper-reach search) needs it ----
	is_ancestor_of_final<-vapply(report$hg, function(h) startsWith(norm_hg(tree_hg), norm_hg(h)), logical(1))
	report<-report[(report$Der>0 | report$Anc>0) & is_ancestor_of_final,]
	report<-unique(rbind(backbone, report))
	report<-report[order(hg_display_sort_key(report$hg)),]

	#---- flag (don't act on) any strictly-shallower position on the accepted
	#chain that reads confirmed ancestral (Der==0, Anc>0) once its own
	#same-position name variants - e.g. a bare vs tilde-suffixed synonym for
	#the identical node - are pooled together. This is a real biological
	#tension worth a human's judgment, not grounds for silently rejecting the
	#deeper call: on a haploid, non-recombining Y chromosome, being derived at
	#a deep SNP while ancestral at one that has to precede it in the same
	#lineage is usually a contradiction, but not always - recurrent mutation
	#or a taxonomic quirk in how ISOGG nests a subclade can produce a real,
	#well-supported instance of exactly this pattern (seen in practice: a
	#deeper call backed by 11 reads across 2 SNPs, contradicting a shallower
	#position backed by 0/many reads across 12 SNPs, that still turned out to
	#be the correct call) ----
	ancestor_contradiction<-NULL
	final_norm<-norm_hg(tree_hg)
	chain_norm<-norm_hg(report$hg)
	same_lineage<-report[chain_norm!=final_norm & startsWith(final_norm, chain_norm),]
	if (dim(same_lineage)[1]>0){
		same_lineage$norm<-norm_hg(same_lineage$hg)
		agg<-aggregate(cbind(Der,Anc) ~ norm, data=same_lineage, sum)
		contra<-agg[agg$Der==0 & agg$Anc>0,]
		if (dim(contra)[1]>0){
			#label each contradicting position with a "~" if every row that fed
			#it was itself tilde-suffixed (ISOGG's own uncertain-position
			#marker) - don't silently normalize that away, since it changes
			#how much weight the evidence deserves (see the rejection rule
			#above, which treats tilde-only contradiction very differently)
			all_tilde<-vapply(contra$norm, function(n) all(grepl('~$', same_lineage$hg[same_lineage$norm==n])), logical(1))
			ancestor_contradiction<-data.frame(hg=paste0(contra$norm, ifelse(all_tilde, '~', '')), Der=contra$Der, Anc=contra$Anc)
		}
	}

	miniseparator<-data.frame(hg='---------------', Der='', Anc='')
	assignment<-data.frame(hg='Haplogroup:', Der=tree_hg, Anc='')

	#insert a blank row between consecutive entries that aren't on the same
	#nested lineage (i.e. a new, unrelated path is starting) - rows are
	#assumed already sorted ancestor-first, so members of the same lineage
	#are already adjacent and only a same-vs-different check against the
	#immediately preceding row is needed
	add_path_breaks<-function(df){
		if (dim(df)[1]<=1) return(df)
		blank<-data.frame(hg='', Der='', Anc='')
		out<-df[1,,drop=FALSE]
		for (i in 2:dim(df)[1]){
			if (!compatible(df$hg[i], df$hg[i-1])) out<-rbind(out, blank)
			out<-rbind(out, df[i,,drop=FALSE])
		}
		out
	}

	df_out<-rbind(report, miniseparator, assignment)
	if (!is.null(ancestor_contradiction)){
		warning_header<-data.frame(hg='WARNING - ancestral above the assigned haplogroup:', Der='', Anc='')
		df_out<-rbind(df_out, miniseparator, warning_header, ancestor_contradiction)
	}
	if (dim(downstream_conflict)[1]>0){
		conflict_header<-data.frame(hg='Downstream conflicts:', Der='', Anc='')
		#compute path-breaks on the plain hg values first (compatible() relies
		#on norm_hg's trailing-"~" stripping, which a tag appended after it
		#would break), then fold the role into the displayed name afterward -
		#the real candidate (Der/Anc evidence of its own) stays unmarked since
		#it's the focus of the row group, while its injected ancestor/
		#descendant context neighbours get an explicit tag, so a reader
		#doesn't have to infer "above/below the branch" purely from position
		conflict_rows<-add_path_breaks(downstream_conflict[,c('hg','Der','Anc')])
		role_tag<-c(candidate='', ancestor=' [ancestor, context]', descendant=' [descendant, context]')
		role_lookup<-setNames(downstream_conflict$role, downstream_conflict$hg)
		conflict_rows$hg<-ifelse(conflict_rows$hg=='', conflict_rows$hg, paste0(conflict_rows$hg, role_tag[role_lookup[conflict_rows$hg]]))
		df_out<-rbind(df_out, miniseparator, conflict_header, conflict_rows)
	}
	if (dim(rejected_contradicted)[1]>0){
		#candidates that otherwise qualified (enough derived support of their
		#own) but were rejected because a confirmed non-tilde marker at the
		#same or an intermediate position reads ancestral - previously
		#invisible in the report even though it's the actual reason a
		#seemingly-good candidate didn't become the answer
		rejected_header<-data.frame(hg='Rejected (ancestral conflict at same/intermediate position):', Der='', Anc='')
		rejected_rows<-data.frame(hg=paste0(rejected_contradicted$hg, '  <-contradicted by- ', rejected_contradicted$reason), Der=rejected_contradicted$Der, Anc=rejected_contradicted$Anc)
		df_out<-rbind(df_out, miniseparator, rejected_header, rejected_rows)
	}
	if (!is.null(viable_deeper_path) && dim(viable_deeper_path)[1]>0){
		viable_header<-data.frame(hg='NOTE: a deeper path may be viable via true tree topology:', Der='', Anc='')
		df_out<-rbind(df_out, miniseparator, viable_header, viable_deeper_path)
	}
	rownames(df_out)<-NULL

	#---- per-marker breakdown behind every haplogroup shown in the summary
	#above, so a Der/Anc count can be traced back to the individual defining
	#SNPs (and their read counts) that produced it ----
	marker_detail<-hgs_table[hgs_table$haplogroup %in% unique(c(report$hg, downstream_conflict$hg)), c('haplogroup','POS','marker','mutation','ANC_count','DER_count','status')]
	marker_detail<-marker_detail[order(hg_display_sort_key(marker_detail$haplogroup), marker_detail$POS),]
	rownames(marker_detail)<-NULL
	colnames(marker_detail)<-c('haplogroup','POS','marker','mutation','Anc_read_count','DER_read_count','status')

	colnames(df_out)<-c('Haplogroup','Der_SNP_count','Anc_SNP_count')

	return(list(summary=df_out, marker_detail=marker_detail))
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
	mat<-data.frame(matrix(ncol=2, nrow=0))
	df_mat<-data.frame(mat)

	derived_hgs<-df_mat

	ancestral_hgs<-df_mat
	counts_df<-counts_df[!is.na(counts_df$hg),]
	counts_df$hg<-as.character(counts_df$hg)

	#a marker's Anc/Der call needs flipping when the tree's REF/ALT direction is the
	#opposite of ISOGG's ancestral/derived direction (hg_strand=='-'); markers explicitly
	#flagged as 'mismatch' (or with no hg_strand at all) are excluded from both counts
	same_orientation<-as.character(counts_df$hg_strand)==as.character(counts_df$status)
	flipped_orientation<-(counts_df$hg_strand=='+' & counts_df$status=='-') | (counts_df$hg_strand=='-' & counts_df$status=='+')

	if (dim(data.frame(table(counts_df$hg[which((same_orientation & counts_df$allele_status=="Der") | (flipped_orientation & counts_df$allele_status=="Anc"))])))[1]>0){
		derived_hgs<-data.frame(table(counts_df$hg[which((same_orientation & counts_df$allele_status=="Der") | (flipped_orientation & counts_df$allele_status=="Anc"))]))
	}
	if (dim(data.frame(table(counts_df$hg[which((same_orientation & counts_df$allele_status=="Anc") | (flipped_orientation & counts_df$allele_status=="Der"))])))[1]>0){
		ancestral_hgs<-data.frame(table(counts_df$hg[which((same_orientation & counts_df$allele_status=="Anc") | (flipped_orientation & counts_df$allele_status=="Der"))]))
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






plotAncDerSNPTree<-function(branch_counts_table,tree){
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
