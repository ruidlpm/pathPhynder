require(phytools)

cat('\n\n',"addAncToTree.R", '\n\n\n')


tmpstr<-system('bash -l',input=c("shopt -s expand_aliases","type pathPhynder"), intern=T)

packpwd<-paste0(gsub('pathPhynder.R','',gsub('\'','',gsub('.*.Rscript ','',tmpstr))),'R')



args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=2) {
    stop("  Arguments needed.\n
        \tusage
        \tRscript addAncToTree.R <input_phylogeny.nwk> <results_folder>
        ", call.=FALSE)
}



tree_file=args[1]
results_folder=args[2]

tree<-read.tree(file=tree_file)
tree<-ladderize(tree)
tree$tip.label<-make.names(tree$tip.label)


readBestNodes<-function(){
	res<-list.files("results_folder/", patter=".best_node_info.txt")
	df<-data.frame(matrix(ncol=5, nrow=0))
	colnames(df) <- c("sample" ,"best_node" ,"position_in_branch" ,"edgeLen" ,"estimated_loc_at_branch")
	
	if (length(res)==0){
		stop("No files with best node info.")
	} else {
		for (i in res){
			tmp<-try(read.table(paste0("results_folder/",i), h=T, stringsAsFactors=F), silent=T)
			if (class(tmp)!="try-error"){
				if (dim(tmp)[1]==1){
					df<-rbind(df,tmp)
				} else {
					cat(paste0("\n",res, " has no data\n"))
				}
			}
		}
	}
	if (dim(df)[1]==0){
		stop("No best node info, cannot add samples to tree.")
	} else {
		cat(paste0("\nRead best node for ", dim(df)[1], " samples.\n"))
		return(df)
	}
}


best_nodes_table<-readBestNodes()
best_nodes_table<-best_nodes_table[order(best_nodes_table$estimated_loc_at_branch, decreasing=F),]



newtree<-tree

for (i in 1:length(best_nodes_table$best_node)){
    print(best_nodes_table[i,])
    # Get descendants of final_node in tree
    desc<-(getDescendants(tree, best_nodes_table$best_node[i]))
    desc_tip_names<-tree$tip.label[tree$tip.label %in% tree$tip.label[desc]]

    if(length(desc_tip_names)==1){
        #if just 1 descendant, then final node is tip
        final_node_in_new_tree<-which(newtree$tip.label==desc_tip_names)
    } else if(length(desc_tip_names)>1){
        #if more than 1 descendant, then need to find MRCA of all descendants
        final_node_in_new_tree<-(findMRCA(newtree, desc_tip_names))
    }
    if (best_nodes_table$estimated_loc_at_branch[i]==0){
        try(newtree<-bind.tip(tree=newtree, position = 0, edge.length=0.0001, where=final_node_in_new_tree,tip.label = best_nodes_table$sample[i]))

    } else if (best_nodes_table$estimated_loc_at_branch[i]>0){
        try(newtree<-bind.tip(tree=newtree, position =best_nodes_table$edgeLen[i]-best_nodes_table$estimated_loc_at_branch[i], edge.length=best_nodes_table$edgeLen[i]-best_nodes_table$estimated_loc_at_branch[i], where=final_node_in_new_tree, tip.label = best_nodes_table$sample[i]))
    }
}


# , height=20, width=20


pdf(file="final_tree.pdf")
plot(newtree, show.tip.label=T,tip.color=ifelse(newtree$tip.label %in% best_nodes_table$sample, yes=2, no=1), cex=0.7)
dev.off()



write.tree(newtree,file=paste0(results_folder,'/final_tree.pdf'))





















# # best_nodes_table<-read.table("results_folder/best_nodes_table.txt", h=T, stringsAsFactors=F)
# #  tmptree<-read.tree('/Users/rm890/testing/na_algo/revisiting_code/tree_sole_only.nwk')

# best_nodes_table<-best_nodes_table[order(best_nodes_table$pos_at_best_edge-best_nodes_table$total_len),]

# newtree<-tmptree
# newtree$edge.length<-tmptree$edge.length

# dup_entr<-best_nodes_table[c('best_node1' ,'best_node2' ,'pos_at_best_edge'  ,'total_len')][duplicated.data.frame(best_nodes_table[c('best_node1' ,'best_node2' ,'pos_at_best_edge'  ,'total_len')]),]

# dup_table<-merge(best_nodes_table, dup_entr)
# dup_table<-dup_table[order(dup_table$pos_at_best_edge-dup_table$total_len),]

# uniq_entr<-best_nodes_table[!best_nodes_table$sample_name %in%dup_table$sample_name ,]
# uniq_entr<-rbind(dup_table[dup_table$pos_at_best_edge==0,],uniq_entr)
# uniq_entr<-uniq_entr[order(uniq_entr$pos_at_best_edge-uniq_entr$total_len),]
# dup_table<-dup_table[dup_table$pos_at_best_edge!=0,]

# dup_table$pos_at_best_edge[which(duplicated(dup_table$pos_at_best_edge))]<-dup_table$pos_at_best_edge[which(duplicated(dup_table$pos_at_best_edge))]-0.0000000001
# uniq_entr<-rbind(dup_table,uniq_entr)
# uniq_entr<-unique(uniq_entr[order(uniq_entr$pos_at_best_edge-uniq_entr$total_len),])



# # try(best_nodes_table$hg_label<-paste(gsub(".intree.txt","",best_nodes_table$sample_name),anchgs$V4[match(gsub(".intree.txt","",best_nodes_table$sample_name), anchgs$V1)], sep="___"))

# for (i in 1:length(uniq_entr$best_node2)){
#     print(i)
#     # Get descendants of final_node in tmptree
#     desc<-(getDescendants(tmptree, uniq_entr$best_node2[i]))
#     desc_tip_names<-tmptree$tip.label[tmptree$tip.label %in% tmptree$tip.label[desc]]

#     if(length(desc_tip_names)==1){
#         #if just 1 descendant, then final node is tip
#         final_node_in_new_tree<-which(newtree$tip.label==desc_tip_names)
#     } else if(length(desc_tip_names)>1){
#         #if more than 1 descendant, then need to find MRCA of all descendants
#         final_node_in_new_tree<-(findMRCA(newtree, desc_tip_names))
#     }
#     if (uniq_entr$pos_at_best_edge[i]==0){
#         newtree<-bind.tip(tree=newtree, position = 0, edge.length=0.001, where=final_node_in_new_tree,tip.label = paste(gsub(".intree.txt","",uniq_entr$sample_name[i]),anchgs$V4[match(gsub(".intree.txt","",uniq_entr$sample_name[i]), anchgs$V1)], sep="___"))
#     } else if (uniq_entr$pos_at_best_edge[i]>0){
#         try(newtree<-bind.tip(tree=newtree, position =(uniq_entr$total_len[i]-uniq_entr$pos_at_best_edge[i]), edge.length=abs(uniq_entr$pos_at_best_edge[i]-uniq_entr$total_len[i]), where=final_node_in_new_tree, tip.label = paste(gsub(".intree.txt","",uniq_entr$sample_name[i]),anchgs$V2[match(gsub(".intree.txt","",uniq_entr$sample_name[i]), make.names(anchgs$V1))],anchgs$V4[match(gsub(".intree.txt","",uniq_entr$sample_name[i]), make.names(anchgs$V1))], sep="___")))
#     }
# }


# new_tree_no_labels<-tree

# for (i in 1:length(uniq_entr$best_node2)){
#     print(i)
#     # Get descendants of final_node in tmptree
#     desc<-(getDescendants(tree, uniq_entr$best_node2[i]))
#     desc_tip_names<-tree$tip.label[tree$tip.label %in% tree$tip.label[desc]]

#     if(length(desc_tip_names)==1){
#         #if just 1 descendant, then final node is tip
#         final_node_in_new_tree<-which(new_tree_no_labels$tip.label==desc_tip_names)
#     } else if(length(desc_tip_names)>1){
#         #if more than 1 descendant, then need to find MRCA of all descendants
#         final_node_in_new_tree<-(findMRCA(new_tree_no_labels, desc_tip_names))
#     }
#     if (uniq_entr$pos_at_best_edge[i]==0){
#         new_tree_no_labels<-bind.tip(tree=new_tree_no_labels, position = 0, edge.length=0.001, where=final_node_in_new_tree,tip.label = gsub(".intree.txt","",uniq_entr$sample_name[i]))
#     } else if (uniq_entr$pos_at_best_edge[i]>0){
#         try(new_tree_no_labels<-bind.tip(tree=new_tree_no_labels, position =(uniq_entr$total_len[i]-uniq_entr$pos_at_best_edge[i]), edge.length=abs(uniq_entr$pos_at_best_edge[i]-uniq_entr$total_len[i]), where=final_node_in_new_tree, tip.label = gsub(".intree.txt","",uniq_entr$sample_name[i])))
#     }
# }



# write.tree(new_tree_no_labels,file=paste0(args[2],'/added_anc_best_node_location_nolabels.nwk'))

