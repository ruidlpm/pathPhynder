suppressWarnings(suppressPackageStartupMessages(library("phytools")))

cat('\n\n',"addAncToTree.R", '\n\n\n')

tmpstr<-system('bash -l',input=c("shopt -s expand_aliases","type pathPhynder"), intern=T)
packpwd<-paste0(gsub('pathPhynder.R','',gsub('\'','',gsub('.*.Rscript ','',tmpstr))),'R')

source(paste0(packpwd,'/functions_pathPhynder.R'))

args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=3) {
    stop("  Arguments needed.\n
        \tusage
        \tRscript addAncToTree.R <input_phylogeny.nwk> <results_folder> prefix
        ", call.=FALSE)
}

tree_file=args[1]
results_folder=args[2]

tree<-read.tree(file=tree_file)
# tree<-ladderize(tree)
tree$tip.label<-make.names(tree$tip.label)


#IMPORTANT
#if the tree has no edge lengths, then add the SNP count at each branch as the edge.length 
#need to get edge lengths if tree has no edge lengths originally
edge_df<-read.table(paste0(args[3],".edge_df.txt"), h=T, stringsAsFactors=F, sep='\t')

if (is.null(tree$edge.length)){
	tree$edge.length<-edge_df$snp_count
}


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

getAncestors<-phytools:::getAncestors

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
    	#sometimes it can't add because there are too many samples there already

    	if (!best_nodes_table$sample[i] %in% newtree$tip.label){
    		final_node_in_new_tree<-getAncestors(newtree, node=final_node_in_new_tree)[1]
        	try(newtree<-bind.tip(tree=newtree, position = 0, edge.length=0.0001, where=final_node_in_new_tree,tip.label = best_nodes_table$sample[i]))
    		print(paste('check',best_nodes_table$sample[i], 'location', sep=' '))
    	}

    }
}

# newtree<-ladderize(newtree)


pdf(file=paste0("results_folder","/final_tree.pdf"), height=estimatePlotDimensions(tree)[[1]], width=estimatePlotDimensions(tree)[[2]])
plot(ladderize(newtree), show.tip.label=T,tip.color=ifelse(newtree$tip.label %in% best_nodes_table$sample, yes=2, no=1), cex=0.15)
dev.off()

write.tree(newtree,file=paste0('results_folder','/final_tree.nwk'))