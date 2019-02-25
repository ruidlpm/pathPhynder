suppressPackageStartupMessages(require(phytools))
suppressPackageStartupMessages(require(scales))
getAncestors<-phytools:::getAncestors


cat('\n\n',"chooseBestPath.R", '\n\n\n')

args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=6) {
    stop("  Arguments needed.\n
        \tusage
        \tRscript chooseBestPath.R <input_phylogeny.nwk> <prefix> intree.txt <results_folder> <maximumTolerance> <out_prefix>
        ", call.=FALSE)
}

# opt$input_tree,opt$prefix,paste0('intree_folder/',sample_name,'.intree.txt'), 'results_folder',  opt$maximumTolerance, sample_name ))

tree_file=args[1]
sites_info_file<-paste0(args[2],".sites.txt")
edge_df_file=paste0(args[2],".edge_df.txt")
results_folder=args[4]
maximumTolerance=args[5]
out_prefix=args[6]
calls_file<-paste0("intree_folder/",args[6], '.intree.txt')

for (testfile in c(tree_file, sites_info_file, edge_df_file, calls_file)){
    if (!file.exists(testfile)) {
        stop(paste(testfile, "- does this file exist?"))
    }
}
















edge_df<-read.table(paste0(args[2],".edge_df.txt"), h=T, stringsAsFactors=F, sep='\t')

sites_info<-read.table(sites_info_file)

tree<-read.tree(file=tree_file)
tree<-ladderize(tree)
tree$tip.label<-make.names(tree$tip.label)

dir.create(args[4], showWarnings = FALSE)



calls<-read.table(calls_file)
colnames(sites_info)<-c('chr','marker','hg','pos','mutation','REF','ALT','branch','status')
colnames(calls)<-c('pos','REF','ALT','REFreads','ALTreads','geno')




#########
#assigning ancient SNPs to branches
########

source('~/in_development/pathPhynder/R/workInProgress.R')

#assignAncientCallsToBranch
derived<-assignAncientCallsToBranch(calls, sites_info)$der
ancestral<-assignAncientCallsToBranch(calls, sites_info)$anc
all_counts<-rbind(derived,ancestral)





# print(all_counts)

#makeSNPStatusOutput
snp_status<-makeSNPStatusOutput(all_counts)

table(all_counts$allele_status)

#makeHaplogroupStatusOutput
# hg_status<-makeHaplogroupStatusOutput(all_counts)
# write.table(hg_status, file=paste0(results_folder,"/",out_prefix,".hg_in_tree_status.txt"), sep='\t',row.names=F, col.names=T, quote=F)

#makeBranchStatusTable
branch_counts <- makeBranchStatusTable(all_counts,edge_df)



#########
#choosing the best node
########

#make paths
paths<-makePaths(tree)

print(maximumTolerance)


#get scores for each path
path_scores<-traversePaths(paths,branch_counts, maximumTolerance)

if (sum(path_scores$total_derived)==0){

	pdf(file=paste0(results_folder,"/",out_prefix,".branch_counts_no_path.pdf"), height=estimatePlotDimensions(tree)[[1]], width=estimatePlotDimensions(tree)[[2]])
	plotAncDerSNPTree(branch_counts)
	dev.off()

	count_all_paths<-makeCountsEveryPath(paths, branch_counts)
	write.table(count_all_paths, file=paste0(results_folder,"/",out_prefix,".all_paths_report.txt"), sep='\t',row.names=F, col.names=T, quote=F)

	stop("No Derived SNPs observed, can't choose path.")
}


#chose best path
best_path<-chooseBestPath(path_scores,branch_counts)

best_node<-best_path[length(best_path)]



#get position in branch (0 if no conflict markers are seen)
position_in_branch<-estimatePositionInBranch(best_path)

best_path_counts<-getCountsforPath(best_path, branch_counts, "nodes")



#make reports
best_path_report<-getCountsforPath(best_path_counts$Edge,branch_counts, "edges")

write.table(best_path_report, file=paste0(results_folder,"/",out_prefix,".best_path_report.txt"), sep='\t',row.names=F, col.names=T, quote=F)

count_all_paths<-makeCountsEveryPath(paths, branch_counts)

write.table(count_all_paths, file=paste0(results_folder,"/",out_prefix,".all_paths_report.txt"), sep='\t',row.names=F, col.names=T, quote=F)



pdf(file=paste0(results_folder,"/",out_prefix,".best_path.pdf"), height=estimatePlotDimensions(tree)[[1]], width=estimatePlotDimensions(tree)[[2]])

	par(mar = rep(2, 4))
if (is.null(position_in_branch)){
	plotBestPathTree(tree,best_path_counts,branch_counts,path_scores, 0, best_node)
} else {
	plotBestPathTree(tree,best_path_counts,branch_counts,path_scores, position_in_branch, best_node)

}
dev.off()



































# # getReadCountForPath<-function(){
# # #return(height, width)
# # }



# # AddAncientToTree<-function(){
# # #return(height, width)
# # }





# plotAncDerSNPTree(branch_counts)
# edgelabels(edge=unique(path_scores$stopped_edges),pch="X", col=1,cex=1)


# path_scores<-traversePaths(paths,branch_counts,20)

# plotAncDerSNPTree(branch_counts)
# edgelabels(edge=unique(path_scores$stopped_edges),pch="X", col=1,cex=1)











# plotAncDerSNPTree(branch_counts)
# edgelabels(edge=unique(path_scores$stopped_edges),pch="X", col=1,cex=1)
# edgelabels(edge=unique(best_path_counts$Edge),pch="B", cex=1, col="green")






# getCountsforPath(paths[[67]],branch_counts, "nodes")



# branch_counts[branch_counts$Edge==130,]$conflict<-100







# pdf(file="test.pdf", height=15, width=10)
# plotAncDerSNPTree(branch_counts)
# edgelabels(edge=unique(path_scores$stopped_edges),pch="X", col=1,cex=1)
# dev.off()
# system("open test.pdf")


# # edgelabels(edge=conflict_branch_counts$Edge,pch=20, col=alpha("red", 0.5), cex=log((conflict_branch_counts$conflict)+1)/2)


# edgelabels(edge=unique(stop_edges))

# queuePaths

























