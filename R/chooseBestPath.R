suppressWarnings(suppressPackageStartupMessages(library(scales)))
suppressWarnings(suppressPackageStartupMessages(library(phytools)))

getAncestors<-phytools:::getAncestors

tmpstr<-system('bash -l',input=c("shopt -s expand_aliases","type pathPhynder"), intern=T)

packpwd<-paste0(gsub('pathPhynder.R','',gsub('\'','',gsub('.*.Rscript ','',tmpstr))),'R')

args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=7) {
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
if (args[7]!='none'){
	call_hgs<-T
} else {
	call_hgs<-F
}


for (testfile in c(tree_file, sites_info_file, edge_df_file, calls_file)){
    if (!file.exists(testfile)) {
        stop(paste(testfile, "- does this file exist?"))
    }
}

edge_df<-read.table(paste0(args[2],".edge_df.txt"), h=T, stringsAsFactors=F, sep='\t')

sites_info<-read.table(sites_info_file, h=T)

tree<-read.tree(file=tree_file)
# tree<-ladderize(tree)
tree$tip.label<-make.names(tree$tip.label)


#IMPORTANT
#if the tree has no edge lengths, then add the SNP count at each branch as the edge.length 
if (is.null(tree$edge.length)){
	tree$edge.length<-edge_df$snp_count
	print('IMPORTANT - The tree has no edge lengths. Adding number of SNPs as edge length.')
}



dir.create(args[4], showWarnings = FALSE)



calls<-read.table(calls_file)

if(dim(calls)[1]==0){
	print('No genotypes were called. Make sure that the reference genome used and bam files have the same contig name!')
}

# colnames(sites_info)<-c('chr','marker','hg','pos','mutation','REF','ALT','branch','status')
colnames(calls)<-c('pos','REF','ALT','REFreads','ALTreads','geno')



source(paste0(packpwd,'/functions_pathPhynder.R'))

#assignAncientCallsToBranch
derived<-assignAncientCallsToBranch(calls, sites_info)$der
ancestral<-assignAncientCallsToBranch(calls, sites_info)$anc
all_counts<-rbind(derived,ancestral)




#makeBranchStatusTable
branch_counts <- makeBranchStatusTable(all_counts,edge_df)


#########
#choosing the best node
########

#make paths
paths<-makePaths(tree)


#get scores for each path

cat('Traversing tree...\n')

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

cat(paste0('Best node: ',best_node), '\n')



#get position in branch (0 if no conflict markers are seen)
position_in_branch<-estimatePositionInBranch(best_path)

best_path_counts<-getCountsforPath(best_path, branch_counts, "nodes")


cat('Writing best path report\n')

#make reports
best_path_report<-getCountsforPath(best_path_counts$Edge,branch_counts, "edges")

write.table(best_path_report, file=paste0(results_folder,"/",out_prefix,".best_path_report.txt"), sep='\t',row.names=F, col.names=T, quote=F)

count_all_paths<-makeCountsEveryPath(paths, branch_counts)

write.table(count_all_paths, file=paste0(results_folder,"/",out_prefix,".all_paths_report.txt"), sep='\t',row.names=F, col.names=T, quote=F)



####################



#makeSNPStatusOutput
snp_status<-makeSNPStatusOutput(all_counts)

write.table(all_counts, file=paste0(results_folder,"/",out_prefix,".all_counts"), sep='\t',row.names=F, col.names=T, quote=F)

table(all_counts$allele_status)

if (call_hgs){
	cat('Writing haplogroup reports\n')

	if (length(unique(all_counts$hg))==1){
		if (is.na(unique(all_counts$hg))){
			print("no hgs in data")
		}
	} else {
	#makeHaplogroupStatusOutput
	hg_status<-makeHaplogroupStatusOutput(all_counts)
	
	write.table(hg_status, file=paste0(results_folder,"/",out_prefix,".hg_in_tree_status.txt"), sep='\t',row.names=F, col.names=T, quote=F)
	
	
	# hg_status_derived<-hg_status[hg_status$derived_count>0,]
	# write.table(hg_status_derived, file=paste0(results_folder,"/",out_prefix,".hg_in_tree_status_derived_only.txt"), sep='\t',row.names=F, col.names=T, quote=F)
	# }
	

	
	df<-all_counts[all_counts$branch %in% best_path_report$Edge,]
	df$hg<-as.character(df$hg)
	df<-df[with(df, order(branch, hg)), ]
	hg_status_derived<-makeHaplogroupStatusOutput(df)
	write.table(hg_status_derived, file=paste0(results_folder,"/",out_prefix,".hg_in_tree_status_along_bestpath.txt"), sep='\t',row.names=F, col.names=T, quote=F)
	
	df$hg_status<-NA
	df$hg_status[which(as.character(df$hg_strand)==as.character(df$status) & df$allele_status=='Der')]<-"Der"	
	df$hg_status[which(as.character(df$hg_strand)==as.character(df$status) & df$allele_status=='Anc')]<-"Anc"	
	df$hg_status[which(as.character(df$hg_strand)!=as.character(df$status) & df$allele_status=='Anc')]<-"Der"	
	df$hg_status[which(as.character(df$hg_strand)!=as.character(df$status) & df$allele_status=='Der')]<-"Anc"	


	# df$allele_status<-NULL


	write.table(df, file=paste0(results_folder,"/",out_prefix,".snp_info_along_path"), sep='\t',row.names=F, col.names=T, quote=F)

	isogg_hgs_out<- paste0("intree_folder/",out_prefix,".hgs.txt")

	if(file.exists(isogg_hgs_out)){
		hg_isogg_status<-determineHG(isogg_hgs_out, best_path_report)
		write.table(hg_isogg_status, file=paste0(results_folder,"/",out_prefix,".isogg_hg_determination.txt"), sep='\t',row.names=F, col.names=T, quote=F)

	} else {
		print('for more detailed haplogroup information, run step 1 with with -G <known haplogroups list>, and then rerun step 2')
	}

	}
}

##################

pdf(file=paste0(results_folder,"/",out_prefix,".best_path.pdf"), height=estimatePlotDimensions(tree)[[1]], width=estimatePlotDimensions(tree)[[2]])

par(mar = rep(2, 4))
if (is.null(position_in_branch)){
	plotBestPathTree(tree,best_path_counts,branch_counts,path_scores, 0, best_node)
} else {
	plotBestPathTree(tree,best_path_counts,branch_counts,path_scores, position_in_branch, best_node)

}

dev.off()



edgeLen<-tree$edge.length[best_path_counts$Edge[best_path_counts$Node2==best_node]]
estimated_loc_at_branch=edgeLen*position_in_branch



best_node_info<-data.frame(sample=out_prefix,best_node=best_node,position_in_branch=position_in_branch,edgeLen=edgeLen, estimated_loc_at_branch=estimated_loc_at_branch )

write.table(best_node_info, file=paste0(results_folder,"/",out_prefix,".best_node_info.txt"), sep='\t',row.names=F, col.names=T, quote=F)


