require(phytools)
require(scales)

# usage:
# Rscript decide.R <input_phylogeny.nwk> <input_prefix> <results_folder> 
# Description: reads tree, aDNA intree files in intree_folder, outputs result to results_folder.

args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=4) {
    stop("  Arguments needed.\n
        \tusage
        \tRscript decide.R <input_phylogeny.nwk> <results_folder> <Maximum Tolarance (int)> <prefix>
        ", call.=FALSE)
} else {
    cat("   Command used:",'\n\n')
    cat(paste("decide.R", args[1],args[2],args[3],args[4], '\n\n'))
}


edge_df<-read.table(paste0(args[4],".edge_df.txt"), h=T, stringsAsFactors=F, sep='\t')

lens<-NULL
for (i in edge_df$Edge){
    tmp<-edge_df[edge_df$Edge==i,]
    lens[i]<-(length(unique(unlist(strsplit(tmp$positions[1], '\\;')))))
}

nums<-cbind(edge_df, number_of_positions=lens)


#read numbers of snps
# nums<-read.table("numbers.txt", h=T, stringsAsFactors=F)

#read tree
tree<-read.tree(file=args[1])
tree<-ladderize(tree)

#replace - with . to prevent downstream problems withs sample ID
tree$tip.label<-make.names(tree$tip.label)

#read anchgs table
try(anchgs<-(read.table("anchgs.txt", sep='\t')))

br<-cbind(data.frame(tree$edge), c(1:length(tree$edge[,1])))
colnames(br)<-c("pos1", "pos2", "branches")

#read table with derived and ancestral marker count at each branch
br_sample_tables<-readRDS(file="results_folder/allele_count_list.RData")

samps<-names(br_sample_tables)

try(hgs<-try(read.table("hgs.txt")))
tmptree<-tree
try(tmptree$tip.label<-paste0(hgs$V4[match(make.names(tree$tip.label), make.names(hgs$V1))],"___",tree$tip.label))

#make all possible paths
getAncestors<-phytools:::getAncestors

#make all possible paths to transverse, including the tip
paths<-list()
for(i in 1:length(tree$tip.label)){
    ancestors<-getAncestors(tree, i, type='all')
    paths[i]<-data.frame(c(rev(ancestors), i))
}
# paths
print(paste0("number of paths = ",length(paths)))
print(paste0("number of tips = ",length(tree$tip.label)))


max_tolerance<-as.numeric(args[3])


samples<-names(br_sample_tables)

best_nodes_table<-NULL




#for each sample
for (samp in 1:length(br_sample_tables)){

    test_not_support<-NULL
    sample_name<-samples[samp]
    print(sample_name)

    countdata<-br_sample_tables[[samp]]

    pos_score<-list()
    neg_score<-list()

    stopped_edges<-NULL
    stopped_nodes<-NULL
    for (path_index in 1:length(paths)){
        path<-paths[[path_index]]

        path_pos_score<-0
        path_neg_score<-0
        stopped_node<-NULL
        stopped_edge<-NULL

        #for each path
        for (pos1_in_path in 1:length(path)){
            # if the position in the path has been stopped before, do not continue
            # break
            if (sum(path[pos1_in_path] %in% stopped_nodes)==0 ){
                entry<-countdata[which(countdata$Node1==path[pos1_in_path] & countdata$Node2==path[pos1_in_path+1]),]

                if(dim(entry)[1]==1){
                    supporting<-NULL
                    notsupporting<-NULL
                    supporting = entry$support
                    notsupporting = entry$notsupport
                    
                    #if the number of nonsupprting alleles is greater or equal to the max tolerance param
                    if (notsupporting>=max_tolerance){
                        #add node to stopped node and stopped edges and switch to a different path
                        test_not_support<-rbind(test_not_support,entry)
                        stopped_node<-path[pos1_in_path + 1]
                        # stopped_node<-entry$Node2[which(entry$Node1==path[pos1_in_path] & entry$Node2==path[pos1_in_path+1])]
                        stopped_nodes<-c(stopped_nodes,stopped_node )
                        stopped_edge<-entry$Edge[which(entry$Node1==path[pos1_in_path] & entry$Node2==path[pos1_in_path+1])]
                        # stopped_edge<-entry$Edge[which(entry$Node1==path[pos1_in_path+1])]
                        stopped_edges<-c(stopped_edges,stopped_edge )
                        break
                    #if the number of nonsupprting alleles is NOT greater or equal to the max tolerance param, add score to path
                    } else {
                        path_pos_score<-c(path_pos_score,supporting)
                        path_neg_score<-c(path_neg_score,notsupporting)
                     }
                }
            }
        }
        #if there is at least one marker supporting branch assignment in any branch, add the counts to the path
        if(length(path_pos_score)>0){
            pos_score[path_index]<-data.frame(path_pos_score)
            neg_score[path_index]<-data.frame(path_neg_score)
        }
    }



#decide best path - the one containing higher number of derived alleles.

sums<-sapply(pos_score, sum)


if (length(which(sums==max(sums)))>1){
    best_paths<-which(sums==max(sums))
    tmp_paths<-NULL
    for (i in 1:length(best_paths)){
        tmp_paths<-c(paths[[which(sums==max(sums))[i]]])
    }
    best_path<-unique(tmp_paths)
    # print(which(sums==max(sums)))
    # break
} else if (length(which(sums==max(sums)))==1) {
    best_path<-paths[[which(sums==max(sums))]]
} else if (length(which(sums==max(sums)))==0){
    cat("insufficient derived alleles")
}


counts_for_best_path<-NULL
counts_for_best_path<-countdata[match(best_path, countdata$Node2),]


counts_for_best_path$number_of_positions<-nums$number_of_positions[match(counts_for_best_path$Edge, nums$Edge)]



counts_for_best_path<-counts_for_best_path[!is.na(counts_for_best_path$Edge),]
counts_for_best_path_with_missing<-counts_for_best_path
colnames(counts_for_best_path_with_missing)[5]<-"conflict"



counts_for_best_path_with_missing<-counts_for_best_path_with_missing[c('Edge','Node1','Node2','support','conflict','number_of_positions','hg')]
print(counts_for_best_path_with_missing)
# counts_for_best_path_with_missing<-counts_for_best_path_with_missing[c('Edge','Node1','Node2','support','conflict','hg')]

write.table(counts_for_best_path_with_missing, file=paste0(args[2],'/',sample_name,'.best_path_all_branches.txt'),quote=F, row.names=F, sep="\t")

counts_for_best_path<-counts_for_best_path[!(counts_for_best_path$support==0 &  counts_for_best_path$notsupport==0),]



counts_for_best_path_toplot<-counts_for_best_path[!(counts_for_best_path$support==0 &  counts_for_best_path$notsupport>0),]


write.table(counts_for_best_path, file=paste0(args[2],'/',sample_name,'.txt'),quote=F, row.names=F, sep="\t")

print(paste(gsub(".intree.txt","",sample_name),anchgs$V2[match(gsub(".intree.txt","",sample_name), make.names(anchgs$V1))],anchgs$V4[match(gsub(".intree.txt","",sample_name), make.names(anchgs$V1))], sep="___"))

if (sum(counts_for_best_path$support>0)==0){
    print("insufficient data")
    next
} else {
    print("OK")
}


best_edge<-NULL
best_node1<-NULL
best_node2<-NULL
best_edge<-counts_for_best_path$Edge[counts_for_best_path$support>0][sum(counts_for_best_path$support>0)]
best_node1<-counts_for_best_path$Node1[counts_for_best_path$support>0][sum(counts_for_best_path$support>0)]
best_node2<-counts_for_best_path$Node2[counts_for_best_path$support>0][sum(counts_for_best_path$support>0)]

der_at_best<-counts_for_best_path$support[counts_for_best_path$Edge==best_edge]
anc_at_best<-counts_for_best_path$notsupport[counts_for_best_path$Edge==best_edge]

pos_at_best_edge<-NULL
if (anc_at_best==0){
    pos_at_best_edge<-0
    total_len<-0
} else if (anc_at_best>0){
    total_len<-tmptree$edge.length[best_edge]
    total_count<-der_at_best+anc_at_best
    fraction_at_best_edge=der_at_best/total_count
    pos_at_best_edge=total_len*fraction_at_best_edge
    print(total_len)
    print(fraction_at_best_edge)
    print(pos_at_best_edge)
}



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



# pdf(file=paste0(args[2],'/',sample_name,'.decision.pdf'), height=10, width=7)
pdf(file=paste0(args[2],'/',sample_name,'.decision.pdf'), height=16, width=8)
plot((tmptree), col='darkgrey',cex=0.15, edge.col=ifelse(countdata$Edge %in% unique(stopped_edges), yes=2, no='lightgrey'), show.tip.label = T, edge.width = 1, tip.color = 0)
par(new=T)
plot((tmptree), cex=0.15, edge.col=ifelse(countdata$Edge %in% unlist(counts_for_best_path_toplot$Edge), yes=3, no=0), show.tip.label = T, edge.width = 1, tip.color = "grey3")
edgelabels(edge=countdata$Edge[countdata$notsupport>0],pch=20, col=alpha("red", 0.5), cex=log((countdata$notsupport[countdata$notsupport>0])+1)/2)
edgelabels(edge=countdata$Edge[countdata$support>0],pch=20, col=alpha("darkgreen", 0.7),cex=log((countdata$support[countdata$support>0])+1)/2)
edgelabels(frame="none",edge=countdata$Edge[countdata$support>0], text=countdata$support[countdata$support>0], cex=0.3)
edgelabels(frame="none",edge=countdata$Edge[countdata$notsupport>0], text=countdata$notsupport[countdata$notsupport>0], cex=0.3)
# nodelabels(node=stopped_nodes, pch=20, col=1)
if(pos_at_best_edge==0){
    try(x <- getphylo_x(tmptree, best_node2))
    try(y <- getphylo_y(tmptree, best_node2))
    points(x,y,  col="black",bg=alpha("yellow", 0.3), pch=21, cex=1) 
} else {
    try(x <- getphylo_x(tmptree, best_node2))
    try(y <- getphylo_y(tmptree, best_node2))
    points(x-(total_len-pos_at_best_edge),y, col="black", bg=alpha("yellow", 0.3), pch=21, cex=1)
}

# try(mtext(paste(anchgs$V1[grep(unlist(strsplit(sample_name, '\\.'))[1], make.names(anchgs$V1))],anchgs$V4[grep(unlist(strsplit(sample_name, '\\.'))[1], make.names(anchgs$V1))], sep="___")))
# try(print(anchgs[grep(unlist(strsplit(sample_name, '\\.'))[1], make.names(anchgs$V1)),]))
mtext(paste(gsub(".intree.txt","",sample_name),anchgs$V2[match(gsub(".intree.txt","",sample_name), make.names(anchgs$V1))],anchgs$V4[match(gsub(".intree.txt","",sample_name), make.names(anchgs$V1))], sep="___"))

# nodelabels()
dev.off()



best_nodes_table$sample_name<-as.character(best_nodes_table$sample_name)
best_nodes_table<-rbind(best_nodes_table,data.frame(sample_name,best_node1, best_node2,pos_at_best_edge, total_len))

closeAllConnections()
}


write.table(file=paste0(args[2],'/best_nodes_table.txt'),best_nodes_table, quote=F, row.names=F, sep='\t')

pdf(file=paste0(args[2],'/best_node_location.pdf'), height=20, width=14)

plot(tmptree, cex=0.3)

for (i in 1:length(best_nodes_table$sample_name)){

if(best_nodes_table$pos_at_best_edge[i]==0){
    try(x <- getphylo_x(tmptree, best_nodes_table$best_node2[i]))
    try(y <- getphylo_y(tmptree, best_nodes_table$best_node2[i]))
    points(x,y,  col="black",bg=alpha("yellow", 0.3), pch=21, cex=1) 
} else {
    try(x <- getphylo_x(tmptree, best_nodes_table$best_node2[i]))
    try(y <- getphylo_y(tmptree, best_nodes_table$best_node2[i]))
    points(x-(best_nodes_table$total_len[i]-best_nodes_table$pos_at_best_edge[i]),y, col="black", bg=alpha("yellow", 0.3), pch=21, cex=1)
}


}
dev.off()



# require(phytools)

# best_nodes_table<-read.table("results_folder/best_nodes_table.txt", h=T, stringsAsFactors=F)
#  tmptree<-read.tree('/Users/rm890/testing/na_algo/revisiting_code/tree_sole_only.nwk')

best_nodes_table<-best_nodes_table[order(best_nodes_table$pos_at_best_edge-best_nodes_table$total_len),]

newtree<-tmptree
newtree$edge.length<-tmptree$edge.length

dup_entr<-best_nodes_table[c('best_node1' ,'best_node2' ,'pos_at_best_edge'  ,'total_len')][duplicated.data.frame(best_nodes_table[c('best_node1' ,'best_node2' ,'pos_at_best_edge'  ,'total_len')]),]

dup_table<-merge(best_nodes_table, dup_entr)
dup_table<-dup_table[order(dup_table$pos_at_best_edge-dup_table$total_len),]

uniq_entr<-best_nodes_table[!best_nodes_table$sample_name %in%dup_table$sample_name ,]
uniq_entr<-rbind(dup_table[dup_table$pos_at_best_edge==0,],uniq_entr)
uniq_entr<-uniq_entr[order(uniq_entr$pos_at_best_edge-uniq_entr$total_len),]
dup_table<-dup_table[dup_table$pos_at_best_edge!=0,]

dup_table$pos_at_best_edge[which(duplicated(dup_table$pos_at_best_edge))]<-dup_table$pos_at_best_edge[which(duplicated(dup_table$pos_at_best_edge))]-0.0000000001
uniq_entr<-rbind(dup_table,uniq_entr)
uniq_entr<-unique(uniq_entr[order(uniq_entr$pos_at_best_edge-uniq_entr$total_len),])



# try(best_nodes_table$hg_label<-paste(gsub(".intree.txt","",best_nodes_table$sample_name),anchgs$V4[match(gsub(".intree.txt","",best_nodes_table$sample_name), anchgs$V1)], sep="___"))

for (i in 1:length(uniq_entr$best_node2)){
    print(i)
    # Get descendants of final_node in tmptree
    desc<-(getDescendants(tmptree, uniq_entr$best_node2[i]))
    desc_tip_names<-tmptree$tip.label[tmptree$tip.label %in% tmptree$tip.label[desc]]

    if(length(desc_tip_names)==1){
        #if just 1 descendant, then final node is tip
        final_node_in_new_tree<-which(newtree$tip.label==desc_tip_names)
    } else if(length(desc_tip_names)>1){
        #if more than 1 descendant, then need to find MRCA of all descendants
        final_node_in_new_tree<-(findMRCA(newtree, desc_tip_names))
    }
    if (uniq_entr$pos_at_best_edge[i]==0){
        newtree<-bind.tip(tree=newtree, position = 0, edge.length=0.001, where=final_node_in_new_tree,tip.label = paste(gsub(".intree.txt","",uniq_entr$sample_name[i]),anchgs$V4[match(gsub(".intree.txt","",uniq_entr$sample_name[i]), anchgs$V1)], sep="___"))
    } else if (uniq_entr$pos_at_best_edge[i]>0){
        try(newtree<-bind.tip(tree=newtree, position =(uniq_entr$total_len[i]-uniq_entr$pos_at_best_edge[i]), edge.length=abs(uniq_entr$pos_at_best_edge[i]-uniq_entr$total_len[i]), where=final_node_in_new_tree, tip.label = paste(gsub(".intree.txt","",uniq_entr$sample_name[i]),anchgs$V2[match(gsub(".intree.txt","",uniq_entr$sample_name[i]), make.names(anchgs$V1))],anchgs$V4[match(gsub(".intree.txt","",uniq_entr$sample_name[i]), make.names(anchgs$V1))], sep="___")))
    }
}


new_tree_no_labels<-tree

for (i in 1:length(uniq_entr$best_node2)){
    print(i)
    # Get descendants of final_node in tmptree
    desc<-(getDescendants(tree, uniq_entr$best_node2[i]))
    desc_tip_names<-tree$tip.label[tree$tip.label %in% tree$tip.label[desc]]

    if(length(desc_tip_names)==1){
        #if just 1 descendant, then final node is tip
        final_node_in_new_tree<-which(new_tree_no_labels$tip.label==desc_tip_names)
    } else if(length(desc_tip_names)>1){
        #if more than 1 descendant, then need to find MRCA of all descendants
        final_node_in_new_tree<-(findMRCA(new_tree_no_labels, desc_tip_names))
    }
    if (uniq_entr$pos_at_best_edge[i]==0){
        new_tree_no_labels<-bind.tip(tree=new_tree_no_labels, position = 0, edge.length=0.001, where=final_node_in_new_tree,tip.label = gsub(".intree.txt","",uniq_entr$sample_name[i]))
    } else if (uniq_entr$pos_at_best_edge[i]>0){
        try(new_tree_no_labels<-bind.tip(tree=new_tree_no_labels, position =(uniq_entr$total_len[i]-uniq_entr$pos_at_best_edge[i]), edge.length=abs(uniq_entr$pos_at_best_edge[i]-uniq_entr$total_len[i]), where=final_node_in_new_tree, tip.label = gsub(".intree.txt","",uniq_entr$sample_name[i])))
    }
}



write.tree(new_tree_no_labels,file=paste0(args[2],'/added_anc_best_node_location_nolabels.nwk'))





# plot(ladderize(newtree), cex=0.5, tip.color = ifelse(newtree$tip.label %in% substring(best_nodes_table$sample_name,1,9), yes=2, no=1))


# get_subtrees<-function(best_nodes_tables){
#     #make subtrees with the samples that 1) belong to the same node; 2) have the same pos_at_best_edge
#     #split by node
#     node_split<-(split(best_nodes_tables,best_nodes_tables$best_node2))
    
#     subtree_list<-list()
#     counter=0
#     for (split_num in 1:length(node_split)){
#        #get split by each node2
#        tmp<-node_split[[split_num]]
#        #split within node by pos_at_best_edge
#         pos_split<-(split(tmp, tmp$pos_at_best_edge))
#         if (length(pos_split)>1){
#             dups<-which(sapply(pos_split, dim)[1,]>1)
#             if (length(dups)>1){
#                 for (dup_num in 1:length(dups)){
#                     print(pos_split[dups][[dup_num]])
#                     subtree_data<-pos_split[dups][[dup_num]]
#                     subtree_len<-dim(subtree_data)[1]
#                     subtree<-pbtree(n=subtree_len)
#                     subtree$tip.label<-subtree_data$sample_name
#                     subtree$edge.length<-subtree_data$total_len-subtree_data$pos_at_best_edge
#                     counter=counter+1
#                     subtree_list[[counter]]<-subtree
#                 }
#             }
#         }
#     }
#     return(subtree_list)
# }
# 
# all_subtrees<-get_subtrees(dup_table)







write.tree(newtree,file=paste0(args[2],'/added_anc_best_node_location.nwk'))

pdf(file=paste0(args[2],'/added_anc_best_node_location.pdf'), height=100, width=20)
plot(ladderize(newtree), cex=0.5, tip.color = ifelse(newtree$tip.label %in% paste(gsub(".intree.txt","",uniq_entr$sample_name),anchgs$V4[match(gsub(".intree.txt","",uniq_entr$sample_name), anchgs$V1)], sep="___"), yes=2, no=1))
dev.off()







# pdf(file="test.pdf")
# plot(tmptree, cex=0.2, edge.col=ifelse(countdata$Edge %in% unique(unlist(visited)), yes=3, no="lightgrey"), show.tip.label = T, edge.width = 1, tip.color = 1)
# edgelabels(frame="none",edge=countdata$Edge[countdata$support>0], text=countdata$support[countdata$support>0], col=3)
# edgelabels(frame="none",edge=countdata$Edge[countdata$notsupport>0], text=countdata$notsupport[countdata$notsupport>0], col=2)
# nodelabels(node=best_path)
# edgelabels(edge=br$branches[br$branches %in% unlist(stopped_path)])
# dev.off()



# plot(tmptree, cex=0.2, edge.col=ifelse(br$branches %in% unlist(dd$Edge), yes=3, no="lightgrey"), show.tip.label = T, edge.width = 1, tip.color = 1)

# par(new=T)
# plot(tmptree, cex=0.2, edge.col=ifelse(br$branches %in% stopped_path, yes=2, no=0), show.tip.label = T, edge.width = 1, tip.color = 0)
# edgelabels(edge=countdata$Edge[countdata$support>0], text=countdata$support[countdata$support>0])
# edgelabels(frame="none",edge=countdata$Edge[countdata$notsupport>0], text=countdata$notsupport[countdata$notsupport>0], col=2)
# dev.off()



# stopped_path=NULL
# continued_path=NULL
# continued_paths=list()
# temp=NULL
# score=list()
# counter=0
# for (path in paths){
#     counter=counter+1
#     score_total=0
#     continued_path=NULL
#     if (sum(path %in% stopped_path)==0){
#     for (pos1_node in 1:length(path)){
#         if (pos1_node<length(path)){
#             entry<-(countdata[which(countdata$Node1==path[pos1_node] & countdata$Node2==path[pos1_node+1]),])
#             supporting = entry$support
#             notsupporting = entry$notsupport
            
#             branch = entry$Edge
# #                if (notsupporting<2 ){
#                if (supporting>0 ){
#                    continued_path=cbind(continued_path, branch)
# #                                   if (notsupporting<10){
# #                    continued_path=cbind(continued_path, branch)
           
#                    score_total=score_total + supporting-notsupporting
#                 } else if (notsupporting >=30000000){
#                    stopped_path=cbind(stopped_path, branch)
#                    score_total=score_total + supporting-notsupporting
#                    continued_paths[[counter]]=continued_path
#                    score[[counter]]<-score_total
#                    break
                   
# #                 } else if (supporting >= 1 && notsupporting ==2){
#                 } else if (supporting >= 1 && notsupporting >=30000000){
#                     continued_path=cbind(continued_path, branch)
#                    score_total=score_total + supporting-notsupporting
#                 } else if (supporting == 0 && notsupporting ==30000000){
#                     stopped_path=cbind(stopped_path, branch)
#                    continued_paths[[counter]]=continued_path
#                     score_total=score_total + supporting-notsupporting
#                    score[[counter]]<-score_total
#                     break
#                 } else {
#                     temp=cbind(temp, branch)
#                 }
#             }
            
#         }
#         continued_paths[[counter]]=continued_path
#         score[[counter]]<-score_total
#     }
# }

# print(samp)
# if (length(continued_paths)==0){
# 	print("insufficient data")
# } else {

# test<-NULL
# for (i in 1:length(continued_paths)){
#     if (length(continued_paths[[i]])==0){
#     print('insufficient data')
# } else {
#     test[i]<-(sum(countdata$support[match(unlist(continued_paths[[i]]), countdata$branches)]))
# }
#     }

    
#         descendants<-tmptree$tip.label[which(test==max(test))]
#     # descendants
#         if (length(descendants)>1){
#             best_node <- findMRCA(tmptree,descendants)
#             best_nodes_per_sample[[samp]]<-best_node
#         } else if (length(descendants)==1){  
#             ###################################
#             # attention tip is not a node
#             best_node<-which(descendants==tmptree$tip.label)
#             best_nodes_per_sample[[samp]]<-best_node
#         }

# print(c(samp, best_node))
# print(descendants)
#     print(sum(countdata$support>0))
    
# if (sum(countdata$support>0)==0){
#     print('insufficient data')
# } else {
    
    
    
# pdf(file=paste0(args[2], '/',samp,'.decision.pdf'), height=15, width=10)


# plot(tmptree, cex=0.2, edge.col=ifelse(br$branches %in% c(unlist(continued_paths),stopped_path), yes=0, no='lightgrey'), show.tip.label = T, tip.color = 0)
# par(new=T)
# plot(tmptree, cex=0.2, edge.col=ifelse(br$branches %in% unlist(continued_paths), yes=3, no=0), show.tip.label = T, edge.width = 1, tip.color = 0)
# par(new=T)
# plot(tmptree, cex=0.2, edge.col=ifelse(br$branches %in% stopped_path, yes=2, no=0), show.tip.label = T)
# par(new=T)
# edgelabels(edge=countdata$branches[countdata$support>0],pch=20, col=alpha("darkgreen", 0.7),cex=log((countdata$support[countdata$support>0])+1)/2)
# edgelabels(edge=countdata$branches[countdata$notsupport>0],pch=20, col=alpha("red", 0.5), cex=log((countdata$notsupport[countdata$notsupport>0])+1)/2)
# edgelabels(edge=countdata$branches[countdata$support>0],text=countdata$support[countdata$support>0], frame = "none", cex=0.3)
# edgelabels(edge=countdata$branches[countdata$notsupport>0],text=countdata$notsupport[countdata$notsupport>0], frame = "none", cex=0.3)
# nodelabels(node=best_node, pch=20, col=alpha("yellow", 1), cex=0.5)
    
# try(mtext(text = paste(anchgs$V1[grep(gsub('^X','',samp), anchgs$V1)], anchgs$V4[grep(gsub('^X','',samp), anchgs$V1)], sep='___')))
#     print(paste(anchgs$V1[grep(gsub('^X','',samp), anchgs$V1)], anchgs$V4[grep(gsub('^X','',samp), anchgs$V1)], sep='___'))
    
#     dev.off()

#     }
# }

# }
# saveRDS(file="best_nodes_per_sample.Rdata", best_nodes_per_sample)





# estimate_where_in_branch_to_bind<-function(x){
#     fraction_supporting<- x$support/(x$support+x$notsupport)
#     place_in_branch<-tmptree$edge.length[mismatch_info$branches]*(1-fraction_supporting)
#     return(place_in_branch)
# }


# best_nodes_table<-matrix(ncol = 3, nrow=length(best_nodes_per_sample))

# for (i in 1:length(best_nodes_per_sample)){
#     print(names(best_nodes_per_sample[i]))
#     line_of_mismatch<-which(br_sample_tables[[i]]$support[br_sample_tables[[i]]$pos2==best_nodes_per_sample[[i]]]>0)
#     mismatch_info<-NULL
#     mismatch_info<-br_sample_tables[[i]][br_sample_tables[[i]]$pos2==best_nodes_per_sample[[i]],][line_of_mismatch,]
#     print(mismatch_info)
#     where_in_branch_to_bind<-0
#     if (length(mismatch_info$support)==0){
#         new_node_to_bind<-as.numeric(as.character(best_nodes_per_sample[[i]]))
#         where_in_branch_to_bind<-0
#     } else if (mismatch_info$support>0){
#     new_node_to_bind<-mismatch_info$pos2
#     where_in_branch_to_bind<-estimate_where_in_branch_to_bind(mismatch_info)
#     } else if (exists(mismatch_info)==F) {
#         print('prob')
#     }
#     best_nodes_table[i,][1]<-names(best_nodes_per_sample[i])
#     best_nodes_table[i,][2]<-new_node_to_bind
#     best_nodes_table[i,][3]<-where_in_branch_to_bind
# }

# best_nodes_table<-data.frame(best_nodes_table)
# colnames(best_nodes_table)<-c('sampleID','final_node','pos')
# best_nodes_table$pos<-as.numeric(as.character(best_nodes_table$pos))
# best_nodes_table$final_node<-as.numeric(as.character(best_nodes_table$final_node))
# best_nodes_table$sampleID<-as.character(best_nodes_table$sampleID)

# best_nodes_table<-best_nodes_table[rev(order(best_nodes_table$pos)),]
# best_nodes_table




# get_supporting_snps<-function(x){
#     br_sample_tables[[x]][br_sample_tables[[x]]$support>0 & br_sample_tables[[x]]$notsupport<5,]
# }

# for (i in names(br_sample_tables)){
#     print(i)
#     print(get_supporting_snps(i))
# #     tmp<-get_supporting_snps(i)
# #     print(tmp[which(length(tmp$hgs_at_branch) > 0),])
# }





# getphylo_x <- function(tree, node) {
#     if(is.character(node)) {
#         node <- which(c(tree$tip.label, tree$node.label)==node)
#     }
#     pi <- tree$edge[tree$edge[,2]==node, 1]
#     if (length(pi)) {
#         ei<-which(tree$edge[,1]==pi & tree$edge[,2]==node)
#         tree$edge.length[ei] + Recall(tree, pi)
#     } else {
#         if(!is.null(tree$root.edge)) {
#             tree$root.edge
#         } else {
#             0
#         }
#     }
# }

# getphylo_y <- function(tree, node) {
#     if(is.character(node)) {
#         node <- which(c(tree$tip.label, tree$node.label)==node)
#     }
#     ci <- tree$edge[tree$edge[,1]==node, 2]
#     if (length(ci)==2) {
#         mean(c(Recall(tree, ci[1]), Recall(tree, ci[2])))
#     } else if (length(ci)==0) {
#         Ntip <- length(tree$tip.label)
#         which(tree$edge[tree$edge[, 2] <= Ntip, 2] == node)
#     } else {
#         stop(paste("error", length(ci)))
#     }
# }




# # fu_samps<-c('Bockstein','Brillenhohle','Burkhardtshohle','Chaudardes1','Cioclovina1','GoyetQ-2','GoyetQ116-1','GoyetQ53-1','HohleFels49','HohleFels79','Iboussieres39','Karelia','Kostenki12','Kostenki14','KremsWA3','Ofnet','Paglicci108','Paglicci133','Pavlov1','Vestonice13','Vestonice14','Vestonice15','Vestonice16','Vestonice43','Villabruna')

# pdf(file="tmp.ycap1.pdf", height=20, width=10)

# plot(tmptree, cex=0.3)
# # nodelabels(node=best_nodes_table$final_node, col=3, cex=3, pch=20,adj = c(best_nodes_table$pos,0))
# # nodelabels(node=best_nodes_table$final_node,text=best_nodes_table$sampleID, col=1, cex=1, frame = 'none', adj = c(0,best_nodes_table$pos))

# for (i in 1:length(best_nodes_table$final_node)){
# try(x <- getphylo_x(tree, best_nodes_table$final_node[i]))
# try(y <- getphylo_y(tree, best_nodes_table$final_node[i]))
# # points(x-best_nodes_table$pos[i],y, bg=3, pch=21, cex=1)

# points(x-best_nodes_table$pos[i],y, bg=alpha(ifelse(best_nodes_table$sampleID[i] %in% samps, yes=2, no='white'), 0.3), pch=21, cex=1)
# # points(x-best_nodes_table$pos[i],y, bg=alpha(1, 0.3), pch=21, cex=1)
# # text(x-best_nodes_table$pos[i],y, col=alpha(best_nodes_table$cols[i], 0.3), pch=20, cex=0.4, paste(best_nodes_table$sampleID[i],anchgs$V2[match(best_nodes_table$sampleID[i],anchgs$V1)], sep='___'), pos=1)

# #     by time
#     # points(x-best_nodes_table$pos[i],y, bg=anchgs$col[match(best_nodes_table$sampleID[i], anchgs$V1)], pch=21, cex=1)
# # text(x-best_nodes_table$pos[i],y, col=alpha(best_nodes_table$cols[i], 0.3), pch=20, cex=0.4, paste(best_nodes_table$sampleID[i],anchgs$V2[match(best_nodes_table$sampleID[i],anchgs$V1)], sep='___'), pos=1)

# }
# dev.off()


# dim(best_nodes_table)

# saveRDS(file="best_nodes_table.Rdata", best_nodes_table)



# sample_labels<-NULL
# newtree<-tmptree
# newtree$edge.length<-tmptree$edge.length
# for (i in 1:length(best_nodes_table$final_node)){
# # if (!(best_nodes_table$sampleID[i] %in% c('RISE664', 'DA345', 'RISE662', 'RISE718'))){
# #     next
# # } else {
    
#     # Get descendants of final_node in tmptree
#     desc<-(getDescendants(tmptree, best_nodes_table$final_node[i]))
#     desc_tip_names<-tmptree$tip.label[tmptree$tip.label %in% tmptree$tip.label[desc]]

#     # Search same descendants in new_tree and get MRCA node of those descendants
#     final_node_in_new_tree<-(findMRCA(newtree, desc_tip_names))
# #     print(final_node_in_new_tree)

#     #if no MRCA is found it probably means that the final node is located at a terminal branch (a tip)
#     if (length(final_node_in_new_tree)==0){
#         #search for descendants of best node in current tree
#         final_node_in_new_tree<-which(newtree$tip.label==desc_tip_names)
#         print(paste(best_nodes_table$sampleID[i],anchgs$V4[match(best_nodes_table$sampleID[i],anchgs$V1)]), sep='___')
        

#     }
 
#     # bind ancient sample to final_node_in_new_tree and update new tree
#     sample_labels[i]<- paste(best_nodes_table$sampleID[i],anchgs$V4[match(best_nodes_table$sampleID[i],anchgs$V1)], sep='___')
#     if (best_nodes_table$pos[i]==0){
# 	    try(newtree<-bind.tip(tree=newtree, position = best_nodes_table$pos[i], edge.length=0.0001, where=final_node_in_new_tree,tip.label = sample_labels[i] ))
# } else {
#     try(newtree<-bind.tip(tree=newtree, position = best_nodes_table$pos[i], edge.length=best_nodes_table$pos[i], where=final_node_in_new_tree, tip.label = sample_labels[i] ))

# }

# }
# # }where=final_node_in_new_tree,position = best_nodes_table$pos[i]
# # , edge.length=best_nodes_table$pos[i],,position = best_nodes_table$pos[i],
# #     
# pdf(file="tmp2_add_ancycap.pdf", height=30, width=10)
# plot(ladderize(newtree), cex=0.3, tip.color = ifelse(newtree$tip.label %in% sample_labels, yes=2, no=1))
# # nodelabels(node=unlist(new_nodes), pch=20, cex=1, col="green")
# dev.off()


# write.tree(newtree, "newtree_added_anc.nwk")
# write.tree(tmptree, 'karmin_tree_no_ancients.nwk')
