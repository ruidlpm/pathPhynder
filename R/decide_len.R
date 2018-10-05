require(phytools)
require(scales)

# usage:
# Rscript decide.R <input_phylogeny.nwk> <input_prefix> <results_folder> 
# Description: reads tree, aDNA intree files in intree_folder, outputs result to results_folder.

args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=3) {
    stop("  Arguments needed.\n
        \tusage
        \tRscript decide.R <input_phylogeny.nwk> <results_folder> <Maximum Tolarance (int)>
        ", call.=FALSE)
} else {
    cat("   Command used:",'\n\n')
    cat(paste("decide.R", args[1], args[2],args[3], '\n\n'))
}


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


max_tolerance<-args[3]


samples<-names(br_sample_tables)

for (samp in 1:length(br_sample_tables)){
sample_name<-samples[samp]
print(sample_name)

countdata<-br_sample_tables[[samp]]

pos_score<-list()
neg_score<-list()
stopped_edges<-NULL

for (path_index in 1:length(paths)){
    path<-paths[[path_index]]
    path_pos_score<-NULL
    path_neg_score<-NULL
    stopped_path<-NULL
    for (pos1_in_path in 1:length(path)){
        entry<-countdata[which(countdata$Node1==path[pos1_in_path] & countdata$Node2==path[pos1_in_path+1]),]
        if(dim(entry)[1]!=0){
        supporting = entry$support
        notsupporting = entry$notsupport
        if (notsupporting<max_tolerance){
            path_pos_score<-c(path_pos_score,supporting)
            path_neg_score<-c(path_neg_score,notsupporting)
        } else if (notsupporting>=max_tolerance){
            path_pos_score<-c(path_pos_score,supporting)
            path_neg_score<-c(path_neg_score,notsupporting)
            stopped_edge<-entry$Edge[which(entry$Node1==path[pos1_in_path] & entry$Node2==path[pos1_in_path+1])]
            if (length(stopped_edge)>0){
                # print(paste("stop at edge",stopped_edge))
                stopped_edges<-c(stopped_edges,stopped_edge )
                break
            }
        }
    }
}


if(length(path_pos_score)>0){
    pos_score[path_index]<-data.frame(path_pos_score)
    neg_score[path_index]<-data.frame(path_neg_score)
    }
}



sums<-sapply(pos_score, sum)

print(length(which(sums==max(sums))))

if (length(which(sums==max(sums)))>1){
    best_paths<-which(sums==max(sums))
    print("TWO best paths")
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
# pos_score[[which(sapply(pos_score, sum)==max(sapply(pos_score, sum)))]]
# neg_score[[which(sapply(pos_score, sum)==max(sapply(pos_score, sum)))]]
counts_for_best_path<-countdata[match(best_path, countdata$Node2),]
counts_for_best_path<-counts_for_best_path[!is.na(counts_for_best_path$Edge),]
counts_for_best_path<-counts_for_best_path[!(counts_for_best_path$support==0 &  counts_for_best_path$notsupport==0),]
counts_for_best_path_toplot<-counts_for_best_path[!(counts_for_best_path$support==0 &  counts_for_best_path$notsupport>0),]



print(counts_for_best_path)


write.table(counts_for_best_path, file=paste0(args[2],'/',sample_name,'.txt'),quote=F, row.names=F, sep="\t")

print(best_path)

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

print(best_edge)
print(best_node1)
print(best_node2)
der_at_best<-counts_for_best_path$support[counts_for_best_path$Edge==best_edge]
anc_at_best<-counts_for_best_path$notsupport[counts_for_best_path$Edge==best_edge]

pos_at_best_edge<-NULL
if (anc_at_best==0){
    pos_at_best_edge<-0
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

# try(x <- getphylo_x(tmptree, best_node1))
# try(y <- getphylo_y(tmptree, best_node1))
# best_node1
# x
# y

# stopped_path<-br$branches[match(stopped_path, br$pos1)]


pdf(file=paste0(args[2],'/',sample_name,'.decision.pdf'), height=10, width=7)
plot((tmptree), col='darkgrey',cex=0.2, edge.col=ifelse(countdata$Edge %in% unique(stopped_edges), yes=2, no='lightgrey'), show.tip.label = T, edge.width = 1, tip.color = 1)
par(new=T)
plot((tmptree), cex=0.2, edge.col=ifelse(countdata$Edge %in% unlist(counts_for_best_path_toplot$Edge), yes=3, no=0), show.tip.label = T, edge.width = 1, tip.color = 1)
edgelabels(edge=countdata$Edge[countdata$notsupport>0],pch=20, col=alpha("red", 0.5), cex=log((countdata$notsupport[countdata$notsupport>0])+1)/2)
edgelabels(edge=countdata$Edge[countdata$support>0],pch=20, col=alpha("darkgreen", 0.7),cex=log((countdata$support[countdata$support>0])+1)/2)
edgelabels(frame="none",edge=countdata$Edge[countdata$support>0], text=countdata$support[countdata$support>0], cex=0.3)
edgelabels(frame="none",edge=countdata$Edge[countdata$notsupport>0], text=countdata$notsupport[countdata$notsupport>0], cex=0.3)

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
dev.off()
closeAllConnections()

}


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
