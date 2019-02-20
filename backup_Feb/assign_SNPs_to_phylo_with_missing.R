# usage:
# Rscript assign_SNPs_to_phylo.R <input_phylogeny.nwk> <input.vcf> <out prefix>

require(phytools, quietly = TRUE)
#require(phangorn, quietly = TRUE)
suppressPackageStartupMessages(require(svMisc))


cat('\n\n',"Assigning SNPs to branches", '\n\n\n')



args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=3) {
    stop("  Arguments needed.\n", call.=FALSE)
} else {
    cat("   Command used:",'\n\n')
    cat(paste("assign_SNPs_to_phylo.R", args[1], args[2],args[3]), '\n\n')
}

#read tree
tree<-read.tree(file=args[1])
tree<-ladderize(tree)

#fix names prevent downstream problems withs sample ID
tree$tip.label<-make.names(tree$tip.label)


#read vcf
vcf<-read.table(args[2], h=T, stringsAsFactors=F)
colnames(vcf)<-make.names(colnames(vcf))


#get genotyoes
genos<- vcf[10:length(vcf)]


#check if any samples are missing and if so, exclude from vcf
miss<- colnames(vcf)[10:length(vcf)][!(colnames(vcf)[10:length(vcf)] %in% tree$tip.label)]

if(length(miss)>0){
    cat(paste0("    Number of individuals: ", dim(genos)[2]),'\n')
    cat(paste0("    Excluded samples in VCF not present in the tree: ", miss),'\n')
    vcf<-vcf[!(colnames(vcf) %in% miss)]
    cat(paste0("    Number of individuals: ", dim(genos)[2]),'\n')
}

#remove if unnecessary
#check if any samples are missing and if so, exclude from vcf
#genos<-genos[colnames(genos) %in% tree$tip.label]


    # if (length(tree$tip.label[!(tree$tip.label %in% colnames(genos))]) >0){
    #     cat(paste0("    Number of individuals: ", dim(genos)[2]),'\n')
    #     cat(paste0("    Excluded ", length(tree$tip.label[!(tree$tip.label %in% colnames(genos))]), " individuals from the tree not in VCF. (See tree_data/", args[3],".indsremoved.txt)"), '\n')
    #     write.table(file=paste0("tree_data/", args[3],".indsremoved.txt"),data.frame(tree$tip.label[!(tree$tip.label %in% colnames(genos))]), quote = F, row.names = F, col.names = F, sep='\t')
    #     cat(paste0("    Wrote new tree to tree_data/tree.", args[3],".indsremoved.txt"), '\n')
    #     cat('\n')
    #     tree<-drop.tip(tree, tree$tip.label[!(tree$tip.label %in% colnames(genos))])
    #     write.tree(file=paste0("tree_data/tree.", args[3],".indsremoved.txt"),tree)

    # }



#Remove missing data
#it would be good to fix this to impute missingness in a phylogenetically aware way
cat(paste0("    Number of SNPs: ", dim(genos)[1]),'\n')


genos[genos=="."]<-NA
genos_cc<-genos[complete.cases(genos), ]
vcf_cc<-vcf[which(complete.cases(genos)),]



#new
genos_miss<-genos[!complete.cases(genos), ]
vcf_miss<-vcf[which(!complete.cases(genos)),]
print(genos_miss)


genos_miss[] <- lapply(genos_miss, function(x) {
    as.numeric(as.character(x))
})
print(dim(genos_miss))





# if (dim(genos)[1]!=dim(genos_cc)[1]){
#     cat(paste0("    Number of SNPs after filtering: ", dim(genos_cc)[1]),'\n')
#     cat(paste0("    Number of SNPs excluded due to missingess: ",dim(genos)[1]-dim(genos_cc)[1], " (written to tree_data/",  args[3],".SNPsExcluded.txt)"),'\n\n\n')
#     write.table(file=paste0("tree_data/", args[3],".SNPsExcluded.txt"),data.frame(vcf$POS[!(vcf$POS %in% vcf_cc$POS)]), quote = F, row.names = F, col.names = F, sep='\t')

#     vcf<-vcf_cc
#     genos<-genos_cc
#     genos<-as.data.frame(genos)
#     genos<-sapply(genos, as.numeric)
#     genos<-as.data.frame(genos)

# }





# if (dim(genos)[1]!=dim(genos_cc)[1]){
#     cat(paste0("    Number of SNPs after filtering: ", dim(genos_cc)[1]),'\n')
#     cat(paste0("    Number of SNPs excluded due to missingess: ",dim(genos)[1]-dim(genos_cc)[1], " (written to tree_data/",  args[3],".SNPsExcluded.txt)"),'\n\n\n')
#     write.table(file=paste0("tree_data/", args[3],".SNPsExcluded.txt"),data.frame(vcf$POS[!(vcf$POS %in% vcf_cc$POS)]), quote = F, row.names = F, col.names = F, sep='\t')
   
#     vcf<-vcf_cc
#     genos<-genos_cc
#     genos<-as.data.frame(genos)
#     genos<-sapply(genos, as.numeric)
#     genos<-as.data.frame(genos)

# }



#get node labels
tree$node.label<-unique(tree$edge[,1])


#make data.frame with branches, and the nodes they link (pos1 and pos2)
br<-cbind(data.frame(tree$edge), c(1:length(tree$edge[,1])))
colnames(br)<-c("pos1", "pos2", "branches")


#iterates through tree from the root to tips
#gets descendants of a given node
#tests this set of descendants all have 1s and all remaining have 0s (der)
#tests this set of descendants all have 0s and all remaining have 1s (der)
#records the positions at each branch which are composed by the ALT (derpos) and REF (ancpos) allele



derpos<-list()
ancpos<-list()
total <- length(br$branches)
# pb <- txtProgressBar(min = 0, max = total, style = 3)

for (branch in br$branches){


    #print(branch)
    subset1<-br[br$branches==branch,]
    nodes<-subset1$pos1[subset1$branches==branch]
    for (node in nodes){
        subset2<-subset1[subset1$pos1==node,]
        ends<-subset2$pos2[subset2$pos1==node]
        # print(c(node,ends, branch))
        desc<-c(tree$tip.label,tree$node.label)[getDescendants(tree,ends)]
        desc<-desc[desc %in% tree$tip.label]
        der<-NULL
        anc<-NULL
        not_assigned<-NULL
        
print(dim(genos_miss))
if (dim(genos_miss)[1]>1){
        for (var in 1:dim(genos_miss)[1]){
            print(var)
            num_miss<-NULL
            num_miss<-sum(is.na(genos_miss[desc][var,]))
            num_miss_whole<-sum(is.na(genos_miss[var,]))

            if (sum(genos_miss[desc][var,], na.rm=T)==(length(desc)-num_miss) & sum(genos_miss[var,], na.rm=T)==(length(desc)-num_miss)){
               der<-c(der,vcf_miss$POS[var])
           }
                if (sum(genos_miss[desc][var,], na.rm=T)==0 & sum(genos_miss[var,], na.rm=T)==(length(colnames(genos_miss))-(num_miss_whole+length(desc)-num_miss))){
       anc<-c(anc,vcf_miss$POS[var])        
   }
        # der<- vcf$POS[which(rowSums(genos[desc], na.rm=T)==length(desc))] & rowSums(genos, na.rm=T)==length(!is.na(genos[desc])))]
        # anc<- vcf$POS[which(rowSums(genos[desc], na.rm=T)==0 & rowSums(genos, na.rm=T)==length(colnames(genos))-length(!is.na(genos[desc])))]
        #print(paste("der",length(der), "|", "anc", length(anc)))
    }


}

    derpos[[branch]]<-unique(der)
    ancpos[[branch]]<-unique(anc)


       genos_miss<- genos_miss[!(vcf_miss$POS %in% c(unlist(derpos), unlist(ancpos))),]

       vcf_miss<- vcf_miss[!(vcf_miss$POS %in% c(unlist(derpos), unlist(ancpos))),]
}

    # cat('\t\t'); setTxtProgressBar(pb, branch)
    # cat('\r',length(unique(sort(c(unlist(derpos), unlist(ancpos))))))
}


    print(derpos)



for (i in 1:length(derpos)){
    if (i>1){
        if (sum(derpos[[i]] %in% derpos[[i-1]])>0){
        derpos[[i]]<-derpos[[i]][!(derpos[[i]] %in% derpos[[i-1]])]
        }
    }
}



for (i in 1:length(ancpos)){
    if (i>1){
        if (sum(ancpos[[i]] %in% ancpos[[i-1]])>0){
        ancpos[[i]]<-ancpos[[i]][!(ancpos[[i]] %in% ancpos[[i-1]])]
        }
    }
}




# close(pb)
# cat('\n')

#save data
dir.create('tree_data', showWarnings = FALSE)

saveRDS(derpos, file=paste0('tree_data/',args[3],".derpos.RData"))
saveRDS(ancpos, file=paste0('tree_data/',args[3],".ancpos.RData"))

write.table(unlist(derpos), file=paste0('tree_data/',args[3],".derpos.txt"), quote = F, row.names = F, col.names = F, sep='\t')
write.table(unlist(ancpos), file=paste0('tree_data/',args[3],".ancpos.txt"), quote = F, row.names = F, col.names = F, sep='\t')


#make positions
pos_to_call<-as.data.frame(unique(sort(c(unlist(derpos), unlist(ancpos)))))
colnames(pos_to_call)<-"pos_to_call"


pos_to_call$chr <- 'Y'
pos_to_call$pos0<-pos_to_call$pos_to_call-1
pos_to_call$REF<-vcf$REF[match(pos_to_call$pos_to_call, vcf$POS)]
pos_to_call$ALT<-vcf$ALT[match(pos_to_call$pos_to_call, vcf$POS)]
pos_to_call$mut<-paste0(pos_to_call$REF,'->', pos_to_call$ALT)
pos_to_call$marker<-NA
pos_to_call$hg<-NA

pos_to_call<-pos_to_call[c('chr', "marker","hg", 'pos_to_call', 'mut', 'REF', 'ALT')]
# chrY    A2607   A00     7029908 C->T    C       T


bedfile_data<-data.frame(pos_to_call$chr, pos_to_call$pos_to_call-1,pos_to_call$pos_to_call)

#write these in pos to call format and in bed format
write.table(file=paste0('tree_data/',args[3],'.sites.txt'),pos_to_call, quote = F, row.names = F, col.names = F, sep='\t')
write.table(file=paste0('tree_data/',args[3],'.sites.bed'),bedfile_data, quote = F, row.names = F, col.names = F, sep='\t')

cat(paste0("\t",dim(unique(pos_to_call))[1]," informative positions for variant calling (written to tree_data/", args[3],".bed)"),'\n\n')

not_added<-vcf$POS[!vcf$POS %in% pos_to_call$pos_to_call]
write.table(file=paste0('tree_data/',args[3],'.not_added.txt'),not_added, quote = F, row.names = F, col.names = F, sep='\t')


#Rscript assign_SNPs_to_phylo.R /Users/rm890/Integrating_Y/karmin_data/RAxML_bestTree.karmin.aln.mGTR.50boot_try2 /Users/rm890/Botai_project/Karmin/try_karmin/karmin_w_changes.vcf 50boot_try2



# chrY    M2684   IJK 7792789 G->A    G   A

# 1   7792789 Y:7792789   A   G   .   .   .   GT  0   0   1   1   1   0   0   1   .   1  1




# 1   7629583 Y:7629583   G   A   .   .   .   GT  0   0   1   1   1   0   0   1   1   1  1
# 1   6753519 Y:6753519   G   A   .   .   .   GT  0   0   1   1   1   0   0   1   1   1  1
# 1   7173143 Y:7173143   A   G   .   .   .   GT  0   0   1   1   1   0   0   1   1   1  1
# 1   7702973 Y:7702973   A   T   .   .   .   GT  0   0   1   1   1   0   0   1   1   1  1
# 1   21571895 Y:21571895  A   G   .   .   .   GT  0   0   1   1   1   0   0   1   1  1   1



# removed these inconsistent positions
# cat ~/software/y_leaf/data/positions.txt| cut -f4 | awk '{print "Y\t"$1-1"\t"$1}' > tmp.bed
# bedtools getfasta -fi /Users/rm890/hs37d5.fa -bed tmp.bed | grep -v \> | paste ~/software/y_leaf/data/positions.txt - |tr "\015" "\n" > test
# # 18111802 19041163 21132987  7763738
# a<-read.table("test")
# tmp<-a[as.character(a$V6)!=as.character(a$V8),]
# a<-read.table("test")
# a$strand<-ifelse(as.character(a$V6)==as.character(a$V8), yes="For", no="Rev")
# a<-a[!(a$V4 %in% tmp$V4[as.character(tmp$V7)!=as.character(tmp$V8)]),]
# write.table(a,file="positions_and_strand.txt")



