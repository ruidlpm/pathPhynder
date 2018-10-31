# usage:
# Rscript assign_SNPs_to_phylo.R <input_phylogeny.nwk> <input.vcf> <out prefix>

require(phytools, quietly = TRUE)
require(phangorn, quietly = TRUE)


cat('\n\n',"Phylogenetically aware imputation", '\n\n\n')



args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=3) {
    stop("\tArguments needed.\n
\tusage:
\tRscript assign_SNPs_to_phylo.R <input_phylogeny.nwk> <input.vcf> <output.vcf", call.=FALSE)
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

vcf<-get_vcf(args[2])

# print(vcf)

#check if any samples are missing and if so, exclude from vcf
miss<- colnames(vcf)[10:length(vcf)][!(colnames(vcf)[10:length(vcf)] %in% tree$tip.label)]

if(length(miss)>0){
    cat(paste0("    Number of individuals: ", dim(vcf)[2]-10),'\n')
    cat(paste0("    Excluded samples in VCF not present in the tree: ", miss),'\n')
    vcf<-vcf[!(colnames(vcf) %in% miss)]
    cat(paste0("    Number of individuals: ", dim(vcf)[2]-10),'\n')
}

#remove if unnecessary
#check if any samples are missing and if so, exclude from vcf
# vcf<-vcf[colnames(genos) %in% tree$tip.label]


if (length(tree$tip.label[!(tree$tip.label %in% colnames(vcf))]) >0){
    cat(paste0("    Number of individuals: ", dim(vcf)[2]-10),'\n')
    cat(paste0("    Excluded ", length(tree$tip.label[!(tree$tip.label %in% colnames(vcf))]), " individuals from the tree not in VCF. (See tree_data/", args[3],".indsremoved.txt)"), '\n')
    write.table(file=paste0("tree_data/", args[3],".indsremoved.txt"),data.frame(tree$tip.label[!(tree$tip.label %in% colnames(vcf[10:dim(vcf)[2]]))]), quote = F, row.names = F, col.names = F, sep='\t')
    cat(paste0("    WARNING: Wrote new tree to tree_data/tree.", args[3],".indsremoved.txt"), '\n')
    cat('\n')
    tree<-drop.tip(tree, tree$tip.label[!(tree$tip.label %in% colnames(vcf))])
    write.tree(file=paste0("tree_data/tree.", args[3],".indsremoved.txt"),tree)

}

#Remove missing data
#it would be good to fix this to impute missingness in a phylogenetically aware way
#I have fixed it, missing data is now being accounted for and used
samples<-colnames(vcf[10:dim(vcf)[2]])

vcf[samples]<-sapply(vcf[samples], as.numeric)

complete_vcf<-vcf
complete_vcf[samples][complete_vcf[samples]=='N']<-NA
complete_vcf<-complete_vcf[complete.cases(complete_vcf[samples]),]

vcf_with_missing<-vcf[!vcf$POS %in% complete_vcf$POS,]

cat(paste0("    Number of SNPs with any missing data: ", dim(vcf)[1]),'\n')
cat(paste0("    Number of SNPs with no missing data: ", dim(complete_vcf)[1]),'\n')

edges<-data.frame(tree$edge)
colnames(edges)<-c('pos1','pos2')
edges$edge<-rownames(edges)

imputed_pos<-NULL
# imputed_vcf<-matrix(nrow=0, ncol=dim(vcf_with_missing)[2])
# imputed_vcf<-data.frame(imputed_vcf)
# colnames(imputed_vcf)<-colnames(vcf_with_missing)
# print(imputed_vcf)
for (edge in edges$edge){
    print(c(edge,'/', length(edges$edge)))
    relevant_node<-edges$pos2[edges$edge==edge]
    desc<-tree$tip.label[getDescendants(tree, relevant_node)][!is.na(tree$tip.label[getDescendants(tree,relevant_node)])]
    nondesc<- samples[!samples %in% desc]
    # vcf_with_missing$na_count_samples <- apply(vcf_with_missing[samples], 1, function(x) sum(is.na(x)))
    vcf_with_missing$na_count_desc <- apply(vcf_with_missing[desc], 1, function(x) sum(is.na(x)))
    # vcf_with_missing$na_count_nondesc <- apply(vcf_with_missing[nondesc], 1, function(x) sum(is.na(x)))
    

    # of the descendants, get which ones are tips
    relevant_tips<-tree$tip.label[unlist(Descendants(tree,node=relevant_node, type='tips'))]
    # get the inner and outer tips
    outer_tip<-relevant_tips[which(relevant_tips %in% tree$tip.label[unlist(Descendants(tree,node=relevant_node, type='children'))])]
    inner_tips<-relevant_tips[which(!relevant_tips %in% tree$tip.label[unlist(Descendants(tree,node=relevant_node, type='children'))])]
           

    #if the length of the outer tip is not 1, pass
    if (!(length(outer_tip)>1 | length(outer_tip)==0)){
        # vcf_with_missing<-vcf_with_missing[!vcf_with_missing$POS %in% imputed_pos,]

        outerALTvcf<-vcf_with_missing[which(vcf_with_missing[[outer_tip]]==1),]
        outerREFvcf<-vcf_with_missing[which(vcf_with_missing[[outer_tip]]==0),]

    if(sum(outerREFvcf$na_count_desc) >0 | sum(outerALTvcf$na_count_desc) >0){

        genos_inner_ALT<-outerALTvcf[which(colnames(outerALTvcf) %in% inner_tips)]
        genos_inner_REF<-outerREFvcf[which(colnames(outerREFvcf) %in% inner_tips)]
        # non_imputable_ALT<-(outerALTvcf$POS[which(genos_inner_ALT == 0)])


        is_present<-function(data, geno){
            excl<-NULL
            for (snp in 1:dim(data)[1]){
                try(tmpsum<-sum(data[snp,]==geno, na.rm=T))
                if (tmpsum>0){
                    excl<-c(excl, snp)
                }
            }
            return(excl)
        }

        # maybe_imputable_ALT<-(outerALTvcf$POS[which(genos_inner_ALT == 1)])
        # maybe_imputable_REF<-(outerREFvcf$POS[which(genos_inner_REF == 0)])

        # maybe_imputable_ALT<-(outerALTvcf$POS[which(genos_inner_ALT == 1)])
        # maybe_imputable_REF<-(outerREFvcf$POS[which(genos_inner_REF == 0)])
        # print(dim(genos_inner_ALT))
        # print(dim(genos_inner_REF))
        #         print((genos_inner_ALT))
        # print((genos_inner_REF))

        # maybe_imputable_ALT<-outerALTvcf[!outerALTvcf$POS %in% outerALTvcf$POS[is_present(genos_inner_ALT, 0)],]
        # maybe_imputable_REF<-outerALTvcf[!outerREFvcf$POS %in% outerREFvcf$POS[is_present(outerREFvcf, 1)],]

# print(genos_inner_REF)

    if (dim(genos_inner_REF)[1]>0){
            maybe_imputable_REF<-outerREFvcf$POS[!outerREFvcf$POS %in% outerREFvcf$POS[is_present(genos_inner_REF, 1)]]
    # print(maybe_imputable_REF)
    } else {
        maybe_imputable_REF<-NULL
    }
    
    if (dim(genos_inner_ALT)[1]>0){
        maybe_imputable_ALT<-outerALTvcf$POS[!outerALTvcf$POS %in% outerALTvcf$POS[is_present(genos_inner_ALT, 0)]]
    # print(maybe_imputable_ALT)
    
    } else {
        maybe_imputable_ALT<-NULL
    }




        # non_imputable_REF<-(outerREFvcf$POS[which(genos_inner_REF == 1)])

        # if(length(maybe_imputable_ALT)>0 & sum(imputed_pos %in% maybe_imputable_ALT)==0 ){
        if(length(maybe_imputable_ALT)>0){
            # print('relevant_node')
            # print(relevant_node)
            # print('relevant_tips')
            # print(relevant_tips)
            # print('outer_tip')
            # print(outer_tip)
            # print('inner_tips')
            # print(inner_tips)

            # print(c('postoimpute',paste(maybe_imputable_ALT)))
            target<-data.frame(outerALTvcf[outerALTvcf$POS %in% maybe_imputable_ALT,])

            target<-target[which(target$na_count_desc>0),]


            for (snp in target$POS ){
                                print(snp)

                tmp<-target[target$POS==snp,]
                # print(tmp)
                # if(!1 %in% tmp[nondesc]){
                    tmp[which(colnames(tmp) %in% inner_tips)][which(is.na(tmp[which(colnames(tmp) %in% inner_tips)]))]<-1
                    target[target$POS==snp,]<-tmp
                    vcf_with_missing[vcf_with_missing$POS==snp,]<-tmp
                    imputed_pos<-c(imputed_pos, snp)
                    # print(tmp[which(colnames(tmp) %in% inner_tips)])

                    # imputed_vcf<-rbind(imputed_vcf, tmp)
                    # print(tmp)
                    # if (snp %in% imputed_vcf$POS){
                    #     imputed_vcf[imputed_vcf$POS==snp,]<-tmp
                    # } else {
                    #     imputed_vcf<-rbind(imputed_vcf, tmp)
                    # }
                # }
            }
        }


        # if(length(maybe_imputable_REF)>0 & sum(imputed_pos %in% maybe_imputable_REF)==0 ){
        if(length(maybe_imputable_REF)>0){
            # print('relevant_node')
            # print(relevant_node)
            # print('relevant_tips')
            # print(relevant_tips)
            # print('outer_tip')
            # print(outer_tip)
            # print('inner_tips')
            # print(inner_tips)

            # print(c('postoimpute',paste(maybe_imputable_REF)))
            target<-data.frame(outerREFvcf[outerREFvcf$POS %in% maybe_imputable_REF,])
                        target<-target[which(target$na_count_desc>0),]

            for (snp in target$POS ){
                print(snp)
                tmp<-target[target$POS==snp,]
                # print(tmp)
                # if(!0 %in% tmp[nondesc]){
                    tmp[which(colnames(tmp) %in% inner_tips)][which(is.na(tmp[which(colnames(tmp) %in% inner_tips)]))]<-0
                    target[target$POS==snp,]<-tmp
                    vcf_with_missing[vcf_with_missing$POS==snp,]<-tmp
                    imputed_pos<-c(imputed_pos, snp)
                    # print(as.vector(rowSums(tmp[which(colnames(tmp) %in% inner_tips)])))

                    # # print(tmp)
                    # if (snp %in% imputed_vcf$POS){
                    #     imputed_vcf[imputed_vcf$POS==snp,]<-tmp
                    # } else {

                    # imputed_vcf<-rbind(imputed_vcf, tmp)
                    # }
                # }
            }
        }

}
}

}
vcf_with_missing[which(vcf_with_missing==NA)]<-'.'
# vcf_with_missing$na_count_samples<-NULL
vcf_with_missing$na_count_desc<-NULL
# vcf_with_missing$na_count_nondesc<-NULL
print(unique(sort(imputed_pos)))
print(length(unique(sort(imputed_pos))))
# print(vcf_with_missing)

all_vcf<-rbind(complete_vcf, vcf_with_missing)

all_vcf<-all_vcf[order(all_vcf$POS),]

write.table(all_vcf,file="imputed.vcf", quote=F, row.names=F, sep='\t')



