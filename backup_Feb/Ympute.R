require(phytools, quietly = TRUE)


cat('\n\n',"Ympute.R", '\n\n\n')




args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=3) {
    stop("\tArguments needed.\n
\tusage:
\tRscript Ympute.R <input_phylogeny.nwk> <input.vcf> <output.vcf>", call.=FALSE)
} else {
    cat("   Command used:",'\n\n')
    cat(paste("Ympute.R", args[1], args[2],args[3]), '\n\n')
}

# vcf<-read.table(args[2], comment.char="", h=T, stringsAsFactors=F)
# vcf<-read.table("to_imp.vcf", comment.char="", h=T)
# vcf<-read.table("torm.vcf", comment.char="", h=T, stringsAsFactors=F)

# vcf[(is.na(vcf))]<-'.'

# write.table(vcf, file=paste0("in_",args[2]), quote=F, sep='\t', row.names=F)

##fileformat=VCFv4.1
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">

a<-read.tree(args[1])
a$tip.label<-make.names(a$tip.label)
root_node<-a$edge[,1][1]

tree<-a

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

a<-tree
root_node<-a$edge[,1][1]

# vcf<-vcf[1,]

# vcf[vcf=="."]<-NA

# a<-read.tree("thirteencladetree.nwk")


getAncestors<-phytools:::getAncestors



report<-data.frame(position=NULL, sample_target=NULL, allele=NULL, imputed=NULL, reason=NULL, stringsAsFactors=F)

impute<-function(genotable, position_in_tree){
	# print(genotable)
	if(position_in_tree=="below"){
		print("imputing below")
		if(ancestral_allele %in% genotable[desc] | derived_allele %in% genotable[nondesc]){
			# report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=derived_allele, imputed=F, reason="snp is inconsistent with tree"))
		} else if(sum(genotable[desc] %in% derived_allele)>=2){
			# print(sum(genotable[nondesc] %in% ancestral_allele))
			genotable[[sample_target]]<-derived_allele
			# report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=derived_allele, imputed=T, reason="consistent"))
		} else if(sum(genotable[desc] %in% derived_allele)<2){
			print("cannot impute singletons/doubleton derived alleles")
			# report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=derived_allele, imputed=F, reason="not imputing doubleton derived allele"))
		} else {
			# report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=derived_allele, imputed=F, reason="not possible"))
		}
	} else if(position_in_tree=="above"){
		# print("imputing above")

		#test if tree is consistent = no ancestral alleles in desc and no derived alleles in nondesc
		if(sum(genotable[desc]==ancestral_allele, na.rm=T)>0| sum(genotable[nondesc]==derived_allele, na.rm=T)>0){
			print("inconsistent w tree")
			report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=derived_allele, imputed=F, reason="snp is inconsistent with tree"))
		} else {
			sis<-getSisters(a, origin_of_derived_allele)
			sis_desc<-getDescendants(a,sis)
			# sis_desc_tips<-NULL
			# for (nodes_and_tips in sis_desc){
			# 	sis_desc_tips<-c(sis_desc_tips,getDescendants(a,nodes_and_tips))
			# }
			# findMRCA(a,unique(sis_desc_tips))


			sis_desc<-a$tip.label[sis_desc]
			sis_desc<-sis_desc[!is.na(sis_desc)]
			ancestor<-getAncestors(a,origin_of_derived_allele)[1]
			subtree<-extract.clade(a,ancestor)
			subtree_tips<-subtree$tip.label
			# print("subtree_tips")
			# print(subtree_tips)

			if (sample_target %in% subtree_tips){
				subtree_tips<-subtree_tips[subtree_tips!=sample_target]
				#if all subtree genos (except the derived allele sample_target) is NA
				if (sum(is.na(genotable[subtree_tips]))==length(subtree_tips)){
					report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=ancestral_allele, imputed=F, reason="all_sibs_NA"))
				} else if (!ancestral_allele %in% genotable[subtree_tips]){
					report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=ancestral_allele, imputed=F, reason="not_imputing_because_there_is_der_in sibs"))
				} else if (ancestral_allele %in% genotable[subtree_tips]){
					genotable[[sample_target]]<-ancestral_allele
					# report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=ancestral_allele, imputed=T, reason="imputing_because_there_is_anc_in sibs"))
				}
			} else if(ancestor==root_node){
				print("is root")
				genotable[[sample_target]]<-ancestral_allele
				# report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=ancestral_allele, imputed=T, reason="ancestral is at root"))
			
			#if sample_target not in subtree and ancestor not in root, need to look to neighbour
			} else if (!sample_target %in% subtree_tips){
				sisanc<-getSisters(a, ancestor)
				sisanc_tips<-getDescendants(a, sisanc)
				sisanc_tip_names<-a$tip.label[sisanc_tips]
				sisanc_tip_names<-sisanc_tip_names[!is.na(sisanc_tip_names)]
				#if the sample_target is in the neighbouring branches to the derived allele
				if (sample_target %in% sisanc_tip_names){
					#are all sibs NA?
					if(sum(is.na(genotable[sisanc_tip_names]))==length(sisanc_tip_names)){
						report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=ancestral_allele, imputed=F, reason="all_sisanc_NA"))
					} else if(ancestral_allele %in% genotable[sisanc_tip_names]){
						genotable[[sample_target]]<-ancestral_allele
						# report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=ancestral_allele, imputed=T, reason="imputing_because_there_is_anc_in_sisanc"))
					}
				} else if (!sample_target %in% sisanc_tip_names){
					genotable[[sample_target]]<-ancestral_allele
					# report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=ancestral_allele, imputed=T, reason="imputing_because_sample_well_above_derived_node"))
				} else {
					# report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=derived_allele, imputed=F, reason="does not match requirements"))

				}
			}
		}
	} else {
		# report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=derived_allele, imputed=F, reason="does not match requirements"))
	}
	return(list(genotable, report))
}









# 			} else {
# 				sisanc<-getSisters(a,ancestor)
# 				sisanc<-a$tip.label[sisanc]
# 				sisanc<-sisanc[!is.na(sisanc)]
# 				print(c("sisanc", sisanc))
# 				if(length(sisanc)==0){
# 					#this means that the sister ancestral is a node, need to check desc
# 					ancestor<-getAncestors(a,origin_of_derived_allele)[2]
# 					sisanc<-getSisters(a,ancestor)
# 					side_desc<-a$tip.label[getDescendants(a,sisanc)]
# 					side_desc<-side_desc[!is.na(side_desc)]
# 					print(side_desc)
# 					if (sample_target %in% side_desc){
# 						if(sum(is.na(genotable[side_desc]))==length(side_desc)){
# 							print("do not impute!")
# 						} else {
# 							genotable[sample_target]<-ancestral_allele
# 							report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=ancestral_allele, imputed=T, reason="consistent1"))
# 						}
# 					} else {
# 						# genotable[sample_target]<-ancestral_allele
# 						report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=ancestral_allele, imputed=T, reason="consistent1"))
# 					}
# 					sisanc<-a$tip.label[sisanc]
# 					sisanc<-sisanc[!is.na(sisanc)]
# 					print(c("sisanc", sisanc))

# 					# genotable[sample_target]<-ancestral_allele
# 					# report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=ancestral_allele, imputed=T, reason="consistent2"))
# 					# print("here")
# 				}
# 				# if (status=="notImputeDerivedAllele - singleton"){

# 				# }

# 				if(length(sisanc)!=0){
# 					if (sum(sisanc %in% a$tip.label)==length(sisanc) & (sample_target %in% sisanc)){
# 						if (sum(is.na(genotable[sisanc]))==length(sisanc)){
# 						print("not imputing this one")
# 						report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=ancestral_allele, imputed=F, reason="impossible1"))
# 						}
# 					} else if (sample_target %in% sis_desc){
# 						if (sum(is.na(genotable[sis_desc]))==length(sis_desc)){
# 							print("all sis near derived_allele are NA, so not imputing those")
# 							report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=ancestral_allele, imputed=F, reason="impossible2"))
# 						}
# 					} else {
# 						# genotable[sample_target]<-ancestral_allele
# 						report<-rbind(report,data.frame(position=snp, sample_target=sample_target, allele=ancestral_allele, imputed=T, reason="consistent2"))
# 					}
# 				}
# 			}
# 		}
# 	}

# 	return(list(genotable, report))
# }




for (snp in vcf$POS){
	cat('\n\n')
	# print(snp)
	line<-vcf[vcf$POS==snp,]
	genos<-line[10:dim(line)[2]]
	try(genos[which(genos=='.')]<-NA)
	try(genos[which(genos=='N')]<-NA)


	# step1 - is the snp polymorphic?
	if(!(1 %in% genos & 0 %in% genos)){
		print("snp is not polymorphic, cannot impute singletons, and it is not informative")
		next
	} else if(sum(is.na(genos))==0){
		print("no missing")
		next
	} else {
		print(snp)

	#step2 - get origin of derived and ancestral allele
	# find node of origin of derived allele
	# find node of origin of ancestral allele (root)
	origin_of_1<-findMRCA(a,which(a$tip.label %in% colnames(genos)[which(genos==1)]))
	origin_of_0<-findMRCA(a,which(a$tip.label %in% colnames(genos)[which(genos==0)]))
	
	print(c("origin_of_1", origin_of_1))
	print(c("origin_of_0", origin_of_0))

		# define ancestral and derived alleles
		if (length(origin_of_1)==0){
			ancestral_allele<-0
			derived_allele<-1
			#this means it is a singleton, der_allele arose in a tip
			origin_of_derived_allele<-which(a$tip.label==colnames(genos)[which(genos==1)])
			status<-"notImputeDerivedAllele - singleton"

		} else if (length(origin_of_0)==0){
			ancestral_allele<-1
			derived_allele<-0
			#this means it is a singleton, der_allele arose in a tip
			origin_of_derived_allele<-which(a$tip.label==colnames(genos)[which(genos==0)])
			status<-"notImputeDerivedAllele"
		} else if (origin_of_1==origin_of_0){
			print("snp is inconsistent with the tree")
			next
			status<-"inconsistent"
		} else if (origin_of_1==root_node){
			ancestral_allele<-1
			derived_allele<-0
			origin_of_derived_allele<-origin_of_0
			status<-"impute1"
		} else if (origin_of_0==root_node){
			ancestral_allele<-0
			derived_allele<-1
			origin_of_derived_allele<-origin_of_1
			status<-"impute2"
		} else if (origin_of_0<origin_of_1){
			origin_of_0<-root_node
			ancestral_allele<-0
			derived_allele<-1
			origin_of_derived_allele<-origin_of_1
			status<-"impute3"
		} else if (origin_of_1<origin_of_0){
			origin_of_1<-root_node
			ancestral_allele<-1
			derived_allele<-0
			origin_of_derived_allele<-origin_of_0
			status<-"impute4"
		}
			print(c("status", status))

		desc<-NULL
		nondesc<-NULL

		# print(c("origin_of_derived_allele", origin_of_derived_allele))

		desc<-a$tip.label[getDescendants(a, origin_of_derived_allele)]
		desc<-desc[!is.na(desc)]
		nondesc<-a$tip.label[!a$tip.label %in% desc]


		# print(c('desc',desc))
		# print(c("origin_of_derived_allele", origin_of_derived_allele))
		# print(c("origin_of_0", origin_of_0))

		miss_samples<-colnames(genos)[which(is.na(genos))]
		for (sample_target in miss_samples){
			# print(c('sample_target', sample_target))
			if(sample_target %in% desc){
				genos<-impute(genos,"below")[[1]]
				report<-impute(genos,"below")[[2]]
			} else if (sample_target %in% nondesc){
				genos<-impute(genos,"above")[[1]]
				report<-impute(genos,"above")[[2]]


			}
		}
	}
	vcf[vcf$POS==snp,][10:dim(line)[2]]<-genos
}

# 	print("here2")

vcf[(is.na(vcf))]<-'.'
write.table(vcf, file=args[3], quote=F, sep='\t', row.names=F)
write.table(report, file="report2.txt", quote=F, sep='\t', row.names=F)
