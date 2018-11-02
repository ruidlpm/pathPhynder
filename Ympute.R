require(phytools)
# vcf<-read.table("missing_13clade_3.vcf", comment.char="", h=T)
vcf<-read.table("tmp", comment.char="", h=T, stringsAsFactors=F)

vcf$ALT<-as.character(vcf$ALT)
vcf$REF<-as.character(vcf$REF)
vcf$ALT[which(vcf$ALT==TRUE)]<-'T'
vcf$REF[which(vcf$REF==TRUE)]<-'T'

# a<-read.tree("thirteencladetree.nwk")
a<-read.tree("/Users/rm890/Integrating_Y/paper_data/Karmin/test.nwk")
a$tip.label<-make.names(a$tip.label)
colnames(vcf)<-make.names(colnames(vcf))
root_node<-a$edge[,1][1]

getAncestors<-phytools:::getAncestors


report<-data.frame(position=NULL, sample=NULL, allele=NULL, imputed=NULL, reason=NULL)

impute<-function(genotable, position_in_tree){
	# print(genotable)
	if(position_in_tree=="below"){
		print("imputing below")
		if(ancestral_allele %in% genotable[desc] | derived_allele %in% genotable[nondesc]){
			print("inconsistent")
			report<-rbind(report,data.frame(position=snp, sample=sample, allele=derived_allele, imputed=F, reason="inconsistent"))
		} else if(sum(genotable[desc] %in% derived_allele)>=2){
			print(sum(genotable[nondesc] %in% ancestral_allele))
			genotable[sample]<-derived_allele
			report<-rbind(report,data.frame(position=snp, sample=sample, allele=derived_allele, imputed=T, reason="consistent"))

		} else if(sum(genotable[desc] %in% derived_allele)<2){
			print("cannot impute singletons/doubletons")
			report<-rbind(report,data.frame(position=snp, sample=sample, allele=derived_allele, imputed=F, reason="doubleton"))

		}
	} else if(position_in_tree=="above"){
		print("imputing above")
		if(ancestral_allele %in% genotable[desc] | derived_allele %in% genotable[nondesc]){
			print("inconsistent")
			report<-rbind(report,data.frame(position=snp, sample=sample, allele=derived_allele, imputed=F, reason="inconsistent"))
		} else {
			sis<-getSisters(a, origin_of_derived_allele)
			sis_desc<-getDescendants(a,sis)
			sis_desc<-a$tip.label[sis_desc]
			sis_desc<-sis_desc[!is.na(sis_desc)]
			if (sample %in% sis_desc){
				if (sum(is.na(genos[sis_desc]))==length(sis_desc)){
					print("all sis near derived_allele are NA, so not imputing those")
					report<-rbind(report,data.frame(position=snp, sample=sample, allele=derived_allele, imputed=F, reason="impossible"))
					# nondesc<-nondesc[!nondesc %in% sis_desc]
				}
			} else {
				genotable[sample]<-ancestral_allele
				report<-rbind(report,data.frame(position=snp, sample=sample, allele=derived_allele, imputed=T, reason="consistent"))
			}
		}
	}
	# print(genotable)
	# print(warnings())
	return(list(genotable, report))
}



for (snp in vcf$POS){
	newgenos<-NULL
	cat('\n\n\n\n')
	print(snp)
	line<-vcf[vcf$POS==snp,]
	genos<-line[10:dim(line)[2]]
	try(genos[which(genos=='.')]<-NA)
	# step1 - is the snp polymorphic?
	if(!(1 %in% genos & 0 %in% genos)){
		print("snp is not polymorphic, cannot impute singletons, and it is not informative")
	} else {
	
	#step2 - get origin of derived and ancestral allele
	# find node of origin of derived allele
	# find node of origin of ancestral allele (root)
	origin_of_1<-findMRCA(a,which(a$tip.label %in% colnames(genos)[which(genos==1)]))
	origin_of_0<-findMRCA(a,which(a$tip.label %in% colnames(genos)[which(genos==0)]))
	
		# define ancestral and derived alleles
		if (length(origin_of_1)==0){
			ancestral_allele<-0
			derived_allele<-1
			origin_of_derived_allele<-which(genos==1)
			status<-"notImputeDerivedAllele - singleton"
		} else if (length(origin_of_0)==0){
			ancestral_allele<-1
			derived_allele<-0
			origin_of_derived_allele<-which(genos==0)
			status<-"notImputeDerivedAllele"
		} else if (origin_of_1==origin_of_0){
			print("snp is inconsistent with the tree")
			next
			status<-"inconsistent"
		} else if (origin_of_1==root_node){
			ancestral_allele<-1
			derived_allele<-0
			origin_of_derived_allele<-origin_of_0
			status<-"impute"

		} else if (origin_of_0==root_node){
			ancestral_allele<-0
			derived_allele<-1
			origin_of_derived_allele<-origin_of_1
			status<-"impute"
		}
		
		desc<-a$tip.label[getDescendants(a, origin_of_derived_allele)]
		desc<-desc[!is.na(desc)]
		nondesc<-a$tip.label[!a$tip.label %in% desc]


		print(length(a$tip.label))
		print(length(nondesc))
		print(length(desc))
		print(ancestral_allele)
		print(derived_allele)


		miss_samples<-colnames(genos)[which(is.na(genos))]
		for (sample in miss_samples){
			print(c('sample', sample))
			if(sample %in% desc){
				genos<-impute(genos,"below")[[1]]
				report<-impute(genos,"below")[[2]]
			} else if (sample %in% nondesc){
				genos<-impute(genos,"above")[[1]]
				report<-impute(genos,"below")[[2]]
			}

		}
	}
	vcf[vcf$POS==snp,][10:dim(line)[2]]<-genos
}

print(report)

# vcf[which(is.na(vcf))]<-'.'
write.table(vcf, file="otherimp.vcf", quote=F, sep='\t', row.names=F)
write.table(report, file="report.txt", quote=F, sep='\t', row.names=F)




# 	# samples downstream of node of origin of derived allele
# 	if (origin2==root_node){
# 		nondesc<-a$tip.label[getDescendants(a, origin)]
# 		nondesc<-nondesc[!is.na(nondesc)]
# 		# samples outside of this tree
# 		desc<-a$tip.label[!a$tip.label %in% nondesc]
# 	} else {
# 		desc<-a$tip.label[getDescendants(a, origin)]
# 		desc<-desc[!is.na(desc)]
# 		# samples outside of this tree
# 		nondesc<-a$tip.label[!a$tip.label %in% desc]
# 	}





# 		miss_samples<-colnames(genos)[which(is.na(genos))]
# 		for (sample in miss_samples){
# 			print(c('sample', sample))







# step3 - assess if sample is a suitable candidate for imputation







# reset excluded_vector
# is sample above or below the node of origin?

# if below, we're dealing with a derived allele
# 	if the tips below do not contain an ancestral allele and if the tips above do not contain a derived allele
# 		if the derived allele not a singleton?
# 			impute derived in desc NA

# 			get sister nodes
# 			for each sister node, 
# 				if tips contain an allele (this has to be ancestral)
# 					impute all NAs above as ancestral
# 				if tips contain only NAs
# 					do not impute and exclude those samples from imputation (add these to a excluded_vector)

# if above, we're dealing with an ancestral allele
# 	if the tips below the node of origin of the derived allele do not contain an ancestral allele and if the tips above do not contain a derived allele	
# 		if the sample in question is not in the excluded vector
# 			impute all nondesc NAs as ancestral







# reasons:
# the site is not polymorphic in the dataset
# the site is not consistent with the tree topology
# the derived allele is a singleton
# it's not possible to know if this sample is derived or ancestral



# position=NA, sample=NA, allele=NA, imputed=NA, reason=NA,
