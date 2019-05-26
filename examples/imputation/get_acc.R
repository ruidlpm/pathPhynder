
get_vcf <- function(vcf_name){
    all_content = readLines(vcf_name)
    skip = all_content[-c(grep("CHROM",all_content))]
    vcf <- read.table(textConnection(skip), stringsAsFactors=F)
    header <- unlist(strsplit(all_content[grep("CHROM", all_content)[length(grep("CHROM", all_content))]], '\t'))
    colnames(vcf) <- make.names(header)
    #if T alleles read as TRUE, convert to character T.
    vcf$REF[vcf$REF==TRUE]<-"T"
    vcf$ALT[vcf$ALT==TRUE]<-"T"
    all_content<-NULL
    return(vcf)
}




vcf<-get_vcf("simulated_data.vcf")
vcf_miss<-get_vcf("simulated_data.miss0.1.vcf")
vcf_imp<-get_vcf("simulated_data.miss0.1.imputed.vcf")



get_position_miss<-function(x){
	vect<-list()
	for (i in 1:dim(x)[1]){
		if(sum(x[10:dim(x)[2]][i,]=='.') >0){
			vect[[i]]<-rep(x$POS[i],sum(x[10:dim(x)[2]][i,]=='.'))
		}
	}
	return(vect)
}



poss<-unlist(get_position_miss(vcf_miss))




res<-data.frame(pos=poss,
miss=as.matrix(vcf_miss[10:dim(vcf_miss)[2]])[which(vcf_miss[10:dim(vcf_miss)[2]]=='.')],
imp=as.matrix(vcf_imp[10:dim(vcf_imp)[2]])[which(vcf_miss[10:dim(vcf_miss)[2]]=='.')],
real=as.matrix(vcf[10:dim(vcf)[2]])[which(vcf_miss[10:dim(vcf_miss)[2]]=='.')])



#if any genotype is incorrectly imputed, print line 
for(i in 1:dim (vcf)[1]){
	if ( sum(vcf[10:dim(vcf)[2]][i,]==1 & vcf_imp[10:dim(vcf)[2]][i,]==0)){
		print(i)
	} else if (sum(vcf[10:dim(vcf)[2]][i,]==0 & vcf_imp[10:dim(vcf)[2]][i,]==1))
		print(i)
}


 

percentage_imputed<-(sum(res$imp!='.')*100)/length(res$miss)



imp_df<-res[res$imp!='.',]


sum(imp_df$imp==imp_df$real)*100/length(imp_df$real)

imputation_accuracy<-sum(imp_df$imp==imp_df$real)*100/length(imp_df$real)


imp_df_der<- imp_df[imp_df$real==1,]



sum(imp_df_der$imp==imp_df_der$real)*100/length(imp_df_der$real)

imputation_accuracy_der<-sum(imp_df_der$imp==imp_df_der$real)*100/length(imp_df_der$real)


imp_df_anc<- imp_df[imp_df$real==0,]


sum(imp_df_anc$imp==imp_df_anc$real)*100/length(imp_df_anc$real)

imputation_accuracy_anc<-sum(imp_df_anc$imp==imp_df_anc$real)*100/length(imp_df_anc$real)


print(data.frame(percentage_imputed, imputation_accuracy_der, imputation_accuracy_anc))




#imp_df[imp_df$imp!=imp_df$real,]







#miss<-vcf_miss[vcf_miss$POS==888543,]
#miss<-miss[10:dim(miss)[2]]


#imp<-vcf_imp[vcf_imp$POS==888543,]
#imp<-imp[10:dim(imp)[2]]



#real<-vcf[vcf$POS==888543,]
#real<-real[10:dim(real)[2]]


#miss[which(miss=='.')]
#imp[which(miss=='.')]
#real[which(miss=='.')]


 #a<-read.tree("simulated_data.nwk")
#plot(a)

#tiplabels(tip=match(names(real)[which(real==1)], a$tip.label), col=3, text=1)

#tiplabels(tip=match(names(real)[which(real==0)], a$tip.label), col=2, text=0)


