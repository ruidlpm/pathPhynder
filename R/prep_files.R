# pathPhynder
# Author: Rui Martiniano
# Contact: rm890 [at] cam.ac.uk
# usage: Rscript prep_files.R <input_phylogeny.nwk> <branches.snp> <out prefix> <haplogroups(optional)>

require(phytools, quietly = TRUE)
options(scipen=999)



cat('\n\n',"Preparing files for sample placement", '\n\n\n')


args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=4) {
    stop("\tArguments needed.\n
    \tusage:
    \tRscript prep_files.R <input_phylogeny.nwk> <branches.snp> <out prefix> <haplogroups(optional)>", call.=FALSE)
} else {
    cat("   Command used:",'\n\n')
    cat(paste("prep_files.R", args[1], args[2],args[3]), '\n\n')
}



tree<-read.tree(args[1])

branches_file<-read.table(args[2], fill=T, stringsAsFactors=F)
branches_file<-branches_file[which(branches_file$V1=='c'):length(branches_file$V1),]
branches_file<-branches_file[1:7]

call_hgs_table<-args[4]

tree<-read.tree(args[1])
results_folder<-args[3]


make_files<-function(branches_file, call_hgs_table){

    chr<-gsub('chr','',as.character(branches_file[branches_file$V1=='c',]$V4))
    variants<-branches_file[branches_file$V1=='V',]

    branches<-branches_file[branches_file$V1=='B',]

    ancestral<-branches_file[branches_file$V1=='A',]

    variants$V2<-as.numeric(as.character(variants$V2))
    branches$V2<-as.numeric(as.character(branches$V2))


    strand<- ancestral$V2
    strand[strand==0]<-'+'
    strand[strand==1]<-'-'


    df<-data.frame(chr=chr,pos=variants$V2, REF=toupper(variants$V5), ALT=toupper(variants$V7), status=strand, branch=branches$V2)

    df$mutation<-paste0(df$REF,'->',df$ALT)

    if (call_hgs_table!='none'){

        #if hg table is passed, then assign haplogroups to tree branches
        
        snps<-read.table(call_hgs_table, h=F, stringsAsFactors=F)
        
        df$marker<-snps$V1[match(df$pos, snps$V3)]
        df$hg<-snps$V2[match(df$pos, snps$V3)]
        df$hg_anc_allele<-snps$V4[match(df$pos, snps$V3)]
        df$hg_der_allele<-snps$V5[match(df$pos, snps$V3)]
        
        
        df$hg_strand<-NA
        df$hg_strand[which(df$REF==df$hg_anc_allele & df$ALT==df$hg_der_allele)]<-'+'
        df$hg_strand[which(df$ALT==df$hg_anc_allele & df$REF==df$hg_der_allele)]<-'-'
        
        
        #excluded hgs whose alleles do not match the vcf alleles
        df2<-df[!is.na(df$hg),]
        pos_to_excl_hg<-df2[is.na(df2$hg_strand),]$pos
                
        df[df$pos %in% pos_to_excl_hg,]$hg_strand<-'mismatch'

    
        df<-rbind(df[c('chr', 'marker', 'hg', 'pos','mutation','REF', 'ALT', 'branch', 'status', 'hg_anc_allele', 'hg_der_allele', 'hg_strand')])
    
    } else {

        df$marker<-NA
        df$hg<-NA
        df$hg_anc_allele<-NA
        df$hg_der_allele<-NA     
        df$hg_strand<-NA

        df<-rbind(df[c('chr', 'marker', 'hg', 'pos','mutation','REF', 'ALT', 'branch', 'status', 'hg_anc_allele', 'hg_der_allele', 'hg_strand')])

    }

    df_bed<-data.frame(df$chr, df$pos-1, df$pos)
    df_chrbed<-data.frame(paste0('chr',df$chr), df$pos-1, df$pos)
    
    dir.create('tree_data', showWarnings = FALSE)
    
    write.table(file=paste0('tree_data/',args[3],'.sites.bed'),df_bed, quote = F, row.names = F, col.names = F, sep='\t')
    write.table(file=paste0('tree_data/',args[3],'.siteschr.bed'),df_chrbed, quote = F, row.names = F, col.names = F, sep='\t')
    write.table(file=paste0('tree_data/',args[3],'.sites.txt'),df, quote = F, row.names = F, col.names = T, sep='\t')
    
    cat(paste0("\t",dim(unique(df_chrbed))[1]," informative positions for variant calling (written to tree_data/", args[3],".sites.bed)"),'\n\n')

return(df)
}




make_edge_df<-function(tree, call_hgs_table, sites_df){

    edge_df <- data.frame(tree$edge)
    colnames(edge_df) <- c("Node1","Node2")

    edge_df$hg <- NA
    edge_df$Edge <- rownames(edge_df)

    varcount<-as.data.frame(table(sites_df$branch))
    colnames(varcount)<-c('Edge','count')

    edge_df$snp_count<- as.numeric(as.character(varcount$Edge[match(edge_df$Edge, varcount$Edge)]))
    edge_df$snp_count[is.na(edge_df$snp_count)]<-0

    if (call_hgs_table!='none'){
            # remove mismatiching snps and snps with no known haplogroup
            sites_df<-sites_df[!is.na(sites_df$hg),]
            sites_df<-sites_df[sites_df$hg_strand!='mismatch',]
           # sites_df<- sites_df[as.character(sites_df$hg_strand)==as.character(sites_df$status),]

        for(i in 1:length(edge_df$Edge)){
        # print(i)
            tmp<-sites_df[sites_df$branch==edge_df$Edge[i],]
            if (dim(tmp)[1]>0){
                edge_df$hg[which(edge_df$Edge==edge_df$Edge[i])]<- paste(unique(tmp$hg[!is.na(tmp$hg)]), collapse=";")
            } else {
                edge_df$hg[which(edge_df$Edge==edge_df$Edge[i])]<-''
            }
        }
    } else {
        edge_df$hg<-NA
    }
    edge_df <- edge_df[c('Edge','Node1','Node2','hg','snp_count')]
    write.table(edge_df,file=paste0("tree_data/",args[3],".edge_df.txt"), quote=F, row.names=F, sep='\t')
    
    cat(paste0("\t","table with branch information written to tree_data/", args[3],".edge_df.txt)"),'\n\n')

}


sites_df<-make_files(branches_file,call_hgs_table)
make_edge_df(tree, call_hgs_table, sites_df)






