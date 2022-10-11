#!/usr/bin/env Rscript

# usage:
# ancient_SNPs_to_branches.R <input_phylogeny.nwk> <prefix>.Rdata <intree_folder> <results_folder>

# Description: reads tree, aDNA intree files in intree_folder, outputs result to results_folder.

suppressPackageStartupMessages(require(phytools))
suppressPackageStartupMessages(require(scales))

cat('\n\n',"ancient_SNPs_to_branches.R", '\n\n\n')

args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=4) {
    stop("  Arguments needed.\n
        \tusage
        \tancient_SNPs_to_branches.R <input_phylogeny.nwk> <prefix>.Rdata <intree_folder> <results_folder>
        ", call.=FALSE)
}


tree_file=args[1]
derpos_file<-paste0(args[2],".derpos.RData")
ancpos_file<-paste0(args[2],".ancpos.RData")


for (testfile in c(tree_file, derpos_file, ancpos_file)){
    if (!file.exists(testfile)) {
        stop(paste(testfile, "- does this file exist?"))
    }
}


edge_df<-read.table(paste0(args[2],".edge_df.txt"), h=T, stringsAsFactors=F, sep='\t')

lens<-NULL
for (i in edge_df$Edge){
    tmp<-edge_df[edge_df$Edge==i,]
    lens[i]<-(length(unique(unlist(strsplit(tmp$positions[1], '\\;')))))
}

nums<-cbind(edge_df, number_of_positions=lens)

tree<-read.tree(file=tree_file)
tree<-ladderize(tree)

#make R readables names to prevent downstream problems withs sample ID
tree$tip.label<-make.names(tree$tip.label)

derpos<-readRDS(file=derpos_file)
ancpos<-readRDS(file=ancpos_file)

dir.create(args[4], showWarnings = FALSE)

br<-read.table(paste0(args[2],".edge_df.txt"), h=T, sep='\t')
br<-br[c('Edge','Node1','Node2','hg')]

file_list<-list.files(args[3], pattern ='intree.txt')
samps<-file_list

cat("\tProcessing",length(samps), "samples", '\n\n\n')


read_in<-function(sample_name){
    info <-file.info(paste0(args[3],"/",sample_name))
    if (info$size==0){
        cat(paste0(args[3],"/",sample_name), "is empty\n\n")
    } else {
        ancvcf<-try(read.table(paste0(args[3],"/",sample_name), stringsAsFactors = F))
        if (inherits(ancvcf, 'try-error')){
            stop("Something wrong with the intree.txt filenames, not being able to read them.")
        } else {
            if(dim(ancvcf)[2]!=6){
                stop("Expected 6 columns in intree.txt, instead found", dim(ancvcf)[2])
            } else if(dim(ancvcf)[2]==6){
                cat(paste0("\tRead ", dim(ancvcf)[1], " positions from file ", sample_name,"\n\n"))
                colnames(ancvcf)[1]<-"POS"
                colnames(ancvcf)[6] <-make.names(sample_name)
                return(ancvcf)
            }
        }
    }
}


#annotating tree branches with derived and ancestral alleles observed at each branch

excluded<-NULL
br_sample_tables<-list()


# for each sample
for(samp in samps){
    ancvcf<-(read_in(samp))
    samp<-make.names(samp)
        if (is.data.frame(ancvcf)==F){
            print("ancient sample intree.txt file is empty")
            excluded<-c(excluded, samp)
            samps<-samps[samps!=samp]
        }  else {
            ancient_sample<-make.names(samp)
            subset_vcf<-data.frame(ancvcf$POS, ancvcf[ancient_sample])
            colnames(subset_vcf)[1]<-'POS'
            #get informative positions overlapping with ancient sample
            all_pos<-unique(c(unlist(derpos), unlist(ancpos)))
            ancient<-subset_vcf[subset_vcf$POS %in% all_pos,]
            colnames(ancient)<-c("POS",ancient_sample)
            ##################################################################
            # get supporting alleles for inclusion in branches
            #  if branch defined by 0 at a marker, and if aDNA sample has 0
            #      or
            #  if branch defined by 1 at a marker, and if aDNA sample has 1
            ##################################################################
            #count supporting alleles for inclusion at each branch
            support_val<-list()
            #for each branch
            for (branch in br$Edge){
                #count number of ancestral alleles '0' at branch defining positions also ancestral '0'
                #count number of derived alleles '1' at branch defining positions also derived '1'
                support_count<-sum(ancient$POS[ancient[ancient_sample]==0] %in% ancpos[[branch]] ) + sum(ancient$POS[ancient[ancient_sample]==1] %in% derpos[[branch]] )
                       
                #add count of suportting positons at each branch for the current sample
                if(support_count>0){
                    support_val[[branch]]<-support_count
                } else if (support_count==0) {
                    support_val[[branch]]<-0
                }
            }
            #count alleles not supporting inclusion at each branch
            notsupport_val<-list()
            #for each branch
            for (branch in br$Edge){
                #count number of ancestral alleles '0' at branch defining positions also ancestral '0'
                #count number of derived alleles '1' at branch defining positions also derived '1'
                notsupport_count<-sum(ancient$POS[ancient[ancient_sample]==1] %in% ancpos[[branch]] ) + sum(ancient$POS[ancient[ancient_sample]==0] %in% derpos[[branch]] )
                #add count of suportting positons at each branch for the current sample
                if(notsupport_count>0){
                    notsupport_val[[branch]]<-notsupport_count
                } else if (notsupport_count==0) {
                    notsupport_val[[branch]]<-0
                }
            }
        #convert to df
        support_count_at_edge<-data.frame(support_val)
        support_count_at_edge<-sapply(support_count_at_edge, as.numeric)
        names(support_count_at_edge)<-names(br$Edge)
        nonsupport_count_at_edge<-data.frame(notsupport_val)
        nonsupport_count_at_edge<-sapply(nonsupport_count_at_edge, as.numeric)
        names(nonsupport_count_at_edge)<-names(br$Edge)
        #make a table with branch, support count (count_at_edge), not support count (mismatchcount_at_edge)
        br_sample<-cbind(br,support=support_count_at_edge, notsupport=nonsupport_count_at_edge)
        br_sample<-br_sample[c("Edge","Node1","Node2","support","notsupport","hg")]
        #add that table to list of tables for all samples
        br_sample_tables[[samp]]<-br_sample
    }
}


plotAncDerSNPTree<-function(sample_name){
    samp<-make.names(sample_name)
    countdata<-br_sample_tables[[samp]]
    cat("SNP count: ",sum(countdata$support),"derived and", sum(countdata$notsupport), "ancestral.\n\n")
    pdf(file=paste0(args[4],"/",samp, '.counts_on_branches.pdf'))
    plot(tree, cex=0.1, edge.col="grey")
    edgelabels(edge=countdata$Edge[countdata$support>0],pch=20, col=alpha("darkgreen", 0.7),cex=log((countdata$support[countdata$support>0])+1)/2)
    edgelabels(edge=countdata$Edge[countdata$notsupport>0],pch=20, col=alpha("red", 0.5), cex=log((countdata$notsupport[countdata$notsupport>0])+1)/2)
    if (length(countdata$Edge[countdata$support>0])!=0){
    edgelabels(edge=countdata$Edge[countdata$support>0],text=countdata$support[countdata$support>0], frame = "none", cex=0.3)
    edgelabels(edge=countdata$Edge[countdata$notsupport>0],text=countdata$notsupport[countdata$notsupport>0], frame = "none", cex=0.3)
    }
    # if (length(countdata$hg[countdata$support>0])!=0){
    # edgelabels(edge=countdata$Edge[countdata$support>0],text=countdata$hg[countdata$support>0], frame = "none", cex=0.1, pos=1)
    # }
    dev.off()
}

invisible(sapply(samps, plotAncDerSNPTree))


writeSampleTables<-function(sample_name){
    samp<-make.names(sample_name)
    countdata<-br_sample_tables[[samp]]
    countdata$number_of_positions<-nums$number_of_positions[match(countdata$Edge, nums$Edge)]
    colnames(countdata)[which(colnames(countdata)=='notsupport')]<-'conflict'
    countdata<-countdata[c('Edge','Node1','Node2','support','conflict','number_of_positions','hg')]
    write.table(countdata, file=paste0(args[4],"/",samp, '.counts_on_branches.txt'),sep='\t', row.names=F, quote=F)
}

invisible(sapply(samps, writeSampleTables))

cat("output txt and pdf files written to", paste0(args[4],"/\n\n"))

saveRDS(br_sample_tables, file=paste0(args[4],"/allele_count_list.RData"))

if (length(excluded)>0){
    cat("Samples Excluded ", excluded,"\n")
}
