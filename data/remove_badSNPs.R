
isogg<-read.table('200803.snps_isogg.txt', h=F)
badsnps<-read.table('badsnps.txt', h=F)
mismatchsnps<-read.table('mismatchsnps.txt', h=F)

notes<-isogg[grep('Notes',isogg$V2, invert=F),]$V1


cat('\n')
cat('\n')


isogg_curated<-isogg[!(isogg$V1 %in% badsnps$V1),]
dim(isogg_curated)

isogg_curated<-isogg_curated[!(isogg_curated$V1 %in% mismatchsnps$V1),]
dim(isogg_curated)

isogg_curated<-isogg_curated[!(isogg_curated$V1 %in% notes),]
dim(isogg_curated)

isogg_curated<-isogg_curated[grep('maybe',isogg_curated$V1, invert=T),]
dim(isogg_curated)


isogg_curated<-isogg_curated[grep('Notes',isogg_curated$V2, invert=T),]



cat('Read' ,dim(isogg)[1],'raw ISOGG SNPs.' ,'\n')
cat('Excluded' ,dim(badsnps)[1],'bad ISOGG SNPs.' ,'\n')
cat('Excluded' ,dim(mismatchsnps)[1],'mismatching ISOGG SNPs.' ,'\n')
cat('Excluded' ,length(notes),' (Notes) ISOGG SNPs.' ,'\n')
cat('Wrote' ,dim(isogg_curated)[1],'curated ISOGG SNPs.' ,'\n\n')





write.table(isogg_curated, file='210513.snps_isogg_curated.txt', quote=F, sep='\t', row.names=F, col.names=F)



