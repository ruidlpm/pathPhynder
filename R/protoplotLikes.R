
##########################################################
# using logs

require(phytools)

ll<-list.files(pattern='posteriors.likes.txt')

for (i in ll){
	tree<-read.tree('~/BigTree/tree_1208forPaper_filt_allSites_210619-miss03-VAR.nwk')
	posterior_tree_order<-read.table(i)
	#b<-read.table('test.9.bam.posteriors.likes.txt')
	
	df<-data.frame(pp=b$x, col=NA)

	

	df_filt<-df[log(df$pp)>-5,]
	df_below<-df[log(df$pp)<(-5),]



	df_filt<-rbind(df_filt,data.frame(pp=1,col=NA),data.frame(pp=exp(-5),col=NA))


	df_below$col<-'grey'

	colPal <- colorRampPalette(c('grey','yellow','orange','red'))
	df_filt$col <- colPal(20)[as.numeric(cut(log(df_filt$pp),breaks = 20))]

	df_filt<-df_filt[1:(dim(df_filt)[1]-2),]

	df<-rbind(df_below,df_filt )

	df<-df[order(as.numeric(rownames(df))),]

	pdf(file=paste0(gsub('likes.txt','', i), 'post.pdf'),height=75, width=25)
	plot(ladderize(tree) ,  show.tip.label=T, edge.col=df$col, cex=0.25, col='grey', edge.width=1+df$pp)
	df_pphigh<-df_filt
	edgelabels(edge=as.numeric(rownames(df_pphigh)),col=df_pphigh$col,bg=df_pphigh$col,pch=21, cex=0.90)
	edgelabels(edge=as.numeric(rownames(df_pphigh)),text=round(log(df_pphigh$pp),3), frame='none', cex=0.2,bg=df_pphigh$col)

	best<-df_pphigh[which(log(df_pphigh$pp)==max(log(df_pphigh$pp))),]

	edgelabels(edge=as.numeric(rownames(best)),col=1,bg=best$col,pch=21, cex=0.90)

	edgelabels(edge=as.numeric(rownames(best)),text=round(log(best$pp),3), frame='none', cex=0.2,bg=best$col)

	
	dev.off()
}
##########################################################



