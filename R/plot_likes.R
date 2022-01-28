######################################################
# plotting functions for phynder sample placement output
# (Rui Martiniano; Bianca de Sanctis 2019/2020)
######################################################

require(phytools, quietly = TRUE)
require(scales)

cat('\n\n',"Plotting phynder results", '\n\n\n')

args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=3) {
    stop("\tArguments needed.\n
    \tusage:
    \tRscript plot_likes.R <input_phylogeny.nwk> <query.phy> <results_folder> ", call.=FALSE)
} else {
    cat("   Command used:",'\n\n')
    cat(paste("plot_likes.R", args[1], args[2],args[3]), '\n\n')
}


######################################################
# 				     Functions
######################################################
plot_likes<-function(tree, branches, clade, suboptimal_branches,sample_name, results_folder, out_prefix){

  	cat(sample_name, '\n')

	df<-rbind(branches, suboptimal_branches)
	df$posterior<-as.numeric(as.character(df$posterior))
	df$score<-as.numeric(as.character(df$score))
	df$col<-NA
	df_filt<-df[log(df$posterior)>-6,]
	df_below<-df[log(df$posterior)<(-6),]

	colPal <- colorRampPalette(c('red','yellow',alpha("green", 0.7)))

    # if no branches above defined minumum threshold, plot clade only
	if (dim(df_filt)[1]==0){
		cat('No edges pass the minimum posterior probability of', exp(-6), '\n')

		pdf(file=paste0(results_folder,'/',out_prefix,'.posteriors.pdf'),height=estimatePlotDimensions(tree)[[1]], width=estimatePlotDimensions(tree)[[2]])
		par(mar = rep(2, 4))
		plot((tree) ,  show.tip.label=T, edge.col='grey', cex=0.25, col='grey', edge.width=1)

		edgelabels(edge=clade$clade,col=1,bg='blue',pch=21, cex=0.90)
		edgelabels(edge=clade$clade,text='C', frame='none', cex=0.2,bg='blue')

		cat(paste0('clade:', clade$clade), '\n\n')

		gradientLegend(valRange=c(0,1),pos=0.1,color =colPal(20), side=1, inside=F)

		dev.off()

	} else {
        #add branches with maximum and minimum score just to make colour scale, and then remove them

        df_filt<-rbind(df_filt,data.frame(type='N',branch='N',posterior=1, score='N',col=NA),data.frame(type='N',branch='N',posterior=exp(-6), score='N',col=NA))

		df_filt$col <- colPal(20)[as.numeric(cut((df_filt$posterior),breaks = 20))]

        df_filt<-(df_filt[df_filt$type!='N',])

		df_filt$branch<-as.numeric(as.character(df_filt$branch))
        df_filt<-df_filt[order(df_filt$branch),]

		to_add<-as.numeric(row.names(data.frame(tree$edge)))[which(!as.numeric(row.names(data.frame(tree$edge))) %in% df_filt$branch)]
		
        df_all<-rbind(df_filt,data.frame(type='N',branch=to_add,posterior=0, score='N',col='grey'))
        df_all<-df_all[order(df_all$branch),]



		pdf(file=paste0(results_folder,'/',out_prefix,'.posteriors.pdf'),height=estimatePlotDimensions(tree)[[1]], width=estimatePlotDimensions(tree)[[2]])
		par(mar = rep(2, 4))
		
        plot((tree) ,  show.tip.label=T, edge.col=df_all$col, cex=0.25, col='grey', edge.width=1)
	   
        df_all$branch<-as.numeric(as.character(df_all$branch))

        # print(df_all)

		edgelabels(edge=df_filt$branch,col=df_filt$col,bg=df_filt$col,pch=21, cex=0.90)
		edgelabels(edge=df_filt$branch,text=format(round(df_filt$posterior,3), nsmall = 3), frame='none', cex=0.2,bg=df_pphigh$col)

		best<-df_filt[df_filt$type=='B',]

		edgelabels(edge=best$branch,col=1,bg=best$col,pch=21, cex=0.90)
		edgelabels(edge=best$branch,text=format(round(best$posterior, 3), nsmall = 3), frame='none', cex=0.2,bg=best$col)


		edgelabels(edge=clade$clade,col=1,bg='blue',pch=21, cex=0.90)
		try(edgelabels(edge=clade$clade,text='C', frame='none', cex=0.2,bg='blue'))

        print(1)

		cat(paste0('branch:', best$branch), '\n')
		cat(paste0('clade:', clade$clade), '\n\n')

		gradientLegend(valRange=c(0,1),pos=0.1,color =colPal(20), side=1, inside=F)

		dev.off()

		}
}


tree.constructor = function(tree,edge.num,new.edge.length,pos, query.name){
  new.tree = bind.tip(tree,tip.label=paste('QUERY',query.name),
                      edge.length=new.edge.length,
                      where=edge.num, 
                      position = pos
                      )
  return(new.tree) 
}


estimatePlotDimensions<-function(tree){
  height=dim(tree$edge)[1]/75
  width=height*(2/3)
  #in case size too small, give a minimum size to prevent issues with plotting
  if (height<5 | width<5){
    height<-5
    width<-5
  }

  sizes<-list(height, width)
  return(sizes)
}


#taken from In plotfunctions: Various Functions to Facilitate Visualization of Data and Analysis
getFigCoords<-function (input = "f")
{
    p <- par()
    x.width = p$usr[2] - p$usr[1]
    y.width = p$usr[4] - p$usr[3]
    x.w = p$plt[2] - p$plt[1]
    y.w = p$plt[4] - p$plt[3]
    if (input == "f") {
        return(c(p$usr[1] - p$plt[1] * x.width/x.w, p$usr[2] +
            (1 - p$plt[2]) * x.width/x.w, p$usr[3] - p$plt[3] *
            y.width/y.w, p$usr[4] + (1 - p$plt[4]) * y.width/y.w))
    }
    else if (input == "p") {
        return(p$usr)
    }
    else if (input == "hp") {
        return(c(0.5 * x.width + p$usr[1], 0.5 * y.width + p$usr[3]))
    }
    else if (input == "hf") {
        return(c(p$usr[1] + (0.5 - p$plt[1]) * (x.width/x.w),
            p$usr[3] + (0.5 - p$plt[3]) * (y.width/y.w)))
    }
    else {
        return(NULL)
    }
}


getCoords<-function (pos = 1.1, side = 1, input = "p")
{
    p <- par()
    if (input == "p") {
        x.width = p$usr[2] - p$usr[1]
        y.width = p$usr[4] - p$usr[3]
        out <- rep(NA, length(pos))
        if (length(side) == 1) {
            side <- rep(side, length(pos))
        }
        out[which(side %in% c(1, 3))] <- pos[which(side %in%
            c(1, 3))] * x.width + p$usr[1]
        out[which(side %in% c(2, 4))] <- pos[which(side %in%
            c(2, 4))] * y.width + p$usr[3]
        return(out)
    }
    else if (input == "f") {
        gfc <- getFigCoords("f")
        x.width = gfc[2] - gfc[1]
        y.width = gfc[4] - gfc[3]
        out <- rep(NA, length(pos))
        if (length(side) == 1) {
            side <- rep(side, length(pos))
        }
        out[which(side %in% c(1, 3))] <- pos[which(side %in%
            c(1, 3))] * x.width + gfc[1]
        out[which(side %in% c(2, 4))] <- pos[which(side %in%
            c(2, 4))] * y.width + gfc[3]
        return(out)
    }
}


gradientLegend<-function (valRange, color = "topo", nCol = 30, pos = 0.5, side = 4,
    length = 0.25, depth = 0.05, inside = TRUE, coords = FALSE,
    pos.num = NULL, n.seg = 3, border.col = NULL, dec = NULL,
    fit.margin = TRUE)
{
    loc <- c(0, 0, 0, 0)
    if (is.null(pos.num)) {
        if (side %in% c(1, 3)) {
            pos.num = 3
        }
        else {
            pos.num = side
        }
    }
    if (length(pos) == 1) {
        pos.other <- ifelse(side > 2, 1, 0)
        if (side %in% c(1, 3)) {
            switch <- ifelse(inside, 0, 1)
            switch <- ifelse(side > 2, 1 - switch, switch)
            loc <- getCoords(c(pos - 0.5 * length, pos.other -
                switch * depth, pos + 0.5 * length, pos.other +
                (1 - switch) * depth), side = c(side, 2, side,
                2))
        }
        else if (side %in% c(2, 4)) {
            switch <- ifelse(inside, 0, 1)
            switch <- ifelse(side > 2, 1 - switch, switch)
            loc <- getCoords(c(pos.other - switch * depth, pos -
                0.5 * length, pos.other + (1 - switch) * depth,
                pos + 0.5 * length), side = c(1, side, 1, side))
        }
    }
    else if (length(pos) == 4) {
        if (coords) {
            loc <- pos
        }
        else {
            loc <- getCoords(pos, side = c(1, 2, 1, 2))
        }
    }
    mycolors <- c()
    if (length(color) > 1) {
        mycolors <- color
    }
    else if (!is.null(nCol)) {
        if (color == "topo") {
            mycolors <- topo.colors(nCol)
        }
        else if (color == "heat") {
            mycolors <- heat.colors(nCol)
        }
        else if (color == "terrain") {
            mycolors <- terrain.colors(nCol)
        }
        else if (color == "rainbow") {
            mycolors <- rainbow(nCol)
        }
        else {
            warning("Color %s not recognized. A palette of topo.colors is used instead.")
            mycolors <- topo.colors(nCol)
        }
    }
    else {
        stop("No color palette provided.")
    }
    vals <- seq(min(valRange), max(valRange), length = length(mycolors))
    if (!is.null(dec)) {
        vals <- round(vals, dec[1])
    }
    im <- as.raster(mycolors[matrix(1:length(mycolors), ncol = 1)])
    ticks <- c()
    if (side%%2 == 1) {
        rasterImage(t(im), loc[1], loc[2], loc[3], loc[4], col = mycolors,
            xpd = T)
        rect(loc[1], loc[2], loc[3], loc[4], border = border.col,
            xpd = T)
        ticks <- seq(loc[1], loc[3], length = n.seg)
        segments(x0 = ticks, x1 = ticks, y0 = rep(loc[2], n.seg),
            y1 = rep(loc[4], n.seg), col = border.col, xpd = TRUE)
    }
    else {
        rasterImage(rev(im), loc[1], loc[2], loc[3], loc[4],
            col = mycolors, xpd = T)
        rect(loc[1], loc[2], loc[3], loc[4], border = border.col,
            xpd = T)
        ticks <- seq(loc[2], loc[4], length = n.seg)
        segments(x0 = rep(loc[1], n.seg), x1 = rep(loc[3], n.seg),
            y0 = ticks, y1 = ticks, col = border.col, xpd = TRUE)
    }
    determineDec <- function(x) {
        out = max(unlist(lapply(strsplit(x, split = "\\."), function(y) {
            return(ifelse(length(y) > 1, nchar(gsub("^([^0]*)([0]+)$",
                "\\1", as.character(y[2]))), 0))
        })))
        return(out)
    }
    labels = sprintf("%f", seq(min(valRange), max(valRange),
        length = n.seg))
    if (is.null(dec)) {
        dec <- min(c(6, determineDec(labels)))
    }
    eval(parse(text = sprintf("labels = sprintf('%s', round(seq(min(valRange), max(valRange), length = n.seg), dec) )",
        paste("%.", dec, "f", sep = ""))))
    if (pos.num == 1) {
        if (fit.margin) {
            lab.height = max(strheight(labels)) * 0.8
            max.pos = getFigCoords()[3]
            if ((max.pos - loc[2]) < lab.height) {
                warning("Increase bottom margin, because labels for legend do not fit.")
            }
        }
        text(y = loc[2], x = ticks, labels = seq(min(valRange),
            max(valRange), length = n.seg), col = border.col,
            pos = 1, cex = 0.8, xpd = T)
    }
    else if (pos.num == 2) {
        if (fit.margin) {
            checkagain = TRUE
            while (checkagain == TRUE) {
                lab.width = (max(strwidth(labels)) + 0.5 * par()$cxy[1]) *
                  0.8
                min.pos = getFigCoords()[1]
                if ((loc[1] - min.pos) < lab.width) {
                  if (!is.null(dec)) {
                    dec = max(c(0, dec - 1))
                    if (dec == 0) {
                      warning("Decimal set to 0 (dec=0), but the labels still don't fit in the margin. You may want to add the color legend to another side, or increase the margin of the plot.")
                      checkagain = FALSE
                    }
                  }
                  else {
                    tmp = max(unlist(lapply(strsplit(labels,
                      split = "\\."), function(x) {
                      return(ifelse(length(x) > 1, nchar(x[2]),
                        0))
                    })))
                    dec = max(c(0, tmp - 1))
                    if (dec == 0) {
                      warning("Decimal set to 0 (dec=0), but the labels still don't fit in the margin. You may want to add the color legend to another side, or increase the margin of the plot.")
                      checkagain = FALSE
                    }
                  }
                  eval(parse(text = sprintf("labels = sprintf('%s', round(seq(min(valRange), max(valRange), length = n.seg), dec) )",
                    paste("%.", dec, "f", sep = ""))))
                }
                else {
                  checkagain = FALSE
                }
            }
        }
        text(y = ticks, x = loc[1], labels = labels, pos = 2,
            cex = 0.8, col = border.col, xpd = T)
    }
    else if (pos.num == 3) {
        if (fit.margin) {
            lab.height = max(strheight(labels)) * 0.8
            max.pos = getFigCoords()[4]
            if ((max.pos - loc[4]) < lab.height) {
                warning("Increase top margin, because labels for legend do not fit.")
            }
        }
        text(y = loc[4], x = ticks, labels = seq(min(valRange),
            max(valRange), length = n.seg), col = border.col,
            pos = 3, cex = 0.8, xpd = T)
    }
    else if (pos.num == 4) {
        if (fit.margin) {
            checkagain = TRUE
            while (checkagain == TRUE) {
                lab.width = (max(strwidth(labels)) + 0.5 * par()$cxy[1]) *
                  0.8
                max.pos = getFigCoords()[2]
                if ((max.pos - loc[3]) < lab.width) {
                  if (!is.null(dec)) {
                    dec = max(c(0, dec - 1))
                    if (dec == 0) {
                      warning("Decimal set to 0 (dec=0), but the labels still don't fit in the margin. You may want to add the color legend to another side, or increase the margin of the plot.")
                      checkagain = FALSE
                    }
                  }
                  else {
                    tmp = max(unlist(lapply(strsplit(labels,
                      split = "\\."), function(x) {
                      return(ifelse(length(x) > 1, nchar(x[2]),
                        0))
                    })))
                    dec = max(c(0, tmp - 1))
                    if (dec == 0) {
                      warning("Decimal set to 0 (dec=0), but the labels still don't fit in the margin. You may want to add the color legend to another side, or increase the margin of the plot.")
                      checkagain = FALSE
                    }
                  }
                  eval(parse(text = sprintf("labels = sprintf('%s', round(seq(min(valRange), max(valRange), length = n.seg), dec) )",
                    paste("%.", dec, "f", sep = ""))))
                }
                else {
                  checkagain = FALSE
                }
            }
        }
        text(y = ticks, x = loc[3], labels = labels, pos = 4,
            col = border.col, cex = 0.8, xpd = T)
    }
}

get_relevant_line <- function(phy_table){
    line_counter<-0
    for (i in readLines('testquery.phy')){
        line_counter=line_counter+1
        if(startsWith(i, 'I')){
            return(line_counter)
        }
    }
}

get_branches<-function(sample, query_results, query_results_index){
    tmp<-query_results_index[which(query_results_index$samples==sample,),]
    tmp2<-query_results[tmp$line:tmp$last_line,]
    return(tmp2)
}


######################################################
#                    Main
######################################################
tree<-read.tree(args[1])
results_folder<-args[3]


query_results<-read.table(args[2], fill=T, stringsAsFactors=F)
query_results<-query_results[1:4]


dir.create(results_folder, showWarnings = FALSE)



samples<-as.character(query_results[query_results$V1=='I',]$V3)


df<-data.frame(samples=as.character(query_results[query_results$V1=='I',]$V3), line=which(query_results$V1=='I'), last_line=c( which(query_results$V1=='I')[2:length(which(query_results$V1=='I'))]-1,dim(query_results)[1]))

new.tree<-tree
for(sample_name in samples){

    res<-get_branches(sample_name,query_results,df)
    branches<- res[res$V1=='B',]
    colnames(branches)<-c('type','branch', 'posterior', 'score')
    clade<- res[res$V1=='C',]
    colnames(clade)<-c('type','clade')
    suboptimal_branches<- res[res$V1=='S',]
    colnames(suboptimal_branches)<-c('type','branch', 'posterior', 'score')
    

    branches$branch<-as.numeric(as.character(branches$branch))
    suboptimal_branches$branch<-as.numeric(as.character(suboptimal_branches$branch))
    clade$clade<-as.numeric(as.character(clade$clade))
    plot_likes(tree, branches, clade,suboptimal_branches,sample_name, results_folder,sample_name )


    descendants<-tree$tip.label[getDescendants(tree, tree$edge[branches$branch,][2])]
    descendants<- descendants[descendants %in% tree$tip.label]

    if(length(descendants)>1){
        updated_branch<-which(new.tree$edge[,2]==getMRCA(new.tree,which(new.tree$tip.label  %in% descendants)))
        updated_node<-new.tree$edge[updated_branch,][2]
        new.tree<-tree.constructor(new.tree,updated_node,1e-3 * min(new.tree$edge.length),new.tree$edge.length[updated_branch]/2,sample_name)
    } else {
        updated_branch<-which(new.tree$edge[,2]==which(new.tree$tip.label  %in% descendants))
        updated_node<-new.tree$edge[updated_branch,][2]
        new.tree<-tree.constructor(new.tree,updated_node,1e-3 * min(new.tree$edge.length),new.tree$edge.length[updated_branch]/2,sample_name)
    }      
}


write.tree(new.tree,file=paste0(results_folder,'/final_tree.phynder.nwk'))

plot_final_tree<-function(new.tree, results_folder){
    pdf(file=paste0(results_folder,'/final_tree.phynder.pdf'),height=estimatePlotDimensions(new.tree)[[1]], width=estimatePlotDimensions(new.tree)[[2]])
    plot((new.tree) ,  show.tip.label=T,tip.color=ifelse(1:length(new.tree$tip.label) %in% grep("QUERY", new.tree$tip.label), yes=2, no=1), cex=0.15)
    dev.off()
}


plot_final_tree(new.tree,results_folder)


