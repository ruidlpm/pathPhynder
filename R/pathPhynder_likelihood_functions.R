# nov 14 2019
# felsenstein likelihood calculator and edge placement optimizer for pathphynder
# bianca de sanctis

# added plotting functions (Rui 22/11/2019)

tree.constructor = function(tree,edge.num,new.edge.length,pos){
  new.tree = bind.tip(tree,tip.label=query.name,
                      edge.length=new.edge.length,
                      where=tree$edge[edge.num,2], 
                      position = pos
                      )
  return(new.tree) 
}

force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
} # this function taken from http://blog.phytools.org/2017/03/forceultrametric-method-for-ultrametric.html 

calculate.prior <- function(tree,type="coalescent"){
  # if the prior is a uniform
  num.edge = length(tree$edge.length)
  
  if(type == "uniform"){
    prior.dist = rep(1/num.edge,num.edge)
  } else if(type == "coalescent"){
    # note: some of these commands are in a very specific order, don't mess with it
    rescale.tree = force.ultrametric(tree,method="extend")  # the other method messes up the edge labels.
    n = coalescent.intervals(tree)$lineages[1]
    expected.tmrca = 2*(1-1/n)
    actual.tmrca = coalescent.intervals(rescale.tree)$total.depth
    rescale.tree$edge.length = rescale.tree$edge.length*(expected.tmrca/actual.tmrca)
    btimes = branching.times(rescale.tree) 
    coal.times = coalescent.intervals(rescale.tree)$interval.length
    skyl = cumsum(coal.times) # = skyline(rescale.tree)$time
    how.many.lineages = coalescent.intervals(rescale.tree)$lineages
    
    prior.dist = rep(0,num.edge)
    for(edgenum in 1:num.edge){
      topnode = rescale.tree$edge[edgenum,1]; bottomnode = rescale.tree$edge[edgenum,2]
      top.to.tip = unname(btimes[as.character(topnode)]) # length from top node to tips
      if(length(Children(rescale.tree,bottomnode)) == 0){bottom.to.tip = 0}
      else{bottom.to.tip = unname(btimes[as.character(bottomnode)])} # length from bottom nodes to tips
      highest.coal.num = which(abs(skyl-top.to.tip) == min(abs(skyl-top.to.tip)))[1] # it might not be exactly equal because of decimal expansions, so this finds the closest number
      lowest.coal.num = which(abs(skyl-bottom.to.tip) == min(abs(skyl-bottom.to.tip)))[1]
      if(bottom.to.tip ==0){lowest.coal.num=0}
      # so we have to add up all of the coalescence times between these two
      # we also want to index somehow so that skyl[0] = 0
      sum=0
      if(highest.coal.num == lowest.coal.num){print(edgenum)}
      for(coalnum in (lowest.coal.num+1):highest.coal.num){
        if(coalnum==1){sum = sum + (1-exp(-skyl[1]))*1/how.many.lineages[1]}
        else{sum = sum + (exp(-skyl[coalnum-1]) - exp(-skyl[coalnum]))*1/how.many.lineages[coalnum]}
        # ^ this is a ridiculous way of ensuring that the cumulative function works
        # recall that exp(-t) is the probability of coalescing before or at time t (it's a cumulative distribution function)
        # so in the case where the bottom cumulant is 0, we need to not subtract anything
      }
      prior.dist[edgenum] = sum
    }
    prior.dist = prior.dist / sum(prior.dist) # now normalize.
    # notice: sum(prior.dist) before normalization should be (more or less) equal to expected.tmrca.
  } else {
    stop('The --prior parameter must be \'coalescent\' or \'uniform\'.')
  }
  prior.dist
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



require(scales)

plot_likes<-function(tree, likes_table, sample_name, results_folder, prior){
  df<-data.frame(pp=posterior.normal.order, col=NA)

  df_filt<-df[log(df$pp)>-6,]
  df_below<-df[log(df$pp)<(-6),]

  df<-data.frame(pp=posterior.normal.order, col=NA)

  df_filt<-df[log(df$pp)>-6,]
  df_below<-df[log(df$pp)<(-6),]
  df_filt<-rbind(df_filt,data.frame(pp=1,col=NA),data.frame(pp=exp(-6),col=NA))
  if (dim(df_below)[1]>0){
      df_below$col<-'grey'
    }
  # colPal <- colorRampPalette(c('grey','yellow','orange','red'))
  colPal <- colorRampPalette(c('red','yellow',alpha("green", 0.7)))
  # df_filt$col <- colPal(20)[as.numeric(cut(log(df_filt$pp),breaks = 20))]
  df_filt$col <- colPal(20)[as.numeric(cut((df_filt$pp),breaks = 20))]
  df_filt<-df_filt[1:(dim(df_filt)[1]-2),]
  df<-rbind(df_below,df_filt )
  df<-df[order(as.numeric(rownames(df))),]
  

  pdf(file=paste0(results_folder,'/',out_prefix,'.posteriors.',prior,'.pdf'),height=estimatePlotDimensions(tree)[[1]], width=estimatePlotDimensions(tree)[[2]])
  par(mar = rep(2, 4))
  plot((tree) ,  show.tip.label=T, edge.col=df$col, cex=0.25, col='grey', edge.width=1+df$pp)
  df_pphigh<-df_filt
  edgelabels(edge=as.numeric(rownames(df_pphigh)),col=df_pphigh$col,bg=df_pphigh$col,pch=21, cex=0.90)
  edgelabels(edge=as.numeric(rownames(df_pphigh)),text=round((df_pphigh$pp),2), frame='none', cex=0.2,bg=df_pphigh$col)

  best<-df_pphigh[which((df_pphigh$pp)==max((df_pphigh$pp))),]

  edgelabels(edge=as.numeric(rownames(best)),col=1,bg=best$col,pch=21, cex=0.90)
  edgelabels(edge=as.numeric(rownames(best)),text=round((best$pp),2), frame='none', cex=0.2,bg=best$col)
  
  gradientLegend(valRange=c(0,1),pos=0.1,color =colPal(20), side=1, inside=F)

  dev.off()

}



plot_final_tree<-function(new.tree, out_prefix, results_folder, prior){
  pdf(file=paste0(results_folder,'/',out_prefix,'final_tree.',prior,'.pdf'),height=estimatePlotDimensions(tree)[[1]], width=estimatePlotDimensions(tree)[[2]])
  plot(ladderize(new.tree) ,  show.tip.label=T)
  dev.off()
}




require(scales)

# plot_likes<-function(tree, likes_table, sample_name, results_folder, prior, pos){
#     new.tree_noLen<-tree
#   new.tree_noLen$edge.length<-NULL
#   df<-data.frame(pp=posterior.normal.order, col=NA)

#   df_filt<-df[log(df$pp)>-6,]
#   df_below<-df[log(df$pp)<(-6),]

#   df<-data.frame(pp=posterior.normal.order, col=NA)

#   df_filt<-df[log(df$pp)>-6,]
#   df_below<-df[log(df$pp)<(-6),]
#   df_filt<-rbind(df_filt,data.frame(pp=1,col=NA),data.frame(pp=exp(-6),col=NA))
#   if (dim(df_below)[1]>0){
#       df_below$col<-'grey'
#     }
#   # colPal <- colorRampPalette(c('grey','yellow','orange','red'))
#   colPal <- colorRampPalette(c('red','yellow',alpha("green", 0.7)))
#   # df_filt$col <- colPal(20)[as.numeric(cut(log(df_filt$pp),breaks = 20))]
#   df_filt$col <- colPal(20)[as.numeric(cut((df_filt$pp),breaks = 20))]
#   df_filt<-df_filt[1:(dim(df_filt)[1]-2),]
#   df<-rbind(df_below,df_filt )
#   df<-df[order(as.numeric(rownames(df))),]
  

#   pdf(file=paste0(results_folder,'/',out_prefix,'.posteriors.',prior,'.pdf'),height=estimatePlotDimensions(tree)[[1]], width=estimatePlotDimensions(tree)[[2]])
#   par(mar = rep(2, 4))
#   plot((new.tree_noLen) ,  show.tip.label=T, edge.col=df$col, cex=0.8, col='lightgrey', edge.width=1, tip.color = 'grey')
#   df_pphigh<-df_filt
#   edgelabels(edge=as.numeric(rownames(df_pphigh)),col=df_pphigh$col,bg=df_pphigh$col,pch=21, cex=3)
#   edgelabels(edge=as.numeric(rownames(df_pphigh)),text=round((df_pphigh$pp),2), frame='none', cex=0.75,bg=df_pphigh$col)

#   best<-df_pphigh[which((df_pphigh$pp)==max((df_pphigh$pp))),]

#   edgelabels(edge=as.numeric(rownames(best)),col=1,bg=best$col,pch=21, cex=3)
#   edgelabels(edge=as.numeric(rownames(best)),text=round((best$pp),2), frame='none', cex=0.75,bg=best$col)
  
#   gradientLegend(valRange=c(0,1),pos=0.1,color =colPal(20), side=1, inside=F)

#   dev.off()

# }





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



