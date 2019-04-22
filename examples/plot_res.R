require("phytools")
tree<-read.tree("small_example_tree")

pdf(file="small_example_tree.pdf")
res<-read.table("tree_data/small_example_tree.edge_df.txt", h=T, sep="\t")
res<-res[!is.na(res$positions),]
plot(tree, direction="down", edge.color="grey",  edge.width=3, tip.color="darkgrey")
nodelabels(node=res$Node2, frame='none', pch=15, col="steelblue")
nodelabels(node=res$Node2, pos=4, frame='none', text=res$Node2, col="steelblue")
edgelabels(edge=res$Edge, frame='none', pch=20, col="indianred")
edgelabels(edge=res$Edge, frame='none',text=res$Edge, pos=2, col="indianred")
dev.off()

