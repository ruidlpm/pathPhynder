require(phytools)
a<-read.tree("simulated_data.nwk")
a$tip.label<-as.numeric(a$tip.label)-1
a$tip.label<-paste("msp_",a$tip.label, sep='')
write.tree(a, file='simulated_data.nwk')
