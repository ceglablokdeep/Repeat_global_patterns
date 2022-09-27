###In R###
data=read.table("gcpercompiled_allclades")
##variety treatment note
colnames(data)=c("clades","status","gc")
new_order <- with(data, reorder(clades , gc, mean , na.rm=T))
jpeg("gc_percent_cladewise_repeats_nonrepeats_comparison.jpeg",width=28,height=18,units="in",res=300)
par(mar=c(5,4,4,1))
myplot <- boxplot(gc ~ status*new_order , data=data  , 
        boxwex=0.4,
        ann=F,
        col=c("seagreen4" , "coral3") , 
        border=c("seagreen3","coral2"),cex.axis=1.25,
        xaxt="n")
title(main="GC% variation in proteins with repeats and proteins without repeats" ,line=1,cex.main=2.3)
title(ylab="GC%", line=2, cex.lab=2.5)        
title(xlab="Clades", line=3, cex.lab=2.5)
my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
my_names <- my_names[seq(1 , length(my_names) , 2)]
##change the number in seq below according to your number of clades. IT should be 2*(number of clades)
axis(1, 
     at = seq(1.5 , 24 , 2), 
     labels = my_names , 
     tick=FALSE , cex.axis=1.75)
     
# Add the grey vertical lines
##change the seq range (in our case it is 30) below according to your number of clades 
for(i in seq(0.5 , 30 , 2)){ 
  abline(v=i,lty=3, col="grey",lwd=1.75)
  }
 
# Add a legend ##Adjust the position of legend by changing x and y below
legend(x=21,y=34, legend = c("Proteins with repeats", "Proteins without repeats"), 
       col=c("seagreen4" , "coral3"),
       pch = 15, pt.cex = 6, cex = 2,  horiz = F, inset = c(0.1, 0.1),box.lwd = 0,box.col = "white",bg = "white")

###adding the significant values for each clade paired comparison
##in yax, I am adding +4 to keep it slightly above the 3rd quantile range of boxplot, you can adjust it according to your plot
xpos=1
for(clade in unique(my_names)){
a1=data[data$clades==clade & data$status=="With_repeats",]
b1=data[data$clades==clade & data$status=="Without_repeats",]
yax=max(quantile(a1$gc,probs=0.75),quantile(b1$gc,probs=0.75))+4
pval=signif(as.numeric(wilcox.test(a1$gc,b1$gc,alternative="greater")[3]),digits=3)
segments(xpos,yax,xpos+1,yax,lty=3,col="purple")
text(xpos+0.5,yax+0.6,labels=pval,cex=1.5)
a1mean=round(mean(a1$gc),digits=2)
b1mean=round(mean(b1$gc),digits=2)
yminax=min(quantile(a1$gc,probs=0.25),quantile(b1$gc,probs=0.25))-2
text(xpos,yminax,bquote(bar("x") == .(a1mean)),cex=1.1)
text(xpos+1,yminax,bquote(bar("x") == .(b1mean)),cex=1.1)
text(xpos,yminax-3,bquote("n" == .(dim(a1)[1])),cex=1.1)
text(xpos+1,yminax-3,bquote("n" == .(dim(b1)[1])),cex=1.1)
xpos=xpos+2
}
dev.off()

###########This script ends here##################
