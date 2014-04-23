relabund <-read.table('coln.final.an.unique_list.0.03.subsample.relabund', header=T)
rownames(relabund) <- relabund[,2]
relabund <- relabund[,4:ncol(relabund)]
relabund <- t(relabund)
relabund <- relabund + 0.0001

coolPlot <- function(group1,group2){
  par(mar=c(5, 5, 0.5, 0.5))
  plot(1, type='n', xlim=c(1e-4,1), ylim=c(1e-4,1), ylab=group1, xlab=group2, xaxt="n", yaxt="n",log='xy')
  points(relabund[,group1], relabund[,group2])
  abline(0,1)
  axis(1, at=c(1e-4,1e-3,1e-2,1e-1,1), labels=c(0,0.001,0.01,0.1,1), las=1, tick=T, cex.axis=1)
  axis(2, at=c(1e-4,1e-3,1e-2,1e-1,1), labels=c(0,0.001,0.01,0.1,1), las=1, tick=T, cex.axis=1)
}