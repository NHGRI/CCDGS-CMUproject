manhattanCaseControl <- function(ds1=ds1,ds2=ds2,pvalCthresh=.01)
{
  library(karyoploteR)
  ymax = -log10(min(min(ds1$pval),min(ds2$pval)))
  ymax = 30
  print('Limist set for DC 0.6 max, remove line 5 to make it dynamic')
  kp <- plotKaryotype(plot.type=4,chromosomes = 'autosomal',labels.plotter = NULL)
  kp <- kpAddChromosomeNames(kp,chr.names = as.character(c(1:22)))
  kpAddLabels(kp, labels = "", srt=90, pos=3, r0=0.5, r1=1, cex=1.8, label.margin = 0.025)
  kpAxis(kp, ymin=0, ymax=ymax, r0=0.5,tick.len = 0)
  pvalsig <- max(max(ds1$pval[p.adjust(ds1$pval,"bonferroni")<pvalCthresh]),0)
  kpDataBackground(kp,color = "mistyrose",r0=0.5,r1=1)
  t <- try(kp <- kpPlotManhattan(kp, data=ds1, highlight = ds1[ds1$grtVal],r0=0.5, r1=1, ymax=ymax,highlight.col = "red",suggestiveline = -log10(pvalsig),genomewideline = -log10(pvalsig)))
  if ("try-error" %in% class(t)) kp <- kpPlotManhattan(kp, data=ds1,r0=0.5, r1=1, ymax=ymax,suggestiveline = -log10(pvalsig),genomewideline = -log10(pvalsig))
  kpAddLabels(kp, labels = "", srt=90, pos=3, r0=0, r1=0.5, cex=1.8, label.margin = 0.025)
  kpAxis(kp, ymin=0, ymax=ymax, r0=0.5, r1=0,tick.len = 0)
  pvalsig <- max(max(ds2$pval[p.adjust(ds2$pval,"bonferroni")<pvalCthresh]),0)
  kpDataBackground(kp,color = "darkseagreen1",r0=0.5,r1=0)
  t <- try(kp <- kpPlotManhattan(kp, data=ds2, highlight = ds2[ds2$grtVal],highlight.col = "green",r0=0.5, r1=0, ymax=ymax, points.col = "2blues",suggestiveline = -log10(pvalsig),genomewideline = -log10(pvalsig)))
  if ("try-error" %in% class(t)) kp <- kpPlotManhattan(kp, data=ds2,r0=0.5, r1=0, ymax=ymax, points.col = "2blues",suggestiveline = -log10(pvalsig),genomewideline = -log10(pvalsig))
  
}
