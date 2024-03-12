library(tidyverse)
setwd("Desktop/Sheila/file/plink/")

indmiss = read.table(file="plink.imiss", header=TRUE)
snpmiss = read.table(file="plink.lmiss", header=TRUE)
# read data into R
summary(indmiss)
summary(snpmiss)

hist(indmiss[,6],main="Histogram individual missingness") #selects column
hist(snpmiss[,5],main="Histogram SNP missingness")

# pdf("histimiss.pdf") #indicates pdf format and gives title to file
hist(indmiss[,6],main="Histogram individual missingness") #selects column
dev.off()
# pdf("histlmiss.pdf")
hist(snpmiss[,5],main="Histogram SNP missingness")
dev.off() # shuts down the current device

maf_freq = read.table("maf.frq", header =TRUE, as.is=T)
# pdf("MAF_distribution.pdf")
hist(maf_freq[,5],main = "MAF distribution", xlab = "MAF")
dev.off()

hwe = read.table (file="plink.hwe", header=TRUE)
hist(hwe[,9],main="Histogram HWE")
# pdf("histhwe.pdf")

het = read.table("R_check.het", head=TRUE)
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
hist(het$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main=
       "Heterozygosity Rate")
# pdf("heterozygosity.pdf")
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
hist(het$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main=
       "Heterozygosity Rate")
dev.off()

het = read.table("R_check.het", head=TRUE)
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
het_fail = subset(het, (het$HET_RATE <
                          mean(het$HET_RATE)-3*sd(het$HET_RATE)) | (het$HET_RATE >
                                                                      mean(het$HET_RATE)+3*sd(het$HET_RATE)));

het_fail$HET_DST = (het_fail$HET_RATE-
                      mean(het$HET_RATE))/sd(het$HET_RATE);

#write.table(het_fail, "fail-het-qc.txt", row.names=FALSE)

relatedness = read.table("pihat_min0.2.genome", header=T)
with(relatedness,plot(Z0,Z1, xlim=c(0,1), ylim=c(0,1), type="n"))
with(subset(relatedness,RT ="PO") , points(Z0,Z1,col=4))
with(subset(relatedness,RT ="UN") , points(Z0,Z1,col=3))
# legend(1,1, xjust=1, yjust=1, legend=levels(relatedness$RT), pch=16, col=c(4,3))

relatedness_zoom = read.table("zoom_pihat.genome", header=T)
with(relatedness_zoom,plot(Z0,Z1, xlim=c(0,0.02), ylim=c(0.98,1),
                           type="n"))
with(subset(relatedness_zoom,RT ="PO") , points(Z0,Z1,col=4))
with(subset(relatedness_zoom,RT ="UN") , points(Z0,Z1,col=3))
# legend(0.02,1, xjust=1, yjust=1, legend=levels(relatedness$RT), pch=16, col=c(4,3))

#pdf("hist_relatedness.pdf")
relatedness = read.table("pihat_min0.2.genome", header=T)
hist(relatedness[,10],main="Histogram relatedness", xlab= "Pihat")
dev.off()

### END of QC ###

vcf = read.table("result.df.tsv", header=TRUE, comment = "")
vcf = vcf %>% filter(REF != 'I' & REF != 'D' & ALT != 'I' & ALT != 'D')
write.table(vcf, file = "filtered.vcf.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")




