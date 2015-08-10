## Differential Expression calls using edgeR
# Based off edgeR lab with input from Upendra, Susan, Julin, Kazu
#This particular script drops low count libraries, 
##also drops libraries with incorrect genotyping as seen in IGV (accessions mixded up)
## subsetted for Col, Sha (no mutants in Col background), sun only, all time points
# 2015-07-30
# Palmer


library(edgeR)

##read in file
allseq <- read.table("SAS_counts_merged_updated.bam.tsv", header = T)
names(allseq)
head(allseq)[1:10]

##turn NA's into 0 (since these represent no expression)
allseq[is.na(allseq)] <- 0
head(allseq)

##get rid of top row with * which is total counts
allseq<-allseq[-1,]

##trim to just samples that have at least 1 million reads...eliminates small libraries
allseq <- allseq[,colSums(allseq[,-1])>1000000]
head(allseq)[1:10]

##need to do this b/c right now the gene names are part of the data frame and aren't counts...make them rownames instead
rownames(allseq)<-allseq$gene

##and now get rid of that column
allseq<-allseq[,-1]

##SUBSET ALLSEQ!!!   This step will remove the libraries that show likely incorrect genotypes
# as based on MDS plot and IGV.  Accessions are messed up.  This does NOT address mistakes in Kazu
# mutant genotypes.  Simply Col vs other accessions mistakes.
# bad libraries are D9L4,C14H25,D16L49,D18L1,C22H25A,E22L16A,C20H1,C20L1,E25H25,D24H16A 

allseqsub<-allseq[,colnames(allseq)!= "D9L4A.merged.bam" & colnames(allseq)!="C14H25A.merged.bam" &
                    colnames(allseq)!= "D16L49A.merged.bam" & colnames(allseq)!= "D18L1A.merged.bam" &
                    colnames(allseq)!= "C22H25A.merged.bam" & colnames(allseq)!= "E22L16A.merged.bam" &
                    colnames(allseq)!= "C20H1A.merged.bam" & colnames(allseq)!= "C20L1A.merged.bam" &
                    colnames(allseq)!= "E25H25A.merged.bam" & colnames(allseq)!= "D24H16A.merged.bam"]
dim(allseqsub)  

## Now subset out just the sun samples (H) for col (23) and sha (22) for all time points

allseqsub <- allseqsub[,grep("(C|D|E)(22|23)(H)(1|4|16|25|49)(A.merged.bam)", colnames(allseqsub))]
dim(allseqsub)

## Dataframe
samples <- data.frame(
  file = names(allseqsub),
  gt   = sub("(C|D|E)(22|23)(H)+(1|4|16|25|49)+A+[[:print:]]+(.bam)","\\2",names(allseqsub)),
  trt  = sub("(C|D|E)(22|23)(H)+(1|4|16|25|49)+A+[[:print:]]+(.bam)","\\3",names(allseqsub)),
  rep  = sub("(C|D|E)(22|23)(H)+(1|4|16|25|49)+A+[[:print:]]+(.bam)","\\1",names(allseqsub)),
  time = sub("(C|D|E)(22|23)(H)+(1|4|16|25|49)+A+[[:print:]]+(.bam)","\\4",names(allseqsub))
)
samples

## prepping for design
gt <- relevel(samples$gt, ref="23")
trt <- relevel(samples$trt, ref = "H")  #not needed for this, but leaving it in for future use on sun/shade
rep <- samples$rep
time <- relevel(samples$time, ref = "1")

# preparing grouping here
group <- factor(paste(samples$gt,samples$trt,sep = "_"))
group
group <- relevel(group, ref="23_H")

# Adding group to the dataframe
samples$group <- group
samples

## this specifies counts (from "data", and the groups here are the different genotypes, "group")
dge <- DGEList(counts = allseqsub, group = group)
dim(dge)

#only keep genes where at least 1 sample has more than 3 counts per million
dge <- dge[rowSums(cpm(dge) > 1) >= 3,] 
dim(dge)

# NORMALIZATION HERE!
dge <- calcNormFactors(dge, method = "TMM") #TMM is the default normalization method, but you could specify a different one
names(dge)
allseqsub[row.names(allseqsub)=="AT4G16780.1",]  ## this is atHb2, again, not needed here, but useful for if you do shade in the future

# cpm (counts per million aka... normalized counts)
countspmill <- cpm(dge, normalized.lib.sizes=TRUE)
head(countspmill)
grep("AT4G16780.1", row.names(countspmill)) ###this is atHb2
countspmill[10466,] # 

# MDS plot
MDS.plot <- plotMDS(dge)
plotMDS(dge, ylim = c(-2.5, 2.5), xlim = c(-7, 10))#, labels=groupall) 
##good grouping by genotype (22 vs 23) with 2 exceptions...
##D22H49A  and D23H16A  look like strong outliers and cluster on top of each other...might want to get rid of them. 
##-->>but when I get rid of them, others seem to get messed up and no longer cluster by genotype??!!
#E23H49 and E22H49 are okay....they don't sit on top of each other, but they aren't in the nice cluster of 22's or 23's
#otherwise, everything clusters nicely.


# calculate dispersion (which is the variation of each gene expression)

##this essentially turns group information into dataframe
Design <- model.matrix(~ gt + time + rep, data = allseqsub)
Design[1:6, ] 

##common assumes all of the genes have the same variance (unlikely! but ideal.)
dge <- estimateGLMCommonDisp(dge, Design) #, verbose = T)  

##trended groups things based on gene expression and calculates variance by "bin"
dge <- estimateGLMTrendedDisp(dge, Design) #, verbose = T) # outputs the many trended dispersion values, for gene sets binned by expression level

##tagwise is assuming each gene has a different variation (most realistic, but most painful....uses info from trended and pushes it to common...most useful to be able to think about everything having the same variance.  Not clear exactly how this works, but ... use this one??  yes, usually, unless you have toooo many factors to do it computationally,  so for you, yeah.)
dge <- estimateGLMTagwiseDisp(dge, Design) # can't do verbose here, it will print out a value for all 14000 genes!

plotBCV(dge) 


## Identification of Differentially Expressed Genes 

## glmFit (instead of exact test)
fit<-glmFit(dge, Design)
colnames(Design)
head(Design[,2])
head(Design[,c(2,31)])
head(Design[,28:30])
lrt<-glmLRT(fit, coef=2)
names(lrt)

##get table with results
toplrt<-topTags(lrt, adjust.method="BH", n = Inf)$table
toplrt<- toplrt[toplrt$FDR<=0.01,]
dim(toplrt)
toplrt$gene<-rownames(toplrt)
toplrt$gene<-substring(as.character(toplrt$gene), 1, 9)
#tells us the number of genes in the model that are sig at FDR<0.01
FDRlrt <- p.adjust(lrt$table$PValue, method="BH")     
length(FDRlrt) # total number of genes
sum(FDRlrt < 0.01) # number of DE genes at FDR < 0.01

sigs<-lrt[lrt$table$PValue < 0.01,]

# annotate the genes
annot <- read.delim("TAIR10_functional_descriptions.txt", quote = "", header = T)

# next we need to chop off the extra characters on the gene names in annot
annot$name <- substring(as.character(annot$name), 1, 9)
names(annot) # each column gives a similar but different description of the gene

# combine the datasets based on the ATG gene numbers; annotate with the short description
toplrt.annot <- merge(x = toplrt, y = annot[, c(1, 3)], by.x = "gene", by.y = 1, all.x = T, all.y = F) 
names(toplrt.annot )

#write out a file of DE genes: 
top4x <- topTags(lrt, n = Inf, adjust.method = "BH")$table

# cut off at at least 4-fold change btw col and sha
top4x <- top4x[abs(top4x$logFC) >= 4 & top4x$FDR <= 0.01,  ] 
dim(top4x) # how many genes are in this list?


rownames(top4x) <- substring(rownames(top4x),1,9)

top4x.annot <- merge(x = top4x, y = annot[, c(1, 3)], by.x = 0, by.y = 1, all.x = T, all.y = F) # combine the datasets based on the ATG gene numbers; annotate with the short description
top4x.annot

write.csv(top4x.annot,"~/Desktop/colshaSUNalltimeEdgeRminusbad4x.csv")

#write out a file of DE genes: 
top2x <- topTags(lrt, n = Inf, adjust.method = "BH")$table
top2x <- top2x[abs(top2x$logFC) >= 2 & top2x$FDR <= 0.01,  ] # cut off at at least 4-fold change btw col and sha
dim(top2x) # how many genes are in this list?
#791 genes with model of ~gt + rep + time,  coefficient for gt

rownames(top2x) <- substring(rownames(top2x),1,9)

top2x.annot <- merge(x = top2x, y = annot[, c(1, 3)], by.x = 0, by.y = 1, all.x = T, all.y = F) # combine the datasets based on the ATG gene numbers; annotate with the short description
top2x.annot

write.csv(top2x.annot,"~/Desktop/colshaSUNalltimeEdgeRminusbad2x.csv")




## TO SUBSET GENES UNDERLYING QTL PEAKS#####
#--  need to have a list of ATG #'s of genes that fall in range of qtl peak.  
# might want to look at a less stringent cutoff (not 4X, etc) to look for DE in that area
# can also just pull out the likely genes (FLC, FRI, etc)
#fun stuff!!  :)   I like thinking about genes!  

##Looking at expression differences in candidate genes near QTL peaks
#this list was made empirically

#Read in Atg #'s for candidate genes
candidates <- read.csv("boltdate candidate summary.csv")
#make sure that the column with the gene names is the same type as the topltr.annot
candidates$ATG.. <- as.character(candidates$ATG..)
typeof(candidates$ATG..)
typeof(toplrt.annot$gene)  #both are character.  

#Search DE results for the atg numbers (use full toptags result, not just the 4x)
toplrt.annot[toplrt.annot$gene %in% candidates$ATG..,]  #only EMF1 comes out as different

##could the outliers be throwing this off?  Seems like FLC, etc should be differentially expressed?
##should I be using a different model (genotype * time instead of + ???)

# visualize the DE data

##this sets the cutoff for significance for de of genes...when you plotSmear, these will be a different color.  
de2 <- decideTestsDGE(lrt, p.value = 0.000001)  

summary(de2)

##this says to only include the rows where there are differentially expressed genes
detags <- rownames(dge)[as.logical(de2)] 


##and it just knows that the de.tags will be a different color
plotSmear(lrt, de.tags=detags)
abline(h=c(-2, 2), col="blue")


#### add the gene names for the top10 DE genes
detags <- rownames(dge)[as.logical(decideTestsDGE(lrt, p.value = 0.000001))]
plotSmear(lrt, de.tags=detags)
abline(h=c(-2, 2), col="blue")
text(x = top4x.annot$logCPM, # x and y give the coordinates to place the labels
     y = top4x.annot$logFC, 
     labels = top4x.annot$Row.names, # labels are the gene names!
     cex = 0.7, # size of text
     adj = c(0,0), # aligns text around the x,y point
     srt = 25) # angle of text
