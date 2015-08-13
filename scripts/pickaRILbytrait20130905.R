testdata<-read.csv("CompiledtraitsforRILselectionshort.csv")
summary(testdata)
head(testdata)
names(testdata)
install.packages("reshape")
library(reshape)
yodata<-melt(testdata, id.vars="Ecotype")
head(yodata)


##Plot the data to try and see which is the most different
library(ggplot2)
ggplot(yodata, aes(x = variable, y = value, color = Ecotype)) +
  geom_point(size = 3) +
  theme(axis.text.x = element_blank())+
  theme(strip.text = element_text(angle = 90, size = 12))+
  theme(strip.background = element_rect(fill = "white"))+
  facet_grid(~variable, scale = "free")
ggsave("allTheTraits_090513.pdf", height = 8, width = 10.5)
                  
ggplot(yodata, aes(x = variable, y = value, color = Ecotype)) +
  geom_point(size = 4, aes(shape=Ecotype, solid=T)) +
  scale_shape_manual(values=c(14:20)) +
  theme(axis.text.x = element_blank())+
  theme(strip.text = element_text(angle = 90, size = 12))+
  theme(strip.background = element_rect(fill = "white"))+
  facet_wrap(~variable, scales = "free", ncol=34)
ggsave("~/Desktop/allTheTraitsV3_090513.pdf", height = 8, width = 10.5)

head(testdata)
colnames(testdata)

##Plot the PCA and distance for JUST the residuals/response variables
testdataresponse<-testdata[,c(1,6,9,12,21,30,35)]
rownames(testdataresponse) <- testdataresponse$Ecotype
eco.dist.response <- dist(scale(testdataresponse[-(6:7),-1]))
loc <- as.data.frame(cmdscale(eco.dist.response))
loc$text <- row.names(loc)
eco.dist.response


pl <- ggplot(data=loc,mapping=aes(x=V1,y=V2,label=text))
pl + geom_text()

###------------------------END Response------------

##Plot the PCA and distance for LOW light phenotypes
testdatalow<-testdata[,c(1,5,8,11,15,16,20,26,29,34)]

rownames(testdatalow) <- testdatalow$Ecotype
eco.dist.low <- dist(scale(testdatalow[-(6:7),-1]))
loc <- as.data.frame(cmdscale(eco.dist.low))
loc$text <- row.names(loc)
eco.dist.low

pl <- ggplot(data=loc,mapping=aes(x=V1,y=V2,label=text))
pl + geom_text()

###------------------------END LOW------------

##Plot the PCA and distance for HIGH light phenotypes
testdatahigh<-testdata[,c(1,2,3,4,7,10,13,14,17,18,19,22,23,24,25,28,31,32,33)]

rownames(testdatahigh) <- testdatahigh$Ecotype
eco.dist.high <- dist(scale(testdatahigh[-(6:7),-1]))
loc <- as.data.frame(cmdscale(eco.dist.high))
loc$text <- row.names(loc)
eco.dist.high


pl <- ggplot(data=loc,mapping=aes(x=V1,y=V2,label=text))
pl + geom_text()

###------------------------END HIGH------------

#####Plot for all the variables
rownames(testdata) <- testdata$Ecotype
eco.dist <- dist(scale(testdata[-(6:7),-1]))
loc <- as.data.frame(cmdscale(eco.dist))
loc$text <- row.names(loc)
eco.dist


pl <- ggplot(data=loc,mapping=aes(x=V1,y=V2,label=text))
pl + geom_text()