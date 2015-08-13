library(qtl)

setwd("/Users/Cody_2/Box Sync/Col_Sha_Paper/data")
colsha_traits <- read.cross("csvsr", genfile ="colsha_map.csv", 
                         phefile="col_sha_blups.csv", genotypes=c("AA","BB"), na.strings = "NA")
head(colsha_traits)
class(colsha_traits)[1] <- "riself"
colsha_traits <- jittermap(colsha_traits)
colsha_traits

# some individuals are missing genotype data.

colsha_traits <- est.rf(colsha_traits)
plot.rf(colsha_traits) 

newmap <- est.map(colsha_traits,verbose=T,error.prob=.01)
plot.map(colsha_traits,newmap) #some compression in this colsha_traits set
colsha_traits

colsha_traits <- replace.map(colsha_traits,newmap) #use new map
plot(colsha_traits) 
plot.map(colsha_traits)

colsha_traits <- sim.geno(colsha_traits,step=1,n.draws=64) 

?calc.genoprob()
#just calculates the probabilities at different locations
colsha_traits <- calc.genoprob(colsha_traits, step = 1)
colsha_traits <- calc.genoprob(colsha_traits, error.prob=0.01)

scanone.perm.imp.1 <- scanone(colsha_traits, method = "imp", pheno.col = 1:18, n.perm = 1000) 
summary(scanone.perm.imp.1) 

perm95.1 <- summary(scanone.perm.imp.1)[1]
perm95.1

scanone.imp.all <- scanone(colsha_traits, pheno.col = 1:17, method = "imp")
plot(scanone.imp.all, bandcol = "gray90", pheno.col = 1:17)
abline(h=perm95.1,lty=2) 
head(scanone.imp.all)

scanone.imp.16 <- scanone(colsha_traits, pheno.col = 17, method = "imp")
plot(scanone.imp.16, bandcol = "gray90")

scanone.imp.15 <- cim(colsha_traits, pheno.col = 15, method = "imp")
plot(scanone.imp.15, bandcol = "gray90")

library(ggplot2)
head(scanone.imp.all)
peak <- max(scanone.imp.all$Area_20131113)
area_1 <- ggplot(scanone.imp.all)
area_1
area_1 <- area_1 +  theme_bw() + geom_line(aes(x = pos, y = Area_20131113), size = 2) +
                        geom_hline(yintercept = 2.58, color = "red", size = 1) +
                        # geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        facet_grid(. ~ chr) +
                        xlab("Genetic Distance (cM)") +
                        ylab("LOD Score") +
                        ggtitle("LOD Curves for QTLs") 
area_1

?ggsave

setwd("/Users/Cody_2/Box Sync/Col_Sha_Paper/output")
ggsave("area_1.pdf", area_1)


lapply(list(df1,df2), function(x)

lapply(list("value1","value2"), function(i) ggplot(x,aes_string(x=i,y="response"))+geom_point() )









rils <- data.frame(RILs = "")
for (trait in varlist) {
  rils <- merge(rils, data.frame(RILs=rownames(coef(models[[trait]])$genotype),                          
                                  placeholder=coef(models[[trait]])$genotype[,1]),
                all.y = T)
  colnames(rils)[length(colnames(rils))] <- trait
}

head(scanone.imp.all)
traitlist <- names(scanone.imp.all)[c(3:19)]
traitlist

for (trait in traitlist) {
peak <- max(scanone.imp.all$trait)
traitplot <- ggplot(scanone.imp.all)
traitplot <- traitplot +  theme_bw() + geom_line(aes(x = pos, y = trait), size = 2) +
                        geom_hline(yintercept = 2.58, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        facet_grid(. ~ chr) +
                        xlab("Genetic Distance (cM)") +
                        ylab("LOD Score") +
                        ggtitle("LOD Curves for QTLs") 
traitplot
}

head(scanone.imp.all)

for(i in traitlist){
  peak <- max(scanone.imp.all$i)
  plt <- ggplot(scanone.imp.all, aes_string(x = "pos", y = i)) +
         geom_line(size = 2) + facet_grid(. ~ chr) +
         theme_bw() + 
         geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
         geom_hline(yintercept = 2.58, color = "red", size = 1) +
         theme(text = element_text(size = 20)) +
         xlab("Genetic Distance (cM)") +
         ylab("LOD Score") +
         ggtitle("LOD Curves for QTLs") 
  # print(plt)
  ggsave(sprintf("%s.pdf", i))
}







