library(lme4)

setwd("/Users/Cody_2/Box Sync/Col_Sha_Paper/data")
colsha  <- read.table("col_sha_biomass_area_flr_combined.csv", 
                      header=TRUE, sep = ",", na.strings = c("-","NA"),
                      stringsAsFactors = FALSE)
head(colsha)

#add shelf 
colsha$Flat_Number
colsha$shelf <- sub("(3)(\\w)(\\d)(\\d+)","\\3", colsha$Flat_Number)
colsha$shelf

colsha <- subset(colsha, treatment == "Sun" )
head(colsha)
tail(colsha)

#the last three digit column value is the true number of these rils
# need to remove the genotypes for the GWASß
colsha$original.number
colsha$genotype <- sub("(13RV)(00)(\\d)","RIL_\\3", colsha$original.number)
colsha$genotype <- sub("(13RV)(0)(\\d)","RIL_\\3", colsha$genotype)
colsha$genotype <- sub("(13RV)(\\d+)","RIL_\\2", colsha$genotype)
colsha$genotype

colsha <- colsha[grep("RIL_", colsha$genotype),]
colsha$genotype

# write table for backup
setwd("/Users/Cody_2/Box Sync/Col_Sha_Paper/data")
write.table(colsha, file="col_sha_ril_phenotypes_clean.csv", sep = ",")

head(colsha)
str(colsha)
colsha$shelf <- as.factor(colsha$shelf)

drybiomass_model1 <- lmer(biomass_dry_leaf ~ (1|shelf) +(1|genotype), 
                          data = colsha, REML = FALSE)
summary(drybiomass_model1)

drybiomass_model2 <- lmer(biomass_dry_leaf ~ (1|shelf),
                          data = colsha, REML = FALSE)
anova(drybiomass_model1,drybiomass_model2)

drybiomass_model3 <- lmer(biomass_dry_leaf ~ (1|genotype),
                          data = colsha, REML = FALSE)
anova(drybiomass_model1,drybiomass_model3)

#floweringtime
boltdays_model1 <- lmer(boltdays ~ (1|shelf) + (1|genotype), 
                          data = colsha, REML = FALSE)
summary(boltdays_model1)
boltdays_model2 <- lmer(boltdays ~ (1|genotype),
                          data = colsha, REML = FALSE)
anova(boltdays_model1,boltdays_model2)


boltdays_model3 <- lmer(boltdays ~ (1|shelf),
                          data = colsha, REML = FALSE)
anova(boltdays_model1,boltdays_model3)


boltdays_model4 <- lmer(boltdays ~ (1|genotype),
                          data = colsha, REML = FALSE)
anova(boltdays_model1,boltdays_model4)

#seems like there is a shelf effect for bolting date
names(colsha)
str(colsha)
#growth 
colsha$area_growth <- colsha$Area_20131113 - colsha$Area_20131108
colsha$perimeter_growth <- colsha$Perimeter_20131113 - colsha$Perimeter_20131108



varlist <- names(colsha)[c(7:21,24:25)]
varlist

#apply same model for all traitsß
models <- lapply(varlist, function(x) {
    lmer(substitute(i ~ (1|shelf) + (1|genotype), list(i = as.name(x))), 
                      data = colsha)
})

#take a look
str(models)

models[[1]]
models[[2]]
#print a few traits to take a look
varlist
coef(models[[1]])$genotype
coef(models[[2]])$genotype
coef(models[[5]])$genotype
#name the model list
names(models) <- varlist

#extract the blups!!!!
#this works, but could be cleaned up
#it extracts random intercept (mean) of each RIL
rils <- data.frame(RILs = "")
for (trait in varlist) {
  rils <- merge(rils, data.frame(RILs=rownames(coef(models[[trait]])$genotype),                          
                                  placeholder=coef(models[[trait]])$genotype[,1]),
                all.y = T)
  colnames(rils)[length(colnames(rils))] <- trait
}

head(rils)
head(rils)
tail(rils)

rils$RILs <- as.character(rils$RILs)
str(rils)

head(rils)
tail(rils)


#format for RQTL
rils.t <- as.data.frame(t(rils))
head(rils.t)
colnames(rils.t)
dim(rils.t)
rils.t[19,] <- rils.t[1,]
rownames(rils.t)[19] <- "id"
rownames(rils.t)
rils.t <- rils.t[-1,]
head(rils.t)
tail(rils.t)

# write table for backup
setwd("/Users/Cody_2/Box Sync/Col_Sha_Paper/data")
write.table(rils.t, file="col_sha_blups.csv", sep = ",")



