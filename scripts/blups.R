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

setwd("/Users/Cody_2/Box Sync/Col_Sha_Paper/data")
write.table(colsha, file="col_sha_ril_phenotypes_clean.csv", sep = ",")

head(colsha)
str(colsha)

# okay started