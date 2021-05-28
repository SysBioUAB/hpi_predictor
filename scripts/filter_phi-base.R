#load phi_base dataset csv format
#phi.base_current <- read.csv("~/Javi/prediction_tool/version_sept_2020_marc/version_definitiva/phi-base_current.csv")
args <- commandArgs(trailingOnly = TRUE)

pkg <- c("dplyr", "data.table")
# Install packages not yet installed
installed_packages <- pkg %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(pkg[!installed_packages])
}
# Packages loading
suppressMessages(invisible(lapply(pkg, library, character.only = TRUE)))

phi.base_current <- read.csv(args[1], head = TRUE)
#disccard last (empty) columns
phi.base_current <- phi.base_current[,1:86]
#retain only Protein.ID and Mutant.Phenotype
phibase <- phi.base_current[, c("Protein.ID", "Mutant.Phenotype")]

#Find a tag for each unique Protein.ID. Some Protein.ID only appear once, so the tag is the assigned tag. Some others appear more than once with discrepant tags. Here, choose the most abundant one. If there are two or more with the most abundance, disccard all entries with this protein.ID as we cannot choose a consensus
dff <- phibase %>% dplyr::count(Protein.ID, Mutant.Phenotype) 
dff2 = dff[order(dff[,'Protein.ID'],-dff[,'n']),]
dff3 <- setDT(dff2)[, if(.N==1) .SD , .(Protein.ID, n)]
dff3 = dff3[!duplicated(dff3$Protein.ID),]
#remove entries with no Protein.ID
dff4 <- subset(dff3, dff3$Protein.ID != "no data found")
chosen_phenotype <- c("increased virulence (hypervirulence)", "loss of pathogenicity", "reduced virulence", "unaffected pathogenicity", "lethal")
dff5 <- dff4[dff4$Mutant.Phenotype %in% chosen_phenotype,]
dff5$mutant_phen_score <- NA
dff5$mutant_phen_score[dff5$Mutant.Phenotype == "increased virulence (hypervirulence)"] <- 0.5
dff5$mutant_phen_score[dff5$Mutant.Phenotype == "loss of pathogenicity"] <- 1
dff5$mutant_phen_score[dff5$Mutant.Phenotype == "unaffected pathogenicity"] <- 0
dff5$mutant_phen_score[dff5$Mutant.Phenotype == "lethal"] <- 1
dff5$mutant_phen_score[dff5$Mutant.Phenotype == "reduced virulence"] <- 0.5
dff5 <- dff5[,-2]
write.table(dff5, args[2], sep = "\t", row.names = FALSE, quote = FALSE)
