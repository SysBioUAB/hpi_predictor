# HUMAN DUALSEQ VALUES

args <- commandArgs(trailingOnly = TRUE)

for(pkg in c("plyr", "tidyverse", "bio3d", "readxl", "stringi", "e1071", "progress", "data.table", "parallel", "Biostrings")){
    suppressMessages(library(pkg, character.only = TRUE))
}


# check args

#args[1]
#args[2]
#args[3]
#args[4]
#args[5]
#is.null(args[5])
#nchar(args[5])
#args[6]
#args[7]
#args[8]
#args[9]
#args[10]



host_dualseqdb_queries <- read.delim(args[1], header=FALSE)
# Write proper colnames
## colnames(querys2) <- c("")

#filter by %identity > %40 and fitness z score < 0.05:

host_dualseqdb_queries_q1 <- subset(host_dualseqdb_queries, host_dualseqdb_queries$V3 > 40 & host_dualseqdb_queries$V7 < 0.05)

#calculate fitness z score mean and standard deviation:


averages <- setDT(host_dualseqdb_queries_q1)[, list(fitness_mean=mean(V6)), V1]
helper <- inner_join(host_dualseqdb_queries_q1, averages, by="V1")
stds <- setDT(host_dualseqdb_queries_q1)[, list(sd=sd(V6)), V1]
helper <- inner_join(helper, stds, by="V1")
helper <- unique(helper[,c(1,2,3,4,5,8,9)])
host_dualseqdb_filtered_values <- unique(helper[,c(1,6,7)])
colnames(host_dualseqdb_filtered_values) <- c("H_ID_query", "H_query_log2FC", "H_query_log2FC_sd")
#write.table(human_dualseqdb_filtered_values, "human_dualseqdb_filtered_values.tsv", sep = "\t", row.names = FALSE, quote =  FALSE)

# BACTERIAL DUALSEQ VALUES

bac_dualseqdb_queries <- read.delim(args[2], header=FALSE)
# Write proper colnames
## colnames(querys2) <- c("")

#filter by %identity > %40 and fitness z score < 0.05:

bac_dualseqdb_queries_q1 <- subset(bac_dualseqdb_queries, bac_dualseqdb_queries$V3 > 40 & bac_dualseqdb_queries$V7 < 0.05)

#calculate fitness z score mean and standard deviation:

averages <- setDT(bac_dualseqdb_queries_q1)[, list(fitness_mean=mean(V6)), V1]
helper <- inner_join(bac_dualseqdb_queries_q1, averages, by="V1")
stds <- setDT(bac_dualseqdb_queries_q1)[, list(sd=sd(V6)), V1]
helper <- inner_join(helper, stds, by="V1")
helper <- unique(helper[,c(1,2,3,4,5,8,9)])
bac_dualseqdb_filtered_values <- unique(helper[,c(1,6,7)])
colnames(bac_dualseqdb_filtered_values) <- c("P_ID_query", "P_query_log2FC", "P_query_log2FC_sd")
#write.table(bac_dualseqdb_filtered_values, "bac_dualseqdb_filtered_values.tsv", sep = "\t", row.names = FALSE, quote =  FALSE)


# BACTERIAL BACFITBASE VALUES


bacfitbase_queries <- read.delim(args[3], header=FALSE)
# Write proper colnames
## colnames(querys2) <- c("")

#filter by %identity > %40 and fitness z score < 0.05:

bacfitbase_queries_q1 <- subset(bacfitbase_queries, bacfitbase_queries$V3 > 40 & bacfitbase_queries$V7 < 0.05)

#calculate fitness z score mean and standard deviation:

averages <- setDT(bacfitbase_queries_q1)[, list(fitness_mean=mean(V6)), V1]
helper <- inner_join(bacfitbase_queries_q1, averages, by="V1")
stds <- setDT(bacfitbase_queries_q1)[, list(sd=sd(V6)), V1]
helper <- inner_join(helper, stds, by="V1")
helper <- unique(helper[,c(1,2,3,4,5,8,9)])
bacfitbase_filtered_values <- unique(helper[,c(1,6,7)])
colnames(bacfitbase_filtered_values) <- c("P_ID_query", "P_query_fitness", "P_query_fitness_sd")
#write.table(bacfitbase_filtered_values, "bacfitbase_filtered_values.tsv", sep = "\t", row.names = FALSE, quote =  FALSE)

# BACTERIAL PHI-BASE VALUES

phibase_queries <- read.delim(args[4], header=FALSE)
dff <- phibase_queries[order(phibase_queries[,'V1'], -phibase_queries[,'V3']),]
pathogen_virulence <- dff[!duplicated(dff$V1),][,c(1,6,7)]
colnames(pathogen_virulence) <- c("P_ID_query", "virulence_tag", "virulence_tag_score")



# load filtered results to add to predicted interactome
if (nchar(args[5]) > 0) {
hi_centrality <- read.table(args[5], header = TRUE)
#colnames(hi_centrality) <- c("H_ID_query", "node_degree")
}
#load predicted interactome
interactome <- read.delim(args[6], header = TRUE)

#merge

interactome_f <- merge(interactome, host_dualseqdb_filtered_values, by="H_ID_query", all = TRUE)

#colnames(PAO1_dualseqdb_filtered_values_GEIM800108) <- c("P_ID_query", "P_query_log2FC", "P_query_log2FC_sd")

interactome_f <- merge(interactome_f, bac_dualseqdb_filtered_values, by = "P_ID_query", all = TRUE)

#colnames(PAO1_bacfitbase_filtered_values_GEIM800108) <- c("P_ID_query", "P_query_fitness", "P_query_fitness_sd")

interactome_f <- merge(interactome_f, bacfitbase_filtered_values, by = "P_ID_query", all = TRUE)

#colnames(hi_centrality) <- c("H_ID_query", "H_query_node_degree")
if (nchar(args[5]) > 0){
interactome_f <- merge(interactome_f, hi_centrality, by = "H_ID_query", all = TRUE)
}

interactome_f <- merge(interactome_f, pathogen_virulence, by = "P_ID_query", all  = TRUE)

interactome_def <- interactome_f[!is.na(interactome_f$H_ID_template),]
interactome_def$norm_H_query_log2FC <- (interactome_def$H_query_log2FC-min(interactome_def$H_query_log2FC, na.rm = TRUE))/(max(interactome_def$H_query_log2FC, na.rm = TRUE)-min(interactome_def$H_query_log2FC, na.rm = TRUE))
interactome_def$norm_P_query_log2FC <- (interactome_def$P_query_log2FC-min(interactome_def$P_query_log2FC, na.rm = TRUE))/(max(interactome_def$P_query_log2FC, na.rm = TRUE)-min(interactome_def$P_query_log2FC, na.rm = TRUE))
interactome_def$norm_P_query_fitness <- (interactome_def$P_query_fitness-max(interactome_def$P_query_fitness, na.rm = TRUE))/(min(interactome_def$P_query_fitness, na.rm = TRUE)-max(interactome_def$P_query_fitness, na.rm = TRUE))


if (nchar(args[5]) > 0){
interactome_def$weighted_score <- rowMeans(interactome_def[,c(35, 37:40)], na.rm=TRUE)
#number of Nas in score
interactome_def$confidence <- apply(interactome_def[,c(35, 37:40)], MARGIN = 1, function(x) sum(is.na(x)))
interactome_def$confidence <- abs(interactome_def$confidence - 5)
#interactome_def$weighted_score <- interactome_def$score/(interactome_def$score_na+1)
#write.table(interactome_def, "/home/rocio/Javi/prediction_tool/version_sept_2020_marc/version_definitiva/toy_example_directory_2/results/Acinetobacter_baumannii_17978_vs_Homo_sapiens_predicted_interactome_def.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

}

if ( nchar(args[5]) == 0){
interactome_def$weighted_score <- rowMeans(interactome_def[,c(33:36)], na.rm=TRUE)
interactome_def$confidence <- apply(interactome_def[,c(33:36)], MARGIN = 1, function(x) sum(is.na(x)))
interactome_def$confidence <- abs(interactome_def$confidence - 4)

}




write.table(interactome_def, args[7], sep = "\t", row.names = FALSE, quote = FALSE)


if (nchar(args[5]) > 0){
interactome_def_only_scores <- interactome_def[, c(1,2,41,42)]
write.table(interactome_def_only_scores, args[8], sep = "\t", row.names = FALSE, quote = FALSE)
}

if (nchar(args[5]) == 0){

    interactome_def_only_scores <- interactome_def[, c(1,2,37,38)]
    write.table(interactome_def_only_scores, args[8], sep = "\t", row.names = FALSE, quote = FALSE)
}


P_log2fc <- unique(interactome_def[, c("P_ID_query", "P_query_log2FC")])
H_log2fc <- unique(interactome_def[, c("H_ID_query", "H_query_log2FC")])
colnames(P_log2fc) <- c("ID_query", "log2fc")
colnames(H_log2fc) <- c("ID_query", "log2fc")
log2fc <- rbind.data.frame(P_log2fc, H_log2fc)
#write.table(log2fc, "/home/rocio/Javi/prediction_tool/version_sept_2020_marc/version_definitiva/toy_example_directory_2/results/Acinetobacter_baumannii_17978_vs_Homo_sapiens_predicted_interactome_def_log2fc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(log2fc, args[9], sep = "\t", row.names = FALSE, quote = FALSE)

P_ID <- cbind.data.frame(unique(interactome_def$P_ID_query), rep("Pathogen", length(unique(interactome_def$P_ID_query))))
H_ID <- cbind.data.frame(unique(interactome_def$H_ID_query), rep("Host", length(unique(interactome_def$H_ID_query))))
colnames(P_ID) <- c("ID", "Organism")
colnames(H_ID) <- c("ID", "Organism")
IDs <- rbind.data.frame(P_ID,H_ID)
#write.table(IDs, "/home/rocio/Javi/prediction_tool/version_sept_2020_marc/version_definitiva/toy_example_directory_2/results/Acinetobacter_baumannii_17978_vs_Homo_sapiens_predicted_interactome_def_IDs.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(IDs, args[10], sep = "\t", row.names = FALSE, quote = FALSE)


