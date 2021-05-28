# HUMAN DUALSEQ VALUES

args <- commandArgs(trailingOnly = TRUE)

for(pkg in c("plyr", "tidyverse", "bio3d", "readxl", "stringi", "e1071", "progress", "data.table", "parallel", "Biostrings", "crayon")){
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
bac_dualseqdb_queries <- read.delim(args[2], header=FALSE)
bacfitbase_queries <- read.delim(args[3], header=FALSE)
phibase_queries <- read.delim(args[4], header=FALSE)
hi_centrality <- read.table(args[5], header = TRUE)[,c(1,3,5)]
interactome <- read.delim(args[6], header = TRUE)


#host_dualseqdb_queries <- read.delim("/home/rocio/copia_predictor_mac/hpi_predictor/results/Homo_sapiens_9606_vs_coxiella_burnetii_pos_custom/false_positive_rate_0_1_predictions/intermediate_files/host_vs_dualseq_values.tsv", header=FALSE)
#bac_dualseqdb_queries <- read.delim("/home/rocio/copia_predictor_mac/hpi_predictor/results/Homo_sapiens_9606_vs_coxiella_burnetii_pos_custom/false_positive_rate_0_1_predictions/intermediate_files/bacteria_vs_dualseq_values.tsv", header=FALSE)
#bacfitbase_queries <- read.delim("/home/rocio/copia_predictor_mac/hpi_predictor/results/Homo_sapiens_9606_vs_coxiella_burnetii_pos_custom/false_positive_rate_0_1_predictions/intermediate_files/bacteria_vs_bacfitbase_values.tsv", header=FALSE)
#phibase_queries <- read.delim("/home/rocio/copia_predictor_mac/hpi_predictor/results/Homo_sapiens_9606_vs_coxiella_burnetii_pos_custom/false_positive_rate_0_1_predictions/intermediate_files/bacteria_vs_bacfitbase_values.tsv", header=FALSE)
#hi_centrality <- read.table("/home/rocio/copia_predictor_mac/hpi_predictor/databases/hi_centrality_9606.tsv", header = TRUE)[,c(1,3,5)]
#interactome <- read.delim("/home/rocio/copia_predictor_mac/hpi_predictor/results/Homo_sapiens_9606_vs_coxiella_burnetii_pos_custom/false_positive_rate_0_1_predictions/interactomes/consensus_interactome_50_percent_of_predictors.tsv", header = TRUE)

#filter by %identity > %40 and fitness z score < 0.05:


if(nrow(host_dualseqdb_queries[complete.cases(host_dualseqdb_queries),])){
host_dualseqdb_queries_q1 <- subset(host_dualseqdb_queries, host_dualseqdb_queries$V3 > 40 & host_dualseqdb_queries$V7 < 0.05)

#calculate fitness z score mean and standard deviation:


averages <- setDT(host_dualseqdb_queries_q1)[, list(fitness_mean=mean(V6)), V1]
helper <- inner_join(host_dualseqdb_queries_q1, averages, by="V1")
stds <- setDT(host_dualseqdb_queries_q1)[, list(sd=sd(V6)), V1]
helper <- inner_join(helper, stds, by="V1")
helper <- unique(helper[,c(1,2,3,4,5,8,9)])
host_dualseqdb_filtered_values <- unique(helper[,c(1,6,7)])
colnames(host_dualseqdb_filtered_values) <- c("H_ID_query", "H_query_log2FC", "H_query_log2FC_sd")
}


# BACTERIAL DUALSEQ VALUES

#filter by %identity > %40 and fitness z score < 0.05:
if(nrow(bac_dualseqdb_queries[complete.cases(bac_dualseqdb_queries),])){

bac_dualseqdb_queries_q1 <- subset(bac_dualseqdb_queries, bac_dualseqdb_queries$V3 > 40 & bac_dualseqdb_queries$V7 < 0.05)

#calculate fitness z score mean and standard deviation:

averages <- setDT(bac_dualseqdb_queries_q1)[, list(fitness_mean=mean(V6)), V1]
helper <- inner_join(bac_dualseqdb_queries_q1, averages, by="V1")
stds <- setDT(bac_dualseqdb_queries_q1)[, list(sd=sd(V6)), V1]
helper <- inner_join(helper, stds, by="V1")
helper <- unique(helper[,c(1,2,3,4,5,8,9)])
bac_dualseqdb_filtered_values <- unique(helper[,c(1,6,7)])
colnames(bac_dualseqdb_filtered_values) <- c("P_ID_query", "P_query_log2FC", "P_query_log2FC_sd")

}

# BACTERIAL BACFITBASE VALUES

#filter by %identity > %40 and fitness z score < 0.05:
if(nrow(bacfitbase_queries[complete.cases(bacfitbase_queries),])){
bacfitbase_queries_q1 <- subset(bacfitbase_queries, bacfitbase_queries$V3 > 40 & bacfitbase_queries$V7 < 0.05)

#calculate fitness z score mean and standard deviation:

averages <- setDT(bacfitbase_queries_q1)[, list(fitness_mean=mean(V6)), V1]
helper <- inner_join(bacfitbase_queries_q1, averages, by="V1")
stds <- setDT(bacfitbase_queries_q1)[, list(sd=sd(V6)), V1]
helper <- inner_join(helper, stds, by="V1")
helper <- unique(helper[,c(1,2,3,4,5,8,9)])
bacfitbase_filtered_values <- unique(helper[,c(1,6,7)])
colnames(bacfitbase_filtered_values) <- c("P_ID_query", "P_query_fitness", "P_query_fitness_sd")
}


# BACTERIAL PHI-BASE VALUES

if(nrow(phibase_queries[complete.cases(phibase_queries),])){

dff <- phibase_queries[order(phibase_queries[,'V1'], -phibase_queries[,'V3']),]
pathogen_virulence <- dff[!duplicated(dff$V1),][,c(1,6,7)]
colnames(pathogen_virulence) <- c("P_ID_query", "virulence_tag", "virulence_tag_score")
}

#merge

if (exists("host_dualseqdb_filtered_values")){
interactome_f <- merge(interactome, host_dualseqdb_filtered_values, by="H_ID_query", all = TRUE)
}

if (exists("bac_dualseqdb_filtered_values")){
interactome_f <- merge(interactome_f, bac_dualseqdb_filtered_values, by = "P_ID_query", all = TRUE)
}

if (exists("bacfitbase_filtered_values")){
interactome_f <- merge(interactome_f, bacfitbase_filtered_values, by = "P_ID_query", all = TRUE)
}

if (nchar(args[5]) > 0){
interactome_f <- merge(interactome_f, hi_centrality, by = "H_ID_query", all = TRUE)
}
if (exists ("pathogen_virulence")){
interactome_f <- merge(interactome_f, pathogen_virulence, by = "P_ID_query", all  = TRUE)
}

interactome_def <- interactome_f[!is.na(interactome_f$H_ID_template),]
if (exists("host_dualseqdb_filtered_values")){
interactome_def$norm_H_query_log2FC <- (interactome_def$H_query_log2FC-min(interactome_def$H_query_log2FC, na.rm = TRUE))/(max(interactome_def$H_query_log2FC, na.rm = TRUE)-min(interactome_def$H_query_log2FC, na.rm = TRUE))
}
if (exists("bac_dualseqdb_filtered_values")){
interactome_def$norm_P_query_log2FC <- (interactome_def$P_query_log2FC-min(interactome_def$P_query_log2FC, na.rm = TRUE))/(max(interactome_def$P_query_log2FC, na.rm = TRUE)-min(interactome_def$P_query_log2FC, na.rm = TRUE))
}
if (exists("bacfitbase_filtered_values")){
interactome_def$norm_P_query_fitness <- (interactome_def$P_query_fitness-max(interactome_def$P_query_fitness, na.rm = TRUE))/(min(interactome_def$P_query_fitness, na.rm = TRUE)-max(interactome_def$P_query_fitness, na.rm = TRUE))
}

norms <- grep("norm", names(interactome_def), value=TRUE)
interactome_def$weighted_score <- rowMeans(interactome_def[,norms], na.rm = TRUE)
interactome_def$confidence <- apply(interactome_def[,norms], MARGIN = 1, function(x) sum(is.na(x)))
interactome_def$confidence <- abs(interactome_def$confidence - length(norms))

write.table(interactome_def, args[7], sep = "\t", row.names = FALSE, quote = FALSE)
interactome_def_only_scores <- interactome_def[,c("P_ID_query", "H_ID_query","weighted_score", "confidence")]
write.table(interactome_def_only_scores, args[8], sep = "\t", row.names = FALSE, quote = FALSE)

