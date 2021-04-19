args <- commandArgs(trailingOnly = TRUE)


# pasar como argumento el working directory desde bash a Rscript
#setwd("/home/rocio/Javi/prediction_tool/version_sept_2020_marc/version_definitiva/toy_example_directory_2/")
setwd(args[9])


#for(pkg in c("tidyverse", "bio3d", "readxl", "stringi", "e1071", "progress", "data.table", "parallel", "Biostrings")){
#  library(pkg, character.only = TRUE)
#}

pkg <- c("plyr", "tidyverse", "bio3d", "readxl", "stringi", "e1071", "progress", "data.table", "parallel", "Biostrings")
# Install packages not yet installed
installed_packages <- pkg %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(pkg[!installed_packages])
}
# Packages loading
suppressMessages(invisible(lapply(pkg, library, character.only = TRUE)))

#args <- commandArgs(trailingOnly = TRUE)

#host_org <- args[1]
host_org <- gsub(args[1], pattern = " ", replacement = "_")
#host_org <- gsub(host_org, pattern = "\\(/| / | \\(| ", replacement =  "_")
#host_org <- gsub(host_org, pattern = "\\)", replacement = "")
#host_org <- gsub(host_org, pattern = "\\.$", replacement = "")

pathogen_org <- gsub(args[2], pattern = " ", replacement = "_")

#pathogen_org <- gsub(pathogen_org, pattern = "\\(/| / | \\(| ", replacement =  "_")

#pathogen_org <- gsub(pathogen_org, pattern = "\\)", replacement = "")

pathogen_org <- gsub(pathogen_org, pattern = "\\.", replacement = "")


#ayuda="ayuda"
#system(sprintf("value=$(echo %s)", ayuda))


#cat(sprintf("pathogen organism %s\n", pathogen_org))

host_dir <- args[3]
#print("args3 ok")
pathogen_dir <- args[4]
#print("args4 ok")
pos_inter_dir <- args[5]
#print("args5 ok")
lagmax <- as.numeric(args[6])
#print("args6 ok")
cutoff <- as.numeric(args[7])
#print("args7 ok")
cutoff_h <- gsub(cutoff, pattern = "\\.", replacement = "_")
results_dir <- args[8]
#print("args8 ok")
new_dir <- paste(results_dir, host_org, "_vs_", pathogen_org, "/", sep = "")

#Distance function
euc.dist <- function(v1,v2){
  max(ccf(unlist(v1),unlist(v2), plot = F, lag.max = lagmax)$acf)
}

cores=detectCores()
print(paste(cores, " CPU cores detected, using ", cores-1, sep = ""))

######################################

#############################
# Build reference profiles #
#############################

# Hay que rehacer el directorio donde esté alojado este archivo y los siguientes, para el ordenador del usuario

#learn_set <- read_excel("/home/rocio/Javi/prediction_tool/version_sept_2020_marc/version_definitiva/toy_example_directory_2/datasets/positive_interactions/phi_data.xlsx")
learn_set <- read.delim(args[10])
learn_set <- filter(learn_set, P_len < 2000 & P_len > 100 & H_len < 2000 & H_len > 100)

p.learn <- learn_set %>% select(P_ID, P_seq) %>% distinct(P_seq, .keep_all = TRUE)
h.learn <- learn_set %>% select(H_ID, H_seq) %>% distinct(H_seq, .keep_all = TRUE)

#Default AAindex indices
#97 Alpha-helix_indices_(Geisow-Roberts,_1980) -- GEIM800101
#101 Beta-strand_indices_(Geisow-Roberts,_1980) -- GEIM800105
#104 Aperiodic_indices_(Geisow-Roberts,_1980) -- GEIM800108
#398 Hydrophobicity_(Zimmerman_et_al.,_1968) -- ZIMJ680101
#401 Isoelectric_point_(Zimmerman_et_al.,_1968) -- ZIMJ680104
#476 Electron-ion_interaction_potential_values_(Cosic,_1994) -- COSI940101

if (file.exists(args[12]) == FALSE) {
    toy_indices <- read.table(args[11], head = FALSE)
} else {
    toy_indices <- read.table(args[12], head = FALSE)
}

descriptors_vec <- vector()

for (i in 1:nrow(toy_indices)) {
  
  descriptors_vec[i] <- as.character(toy_indices[i,1])
}

#print(host_org)
#print(pathogen_org)
#print(host_dir)
#print(pathogen_dir)
#print(pos_inter_dir)
#organisms <- read.delim("/home/rocio/Javi/prediction_tool/version_sept_2020_marc/version_definitiva/toy_example_directory_2/results/organisms", header  = FALSE)

#organisms_vec <- vector()

#for (i in 1:nrow(organisms)) {

#  organisms_vec[i] <- as.character(organisms[i,1])
#}





#indices <- c("GEIM800101", "GEIM800105", "GEIM800108", "ZIMJ680101", "ZIMJ680104")
#i=1
#indices <- "ZIMJ680104"


p.learn.index <- list()
h.learn.index <- list()
for(i in 1:length(descriptors_vec)){
 
    start_time <- Sys.time()
    helper_r_i_h <- paste(pos_inter_dir, "host/", descriptors_vec[i], "_host.tsv", sep = "")
    helper_r_i_p <- paste(pos_inter_dir, "pathogen/", descriptors_vec[i], "_pathogen.tsv", sep = "")
    if (file.exists(helper_r_i_p) == FALSE){
        cat(sprintf("Building the pathogen reference database"))

        max.val <- max(aa.index[[descriptors_vec[i]]]$I)
        min.val <- min(aa.index[[descriptors_vec[i]]]$I)
        calculation <- lapply(1:length(unlist(p.learn$P_seq)), function(x){(aa2index(aa = as.vector(stri_list2matrix(str_extract_all(p.learn$P_seq[x], boundary("character")), byrow = TRUE)), index = descriptors_vec[i], window = 9)-min.val)/(max.val-min.val)})
        calculation <- lapply(1:length(calculation), function(x){unlist(calculation[x])-mean(unlist(calculation[x]))})
        calculation <- lapply(1:length(calculation), function(x){unlist(calculation[x])*hamming.window(length(unlist(calculation[x])))})
        calculation <- lapply(1:length(calculation), function(x){c(unlist(calculation[x]), rep(0, 4096-length(unlist(calculation[x]))))})
        p.learn.index[[i]] <- calculation
        print(helper_r_i_p)
        write.table(calculation, helper_r_i_p, sep = "\t", row.names = FALSE, quote = FALSE)       
        print("Done!")
        end_time <- Sys.time()
        print(end_time - start_time)
    } else {
#        p.learn.index[[i]] <- read.delim(helper_r_i_p)
    #print("Retrieving from database")
        cat(sprintf("Retrieving %s from pathogen database\n", descriptors_vec[[i]]))
        p.learn.index[[i]] <- read.delim(helper_r_i_p)
    }

    start_time <- Sys.time()
    if (file.exists(helper_r_i_h) == FALSE){
        cat(sprintf("Building the host reference database"))

        max.val <- max(aa.index[[descriptors_vec[i]]]$I)
        min.val <- min(aa.index[[descriptors_vec[i]]]$I)
        calculation <- lapply(1:length(unlist(h.learn$H_seq)), function(x){(aa2index(aa = as.vector(stri_list2matrix(str_extract_all(h.learn$H_seq[x], boundary("character")), byrow = TRUE)), index = descriptors_vec[i], window = 9)-min.val)/(max.val-min.val)})
        calculation <- lapply(1:length(calculation), function(x){unlist(calculation[x])-mean(unlist(calculation[x]))})
        calculation <- lapply(1:length(calculation), function(x){unlist(calculation[x])*hamming.window(length(unlist(calculation[x])))})
        calculation <- lapply(1:length(calculation), function(x){c(unlist(calculation[x]), rep(0, 4096-length(unlist(calculation[x]))))})
        h.learn.index[[i]] <- calculation
        write.table(calculation, helper_r_i_h, sep = "\t", row.names = FALSE, quote = FALSE)
        print("Done!")
        rm(calculation); gc()
        end_time <- Sys.time()
        print(end_time - start_time)
    } else {
#        h.learn.index[[i]] <- read.delim(helper_r_i_h)
#        print("Retrieving from database")
        cat(sprintf("Retrieving %s from host database\n", descriptors_vec[[i]]))
        h.learn.index[[i]] <- read.delim(helper_r_i_h)
    }

}

#################################################
# Load and calculate data for host and pathogen #
#################################################

help_p <- paste(pathogen_dir, pathogen_org, "/", pathogen_org, "_proteome.tsv", sep = "")
pathogen <- read.delim(help_p, header = TRUE)
colnames(pathogen) <- c("P_ID", "Pathogen", "Pathogen_taxon_ID", "P_len", "P_seq")
pathogen <-  filter(pathogen, P_len < 2000 & P_len > 100)
p.proteome <- pathogen %>% select(P_ID, P_seq) %>% distinct(P_seq, .keep_all = TRUE)


help_h <- paste(host_dir, host_org, "/", host_org, "_proteome.tsv", sep = "")
host <- read.delim(help_h, header = TRUE)
colnames(host) <- c("H_ID", "Host", "Host_taxon_ID", "H_len", "H_seq")
#colnames(host) <- c("H_ID", "H_AC", "H_len", "H_seq")
host <- filter(host, H_len < 2000 & H_len > 100)
h.proteome <- host %>% select(H_ID, H_seq) %>% distinct(H_seq, .keep_all = TRUE)


p.proteome.index <- list()
h.proteome.index <- list()

for(i in 1:length(descriptors_vec)){
start_time <- Sys.time()
#print("Building the pathogen proteome database")


  helper_i_p <- paste(pathogen_dir, pathogen_org, "/descriptors/p_results_lagmax_", lagmax, "_cutoff_", cutoff_h, "_",  descriptors_vec[i], ".txt", sep = "")
    if (file.exists(helper_i_p) == FALSE){
        cat(sprintf("Building the pathogen proteome database"))

#help_d_p <- paste("home/rocio/Javi/prediction_tool/version_sept_2020_marc/version_definitiva/toy_example_directory_2/datasets/pathogens/", organisms_vec[2], "/descriptors/", "p_", sep = "" )

        calculation <- lapply(1:length(unlist(p.proteome$P_seq)), function(x){aa2index(aa = as.vector(stri_list2matrix(str_extract_all(p.proteome$P_seq[x], boundary("character")), byrow = TRUE)), index = descriptors_vec[i], window = 9)})
        calculation <- lapply(1:length(calculation), function(x){unlist(calculation[x])-mean(unlist(calculation[x]))})
        calculation <- lapply(1:length(calculation), function(x){unlist(calculation[x])*hamming.window(length(unlist(calculation[x])))})
        calculation <- lapply(1:length(calculation), function(x){c(unlist(calculation[x]), rep(0, 4096-length(unlist(calculation[x]))))})
        p.proteome.index[[i]] <- calculation

        cat(sprintf("Done!"))
        rm(calculation); gc()
        end_time <- Sys.time()
        print(end_time - start_time)
    } else {
        cat(sprint("Descriptor already in pathogen database"))
    }

    start_time <- Sys.time()
#print("Building the human proteome database")
    helper_i_h <- paste(host_dir, host_org, "/descriptors/h_results_lagmax_", lagmax, "_cutoff_", cutoff_h, "_", descriptors_vec[i], ".txt", sep = "")
    if (file.exists(helper_i_h) == FALSE){
        cat(sprintf("Building the host proteome database"))
        print(helper_i_h)
        calculation <- lapply(1:length(unlist(h.proteome$H_seq)), function(x){aa2index(aa = as.vector(stri_list2matrix(str_extract_all(h.proteome$H_seq[x], boundary("character")), byrow = TRUE)), index = descriptors_vec[i], window = 9)})
        calculation <- lapply(1:length(calculation), function(x){unlist(calculation[x])-mean(unlist(calculation[x]))})
        calculation <- lapply(1:length(calculation), function(x){unlist(calculation[x])*hamming.window(length(unlist(calculation[x])))})
        calculation <- lapply(1:length(calculation), function(x){c(unlist(calculation[x]), rep(0, 4096-length(unlist(calculation[x]))))})
        h.proteome.index[[i]] <- calculation

        rm(calculation); gc()
        cat(sprintf("Done!"))
        end_time <- Sys.time()
        print(end_time - start_time)
    } else {
        cat(sprintf("Descriptor already in host database"))
    }

}

##############################
### Calculate interactomes ###
##############################


##### Profile similarity in pathogen proteins
start_time <- Sys.time()
#print("Calculating the pathogen similarity profiles")


for (i in 1:length(descriptors_vec)){
    
    helper_i_p <- paste(pathogen_dir, pathogen_org, "/descriptors/p_results_lagmax_", lagmax, "_cutoff_", cutoff_h, "_",  descriptors_vec[i], ".txt", sep = "")
    if (file.exists(helper_i_p) == FALSE){
#        print("Calculating the pathogen similarity profiles")
        cat(sprintf("Calculating the pathogen similarity profiles for descriptor %s \n", descriptors_vec[[i]]))
        p.matrix <- data.table()
#pb <- progress_bar$new(total = length(p.proteome.index[[i]]))

        total <- length(p.proteome.index[[i]])
        pb <- txtProgressBar(min = 0, max = total, style = 3)


        for(x in 1:length(p.proteome.index[[i]])) {
  #pb$tick()
            calculation <- mclapply(1:length(p.learn.index[[i]]), mc.cores=cores-1, function(y) euc.dist(p.proteome.index[[i]][x],p.learn.index[[i]][y]))
            dist.result <- as.data.table(cbind(rep(as.character(p.proteome$P_ID[x]),length(p.learn.index[[i]])), as.character(p.learn$P_ID), unlist(calculation)))
            dist.result <- dist.result[as.numeric(as.character(dist.result$V3))>cutoff,]
            p.matrix <- rbind(p.matrix, dist.result)
#            if (total < 300){
#                setxtProgressBar(pb, x)
#            }
#            if (total >= 300){
#               if (x %% 3 == 0){
            setTxtProgressBar(pb, x)
#               }
#          }
        }
        close(pb)
        p.matrix <- cbind(rep(i, nrow(p.matrix)), p.matrix)
        write.table(p.matrix, helper_i_p, quote = FALSE, sep = "\t", row.names = FALSE)
        rm(p.matrix); gc()
    } else {
#        print("Descriptor already in database")
        cat(sprintf("Profile similarity of %s in pathogen already in database \n", descriptors_vec[[i]]))
    }
}

cat(sprintf("Done!")
end_time <- Sys.time()
print(end_time - start_time)


##### Profile similarity in human proteins
start_time <- Sys.time()

for (i in 1:length(descriptors_vec)){
#    print("Calculating the human similarity profiles")

    helper_i_h <- paste(host_dir, host_org, "/descriptors/h_results_lagmax_", lagmax, "_cutoff_", cutoff_h, "_", descriptors_vec[i], ".txt", sep = "")
    if (file.exists(helper_i_h) == FALSE){
#        print("Calculating the pathogen similarity profiles")
        cat(sprintf("Calculating the host similarity profiles for descriptor %s \n", descriptors_vec[[i]]))
        h.matrix <- data.table()
#pb <- progress_bar$new(total = length(h.proteome.index[[i]]))

        total <- length(h.proteome.index[[i]])
        pb <- txtProgressBar(min = 0, max = total, style = 3)

        for(x in 1:length(h.proteome.index[[i]])) {
  #pb$tick()
            calculation <- mclapply(1:length(h.learn.index[[i]]), mc.cores=cores-1, function(y) euc.dist(h.proteome.index[[i]][x],h.learn.index[[i]][y]))
            dist.result <- as.data.table(cbind(rep(as.character(h.proteome$H_ID[x]),length(h.learn.index[[i]])), as.character(h.learn$H_ID), unlist(calculation)))
            dist.result <- dist.result[as.numeric(as.character(dist.result$V3))>cutoff,]
            h.matrix <- rbind(h.matrix, dist.result)
#            if (total < 300){
#                setxtProgressBar(pb, x)
#            }
#            if (total >= 300){
#                if (x %% 3 == 0){ 
            setTxtProgressBar(pb, x)
#               }   
#            }   
        }
        close(pb)
        h.matrix <- cbind(rep(i, nrow(h.matrix)), h.matrix)
        write.table(h.matrix, helper_i_h, quote = FALSE, sep = "\t")
        rm(h.matrix); gc()
    } else {
#        print("Descriptor already in database")
        cat(sprintf("Profile similarity of %s in host already in database \n", descriptors_vec[[i]]))
    }
}
cat(sprintf("Done!"))
end_time <- Sys.time()
print(end_time - start_time)


cat(sprintf("Starting with the analysis script"))

#system(sprintf("mkdir %s", new_dir))

system(sprintf("if [ ! -d %s ]; then mkdir %s; fi", new_dir, new_dir))



#sprintf("if [ ! -d %s ]; then mkdir %s; fi", new_dir)
# Add reference for amino acid U in the BLOSUM matrix
data("BLOSUM62")
BLOSUM62 <- rbind(BLOSUM62, rep(-5, 25))
BLOSUM62 <- cbind(BLOSUM62, c(rep(-5, 25), 1))
colnames(BLOSUM62)[26] <- "U"; rownames(BLOSUM62)[26] <- "U"

# Define alignment function to calculate the aligment score for interactome entries
align2 <- function(v1, v2){
  pairwiseAlignment(v1, v2, substitutionMatrix = BLOSUM62, gapOpening = 0, gapExtension = 8)@score
}


# "a" and "b" parameters to define lines for FPR = (0.01, 0.001, 0.0001)

false_positive <- as.numeric(args[13])


if (false_positive == 0.1)  { 
  a_p = -0.00040613
  b_p = 0.754
  a_h = -0.00034715
  b_h = 0.732
}


if (false_positive == 0.01)  { 
    a_p = -0.00040613
    b_p = 0.817
    a_h = -0.00034715
    b_h = 0.799
}

if (false_positive == 0.001)  {  
    a_p = -0.00040613
    b_p = 0.869
    a_h = -0.00034715
    b_h = 0.853
}

if (false_positive == 0.0001)  {  
    a_p = -0.00040613
    b_p = 0.905
    a_h = -0.00034715
    b_h = 0.902
}


if (false_positive == 0.00001)  {  
  a_p = -0.00040613
  b_p = 0.93
  a_h = -0.00034715
  b_h = 0.942
}



false_positive_h <-  gsub(false_positive, pattern = "\\.", replacement = "_")


#new_dir <- paste(results_dir, host_org, "_vs_", pathogen_org, "/", sep = "")


# create subdirectory for the chosen false positive rate
false_positive_dir <-paste(new_dir, "false_positive_rate_", false_positive_h, "_predictions/", sep = "")


system(sprintf("if [ ! -d %s ]; then mkdir %s; fi", false_positive_dir, false_positive_dir))

#create two new subdirectories
interactomes_dir <- paste(false_positive_dir, "interactomes/", sep = "")
intermediate_files_dir <- paste(false_positive_dir, "intermediate_files/", sep = "")

system(sprintf("if [ ! -d %s ]; then mkdir %s; fi", interactomes_dir, interactomes_dir))
system(sprintf("if [ ! -d %s ]; then mkdir %s; fi", intermediate_files_dir, intermediate_files_dir))


cat(sprintf("Creating interactome........."))
start_time <- Sys.time()

#Load and create interactomes in a loop, one for each previously descriptor


interactome_names <- list()
for (i in 1:length(descriptors_vec)){
# Load prediction results
    helper_name <- paste(interactomes_dir, "descriptor_", descriptors_vec[i], "_interactome.txt", sep = "")
    if (file.exists(helper_name) == FALSE){
        helper_i_p <- paste(pathogen_dir, pathogen_org, "/descriptors/p_results_lagmax_", lagmax, "_cutoff_", cutoff_h, "_",  descriptors_vec[i], ".txt", sep = "")  
        helper_i_h <- paste(host_dir, host_org, "/descriptors/h_results_lagmax_", lagmax, "_cutoff_", cutoff_h, "_", descriptors_vec[i], ".txt", sep = "")
        p_results <- read.delim(helper_i_p)[,-1]
        colnames(p_results) <- c("P_ID_query", "P_ID_template", "P_dist")
        h_results <- read.delim(helper_i_h)[,-1]
        colnames(h_results) <- c("H_ID_query", "H_ID_template", "H_dist")

# Load reference datasets

        help_p <- paste(pathogen_dir, pathogen_org, "/", pathogen_org, "_proteome.tsv", sep = "")
        pathogen <- read.delim(help_p, header = TRUE)
        colnames(pathogen) <- c("P_ID_query", "Pathogen", "Pathogen_taxon_ID", "P_len_query", "P_seq_query")
        help_h <- paste(host_dir, host_org, "/", host_org, "_proteome.tsv", sep = "")
        host <- read.delim(help_h, header = TRUE)
        colnames(host) <- c("H_ID_query","Host", "Host_taxon_ID", "H_len_query", "H_seq_query")


        learn_set <- read.delim(args[10])
        colnames(learn_set) <- c("Strain_template", "Taxonomy_template", "P_ID_template", "P_AC_template", "P_seq_template", 
                         "P_len_template", "H_ID_template", "H_AC_template", "H_seq_template", "H_len_template", 
                         "Method", "Pubmed_ID")

#### Generate the final interactomes, one for each descriptor


        h_results <- as_tibble(h_results)
        p_results <- as_tibble(p_results)


        p_filtered <- inner_join(p_results, pathogen, by="P_ID_query")
        p_filtered <- filter(p_filtered, P_len_query > 150) 
        p_filtered$P_dis_cutoff <- (a_p*p_filtered$P_len_query) + b_p  #This is an empiric function derived from the relation between length and similarity values
        p_filtered <- filter(p_filtered, P_dist > P_dis_cutoff)

        h_filtered <- inner_join(h_results, host, by="H_ID_query")
        h_filtered <- filter(h_filtered, H_len_query > 150)
        h_filtered$H_dis_cutoff <- (a_h*h_filtered$H_len_query)+ b_h
        h_filtered <- filter(h_filtered, H_dist > H_dis_cutoff)


# Join all data to create the interactome
        interactome <- inner_join(p_filtered, learn_set, by="P_ID_template")
        interactome <- inner_join(interactome, h_filtered, by="H_ID_template")

# Include conservation analysis
        interactome$P_align <- mapply(align2,interactome$P_seq_query, interactome$P_seq_template)
        interactome$H_align <- mapply(align2,interactome$H_seq_query, interactome$H_seq_template)
        helper_name <- paste(interactomes_dir, "descriptor_", descriptors_vec[i], "_interactome.txt", sep = "")
#        helper_name <- paste(new_dir, host_org, "_vs_", pathogen_org, "_", descriptors_vec[i], "_false_positive_", false_positive_h, "_interactome.txt", sep = "")
        interactome_length <- nrow(interactome)
        if ( interactome_length < 1) {
        cat(sprintf("No entries generated for predicted interactome with %s  \n", descriptors_vec[[i]]))
       
        }





        else {
#        interactome_names[i] <- helper_name
        interactome_names <- append(interactome_names, helper_name)
        write.table(interactome, helper_name, quote = FALSE, sep = "\t", row.names = FALSE)

        cat(sprintf("Done!"))
        end_time <- Sys.time()
        print(end_time - start_time)
        cat(sprintf("..........Congratulations, your interactome is ready for analysis!"))
        }
    }

    else     {
  
    cat(sprintf("Interactome already in database"))
#    interactome_names[i] <- helper_name
    interactome_names <- append(interactome_names, helper_name)
    }
  
}

length_interactome_names <- length(interactome_names)
cat(sprintf("NUMBER OF GENERATED INTERACTOMES IS %s  \n", length_interactome_names))



#----------------Combining interactomes

if (length_interactome_names > 1) {  
    interactome_names <- matrix(unlist(interactome_names))
#write.table(interactome_names, "interactome_names.txt", sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

    consensus_predictors_p <- as.numeric(args[14])
    consensus_predictors_p_h <- consensus_predictors_p*100
    consensus_predictors <- round(consensus_predictors_p*length_interactome_names)

    if (length_interactome_names == 2) {

        consensus_predictors_p <- as.numeric(args[14])
        consensus_predictors_p_h <- consensus_predictors_p*100
        consensus_predictors <- 2


    }


    interactomas <- list()
    for (i in 1:length(interactome_names)){
        interactomas[[i]] <- read.delim(interactome_names[i], header = TRUE)
    }

    length_interactomas <- length(interactomas)

#if ( length_interactomas > 1 ){
#Remove interactomes with 0 predicted interactions
    interactomas <- interactomas[sapply(interactomas, nrow)>0]
    interactomes_combined <- bind_rows(interactomas, .id = "descrip")
    interactomes_combined <- unique(interactomes_combined[,c(1,2,21)])
    helper <- ddply(interactomes_combined, .(H_ID_query, P_ID_query), nrow)
    colnames(helper) <- c("H_ID_query", "P_ID_query", "supporting_predictions")
    all_interactomes <- bind_rows(interactomas, .id = "descrip")
    all_interactomes <- unique(all_interactomes)
    merged_interactomes <- merge(all_interactomes, helper, by = c("H_ID_query", "P_ID_query"))
    unique_predicted_interactome <- merged_interactomes[!duplicated(merged_interactomes[c("H_ID_query","P_ID_query")]),]
#write.table(unique_predicted_interactome, "unique_predicted_interactome.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
#here ask how many supporting predictions you want to set as a threshold

    predicted_interactome_by_some_pred <- subset(unique_predicted_interactome, unique_predicted_interactome$supporting_predictions > consensus_predictors)[,-c(3,5,22,27,28)]
    helper_n <- paste(interactomes_dir, "consensus_interactome_", consensus_predictors_p_h, "_percent_of_predictors.tsv", sep = "")
#    helper_n <- paste(new_dir, host_org, "_vs_", pathogen_org, "_predicted_interactome_false_positives_", false_positive_h, "_by_", consensus_predictors_p_h, "_percent_of_predictors.tsv", sep = "")
    tempfile <- paste(results_dir, "/temp_file", sep = "")
    tempfile_subdir <- paste(results_dir, "/temp_file_subdir", sep = "")
    results_subdir <- as.data.frame(rbind(interactomes_dir, intermediate_files_dir))
    write.table(results_subdir, tempfile_subdir, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    write.table(helper_n, tempfile, row.names = FALSE, quote = FALSE, col.names = FALSE)
    if (file.exists(helper_n) == FALSE){
        write.table(predicted_interactome_by_some_pred, helper_n, sep = "\t", row.names = FALSE, quote = FALSE)
    }


} else if (length_interactome_names <= 1) {
#} else {
   cat(sprintf("hola"))
    helper_n <- "No consensus interactome generated"
    tempfile <- paste(results_dir, "/temp_file", sep = "")
    write.table(helper_n, tempfile, row.names = FALSE, quote = FALSE, col.names = FALSE)

#}


}
#else if (length_interactome_names == 0) {
#    helper_n <- "No interactome generated"
#    write.table(helper_n, tempfile, row.names = FALSE, quote = FALSE, col.names = FALSE)

#}

#system(sprintf("value=$(echo %s)", helper_n))
