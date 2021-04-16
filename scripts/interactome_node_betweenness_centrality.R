#Download from STRING:


      # 1). protein network data (scored links between proteins)
## wget https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz

      # 2). st of STRING proteins incl. their display names and descriptions
## wget https://stringdb-static.org/download/protein.aliases.v11.0/9606.protein.aliases.v11.0.txt.gz

# Filter 1), combined score > 900

## awk '{ if($3 >= 900) {print $1, $2} }' /home/rocio/Javi/prediction_tool/version_sept_2020_marc/version_definitiva/toy_example_directory_2/databases/9606.protein.links.v11.0.txt > /home/rocio/Javi/prediction_tool/version_sept_2020_marc/version_definitiva/toy_example_directory_2/databases/protein_interactions_ensp
## grep 'BLAST_UniProt_AC Ensembl_HGNC_UniProt_ID(supplied_by_UniProt) Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)_AC Ensembl_UniProt Ensembl_UniProt_AC' /home/rocio/Javi/prediction_tool/version_sept_2020_marc/version_definitiva/toy_example_directory_2/databases/9606.protein.aliases.v11.0.txt | awk -F '\t' -v OFS='\t' '{print $1, $2}' > /home/rocio/Javi/prediction_tool/version_sept_2020_marc/version_definitiva/toy_example_directory_2/databases/protein_aliases

pkg <- c("plyr", "igraph")
# Install packages not yet installed
installed_packages <- pkg %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(pkg[!installed_packages])
}


invisible(lapply(pkg, library, character.only = TRUE))

args <- commandArgs(trailingOnly = TRUE)





protein_interactions <- read.csv(args[1], sep="")
protein_info <- read.table(args[2], sep="\t", header=FALSE, comment.char="#", na.strings=".", stringsAsFactors=FALSE, quote="", fill=FALSE)

prot_names1 <- protein_info
prot_names2 <- protein_info

colnames(prot_names1) <- c("protein1", "new_protein1")
colnames(prot_names2) <- c("protein2", "new_protein2")

prot_helper <- join(protein_interactions, prot_names1, by = "protein1")
prot_helper2 <- join(prot_helper, prot_names2, by = "protein2")
prot_inter_over900 <- prot_helper2[, c("new_protein1", "new_protein2")]
prot_inter_over900 <- prot_inter_over900[complete.cases(prot_inter_over900), ]



#REMOVE NAS USING COMPLETE CASES
#INTERLEAVE THE TWO COLUMNS INTO ONE 

prot_inter_over_900 <- prot_inter_over_900[complete.cases(prot_inter_over_900),]
data <- c(rbind(prot_inter_over_900$new_protein1, prot_inter_over_900$new_protein2))

g4 <- graph( data, directed = FALSE)
degree.cent <- centr_degree(g4, mode = "all")
bbb <- betweenness(g4)
hi_centrality <- cbind.data.frame(unique(data), degree.cent$res, bbb)
hi_centrality$norm_degree <- (hi_centrality$`degree.cent$res` - min(hi_centrality$`degree.cent$res`))/(max(hi_centrality$`degree.cent$res`) - min(hi_centrality$`degree.cent$res`))
#hi_centrality$norm_betweenness <- hi_centrality$bbb/max(hi_centrality$bbb)
hi_centrality$norm_betweenness <- (hi_centrality$bbb - min(hi_centrality$bbb))/(max(hi_centrality$bbb) - min(hi_centrality$bbb))
rownames(hi_centrality) <- c()
colnames(hi_centrality) <- c("H_ID_query", "degree_centrality", "betweenness_centrality", "norm_degree_centrality", "norm_betweenness_centrality")
write.table(hi_centrality, args[3])
