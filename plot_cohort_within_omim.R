### Master ontology created

library(ontologyIndex)
library(ontologySimilarity)
library(Rcpp)
library(sets)
library(ggplot2)
library(umap)
library(ggrepel)

# THE NEXT BLOCK IS A CONFIG, CHANGE IT TO YOUR PATHS
# Path to hp.obo from HPO
hpo <- get_ontology("/path_to_downloaded_files/hp.obo")
# Path to the file with the cohort
data_for_plotting <- read.table("/path_to_database/database.tsv", header=T, sep="\t", comment.char="&", fill=T, quote="\"")
# A column with the sample IDs
index_of_IDs = which(colnames(data_for_plotting) == "X.name")
# A column with the 
index_of_hpos = which(colnames(data_for_plotting) == "disease_details_HPO_term_id")
# HPOs to disease cathegory - COMMENT if you use manual annotated disease groups
disease_cathegory_doc <- read.table("/path_to_generated_with_create_hpo_list_for_annotationR/HPOs.disease.tsv", comment.char="?")
number_of_clusters = 60






# CODE - DON'T CHANGE ANYTHING IF YOU JUST WANT A PLOT
blacklist_hpos = c("HP:0000006", "HP:0000007", "HP:0001417", "HP:0001419", "HP:0001423", "HP:0001428", "HP:0001450", "HP:0003593")
gene_to_pheno <- read.table("/path_to_downloaded_files/genes_to_phenotype.txt", sep="\t", quote="")
gene_to_pheno = gene_to_pheno[-which(gene_to_pheno[,3] %in% blacklist_hpos), ]

only_phenotypes_hpo = get_descendants(hpo, "HP:0000118", exclude_roots = FALSE)


patients_to_show <- unique(data_for_plotting[,index_of_IDs])


index_of_manual_annotation = which(colnames(data_for_plotting) == "disease_group")
annotate_cathegory_manually = function(name) {
  disease_for_annotation = data_for_plotting[which(data_for_plotting[,index_of_IDs] == name),index_of_manual_annotation]
  if (disease_for_annotation == "n/a") {
    disease_for_annotation = "Other"
  }
  return(disease_for_annotation)
}

annotate_cathegory = function(HPOs) {
  if (length(intersect(HPOs, unique(disease_cathegory_doc[,1]))) > 0) {
    diseases <- sort(unique(disease_cathegory_doc[which(disease_cathegory_doc[,1] %in% HPOs), 2]))
  } else {
    diseases = "Other"
  }
  return(paste0(diseases, collapse="_"))
}

list_of_phenotypes_patients = list()
list_of_genes = list()
number_of_patients = list()
disease_category = list()

for (name in patients_to_show) {
  hpo_numbers = data_for_plotting[which(data_for_plotting$X.name == name), index_of_hpos]
  hpos_list <- strsplit(hpo_numbers, "; ")[[1]]
  hpo_numbers = c()
  for (elem in hpos_list) {
    hpo_num = strsplit(elem, " - ")[[1]][1]
    if (hpo_num %in% hpo$id) {
      hpo_numbers = c(hpo_numbers, hpo_num)
    }
  }
  
  
  hpo_numbers = hpo_numbers[which(hpo_numbers %in% hpo$id & !hpo_numbers %in% blacklist_hpos & hpo_numbers %in% only_phenotypes_hpo)]
  hpo_numbers = sort(hpo_numbers)
  
  disease_cathegory_current = annotate_cathegory_manually(name)
  
  disease_cathegory_current = annotate_cathegory(hpo_numbers)
  
  hpo_codes = paste0(hpo_numbers, collapse="|")
  # comment IF if all the patients needed
  if (disease_cathegory_current != "Other" & !hpo_codes %in% names(number_of_patients)) {
    disease_category[[name]] = disease_cathegory_current
    list_of_phenotypes_patients[[name]] = hpo_numbers
    number_of_patients[[hpo_codes]] = 1
  } else if (hpo_codes %in% names(number_of_patients)) {
    number_of_patients[[hpo_codes]] = number_of_patients[[hpo_codes]] + 1
  }
}




overall_disease = c()
for (gene in unique(gene_to_pheno[,2])) {
  diseases = unique(gene_to_pheno[gene_to_pheno[,2] == gene, 9])
  overall_disease = c(overall_disease, diseases)
}
overall_disease = unique(overall_disease)

ids_to_coords_on_the_plot <- matrix(0, ncol=2)

hpos_unique_for_disease <- c()
list_of_phenotypes_for_plotting <- list()
for (disease in overall_disease) {
  set_of_hpos <- gene_to_pheno[which(gene_to_pheno[,9] == disease & !gene_to_pheno[,6] %in% c("HP:0040284") & gene_to_pheno[,3] %in% only_phenotypes_hpo ),3]
  if (length(set_of_hpos) == 0) set_of_hpos <- gene_to_pheno[which(gene_to_pheno[,9] == disease  & gene_to_pheno[,3] %in% only_phenotypes_hpo ),3]
  if (length(set_of_hpos) > 0) {
    set_of_hpos = set_of_hpos[set_of_hpos %in% hpo$id & set_of_hpos %in% only_phenotypes_hpo]
    sorted_hpos <- sort(set_of_hpos)
    if (!paste0(sorted_hpos, collapse = "_") %in% hpos_unique_for_disease) {
      list_of_phenotypes_for_plotting[[disease]] = set_of_hpos
      hpos_unique_for_disease <- c(hpos_unique_for_disease, paste0(sorted_hpos, collapse = "_"))
    }
    ids_to_coords_on_the_plot <- rbind(ids_to_coords_on_the_plot, c(disease, paste0(sorted_hpos, collapse = "_")))
  } else {
    break
  }
}
num_of_diseases = nrow(ids_to_coords_on_the_plot)
num_of_diseases_dots <- length(list_of_phenotypes_for_plotting)

hpos_unique_for_dot <- c()
for (name in names(list_of_phenotypes_patients)) {
  set_of_hpos <- list_of_phenotypes_patients[[name]]
  sorted_hpos <- sort(set_of_hpos)
  if (!paste0(sorted_hpos, collapse = "_") %in% hpos_unique_for_dot & !paste0(sorted_hpos, collapse = "_") %in% hpos_unique_for_disease) {
    list_of_phenotypes_for_plotting[[name]] = set_of_hpos
    hpos_unique_for_dot <- c(hpos_unique_for_dot, paste0(sorted_hpos, collapse = "_"))
  }
  ids_to_coords_on_the_plot <- rbind(ids_to_coords_on_the_plot, c(name, paste0(sorted_hpos, collapse = "_")))
}
num_of_patients = nrow(ids_to_coords_on_the_plot) - num_of_diseases
num_of_patients_dots <- length(list_of_phenotypes_for_plotting) - num_of_diseases_dots



information_content <- descendants_IC(hpo)
master_sim_mat_tmp <- get_sim_grid(ontology=hpo, information_content=information_content,  term_sets=list_of_phenotypes_for_plotting)

master_sim_mat = master_sim_mat_tmp



set.seed(2)
custom.settings = umap.defaults
custom.settings$input = "dist"
custom.settings$n_components = 2
custom.settings$n_neighbors = 25
custom.settings$min_dist = 0.1
custom.settings$spread = 4




res_umap <- umap(as.matrix(max(master_sim_mat) - master_sim_mat), config=custom.settings)

coords_for_ids <- data.frame()
for (i in 1:nrow(ids_to_coords_on_the_plot)) {
  id_to_check = min(which(ids_to_coords_on_the_plot[,2] == ids_to_coords_on_the_plot[i,2]))
  coords <- res_umap$layout[which(row.names(res_umap$layout) == ids_to_coords_on_the_plot[id_to_check,1]),]
  patient_cathegory = "Disease"
  if (!ids_to_coords_on_the_plot[id_to_check,1] %in% names(disease_category)) {
    disease_or_pateint = "Disease"
  } else {
    disease_or_pateint = "Patient"
    patient_cathegory = disease_category[[ids_to_coords_on_the_plot[id_to_check,1]]]
  }
  coords_for_ids = rbind(coords_for_ids, c(ids_to_coords_on_the_plot[id_to_check,1], coords[1], coords[2], disease_or_pateint, patient_cathegory))
}
coords_for_ids = coords_for_ids[-1,]
colnames(coords_for_ids) = c("ID", "x", "y", "Patient_or_disease", "Disease_cathegory")
coords_for_ids$x = as.numeric(coords_for_ids$x)
coords_for_ids$y = as.numeric(coords_for_ids$y)











result_only_disease = coords_for_ids[which(coords_for_ids$Disease_cathegory == "Disease"),]
result_without_disease = coords_for_ids[which(coords_for_ids$Disease_cathegory != "Disease"),]

set.seed(1)
clusters = kmeans(result_only_disease[,2:3], centers=number_of_clusters, iter.max = 5000)

combinedMatrix = rbind(clusters$centers,result_without_disease[,2:3])
distances = as.matrix(dist(combinedMatrix))

number_of_hpos = 6
matrix_of_clusters = matrix(0, nrow=0, ncol=number_of_hpos + 1)
for (cl in 1:max(clusters$cluster)) {
  number_of_samples = length(which(clusters$cluster == cl))
  disease_name = result_only_disease[which(clusters$cluster == cl), 1]
  disease_HPO = gene_to_pheno[which(gene_to_pheno[,9] %in% disease_name), c(4,9)]
  disease_HPO = unique(disease_HPO)
  hpos_tab = disease_HPO[,1]
  tbl = table(hpos_tab)
  tbl = tbl[which(tbl >= sort(tbl, decreasing = T)[number_of_hpos])]
  tbl = sort(tbl, decreasing = T)[1:(number_of_hpos - 1)]
  number_of_hits = sapply(1:nrow(result_without_disease), function(i) {
    which.min(distances[1:number_of_clusters,number_of_clusters + i])
  })
  vector_of_hpos_to_write_out = paste0(100 * round(tbl / length(unique(disease_name)), 2), "% ", names(tbl))
  matrix_of_clusters = rbind(matrix_of_clusters, c(cl, vector_of_hpos_to_write_out, sum(number_of_hits == cl)))
}

colnames(matrix_of_clusters) = c("Cluster","1st common HPO", "2nd common HPO", "3rd common HPO", "4th common HPO", "5th common HPO", "Number of patients in the cluster")
write.table(file="./images/clusters.txt", matrix_of_clusters, quote=F, row.names=F, col.names=T, sep="\t")




clusters_for_diseases <- read.table("./images/clusters.txt", sep="\t", header=T)

categories = c()
for (i in 1:nrow(clusters_for_diseases)) {
  categories <- c(categories, clusters_for_diseases$X1st.common.HPO[i])
}
cluster_labels = cbind.data.frame(clusters$centers, categories)


p <- ggplot() + theme_bw() + theme(legend.position="bottom") +
  geom_point(data = subset(coords_for_ids, Disease_cathegory == 'Disease'),
             aes(x = x, y = y), alpha=0.9, color="lightgrey") +
  geom_point(data = subset(coords_for_ids, Disease_cathegory != 'Disease'),
             aes(x = x, y = y, color = Disease_cathegory), alpha=0.5, size=2) +  theme(
               legend.title = element_text( size = 8),
               legend.text = element_text(size = 8)
             ) +
  geom_label_repel(data=cluster_labels, aes(x=x, y=y, label=categories), size=3, color="black", show.legend = FALSE) + 
  theme(text=element_text(family="Helvetica", size=12)) + xlab("Umap 1") + ylab("Umap 2")



svg(file="./images/cohort_with_cluster_labels.svg", height=12, width=12)
print(p)
dev.off()


p <- ggplot() + theme_bw() + theme(legend.position="bottom") +
  geom_point(data = subset(coords_for_ids, Disease_cathegory == 'Disease'),
             aes(x = x, y = y), alpha=0.9, color="lightgrey") +
  geom_point(data = subset(coords_for_ids, Disease_cathegory != 'Disease'),
             aes(x = x, y = y, color = Disease_cathegory), alpha=0.5, size=2) +  theme(
               legend.title = element_text( size = 8),
               legend.text = element_text(size = 8)
             )  + 
  theme(text=element_text(family="Helvetica", size=12)) + xlab("Umap 1") + ylab("Umap 2")



svg(file="./images/cohort_without_cluster_labels.svg", height=12, width=12)
print(p)
dev.off()


