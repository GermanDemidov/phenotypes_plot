library(ontologyIndex)

data_for_plotting <- read.table("/mnt/users/ahdemig1/forBenita/pheno_plot_22_07_27/outfile_annotated2_cp.txt", header=T, sep="\t", comment.char="&", fill=T)

NGSD_database <- read.table("/mnt/users/ahdemig1/forJoo/NGSD150722_only_research.txt", header=T, comment.char = "?", sep="\t", fill = TRUE, quote = "\"")

hpo <- get_ontology("/mnt/users/ahdemig1/forBenita/height_of_hpos/hp.obo")

blacklist_hpos = c("HP:0000006", "HP:0000007", "HP:0001417", "HP:0001419", "HP:0001423", "HP:0001428", "HP:0001450", "HP:0003593")

only_phenotypes_hpo = get_descendants(hpo, "HP:0000118", exclude_roots = FALSE)

set_of_branches = c("HP:0001257", "HP:0001332", "HP:0001251")
set_of_hpos_to_select <- c()
table_to_output <- matrix(nrow=0, ncol=2)
for (elem in set_of_branches) {
  name_id = which(hpo$id == elem)
  name = hpo$name[name_id]
  set_of_hpos_to_select <- c(set_of_hpos_to_select, get_descendants(hpo, elem, exclude_roots = FALSE))
  for (hpo_elem in get_descendants(hpo, elem, exclude_roots = FALSE)) {
    table_to_output <- rbind(table_to_output, c(hpo_elem, name))
  }
}
write.table(table_to_output, file="./HPOs.disease.tsv", sep="\t", row.names = F, quote=F, col.names = F)