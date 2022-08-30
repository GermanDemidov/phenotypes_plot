# phenotypes_plot

Scripts for making nice plots of cohorts of patients based on HPO terms

# Examples

Grey dots = diseases from OMIM or ORPHA ( http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt ). Colored dots = patients from your cohort. 

## Cohort plot based on predefined disease group annotation

## Cohort plot based on automatically generated with HPO terms disease group annotation

# What do you need

Install the following packages for you R:

```R
install.packages("ontologyIndex")
install.packages("ontologySimilarity")
install.packages("Rcpp")
install.packages("sets")
install.packages("ggplot2")
install.packages("umap")
install.packages("ggrepel")
```

Download Human Phenotype Ontology obo file: https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo .

Download genes_to_phenotype from https://hpo.jax.org/app/download/annotation .

You should have a text file with your cohort, with one column denoting unique sample IDs, and another column providing HPOs per patient, as "HP:0000238 - Hydrocephalus; HP:0000505 - Visual impairment". If your HPOs are provided differently, you have to modify the following code:

```
  hpo_numbers = data_for_plotting[which(data_for_plotting$X.name == name), index_of_hpos]
  hpos_list <- strsplit(hpo_numbers, "; ")[[1]]
  hpo_numbers = c()
  for (elem in hpos_list) {
    hpo_num = strsplit(elem, " - ")[[1]][1]
    if (hpo_num %in% hpo$id) {
      hpo_numbers = c(hpo_numbers, hpo_num)
    }
  }
```

# How to run

Modify your `plot_cohort_within_omim.R` file (the CONFIG block, before the CODE block) and put your paths there. Specify which column name indicates sample names (`index_of_IDs = which(colnames(data_for_plotting) == "X.name")`) and which column name is for HPO terms (`index_of_hpos`).

If you have manually annotated cathegories (shown as color dots on the final plot), comment `disease_cathegory_doc` and use function `annotate_cathegory_manually`. Otherwise, you need to run `create_hpo_list_for_annotation.R` prior to `plot_cohort_within_omim.R`.

In `create_hpo_list_for_annotation.R` you need to specify HPOs which define your disease group (e.g., example, ataxia, spasticity and dystonia would be `set_of_branches = c("HP:0001257", "HP:0001332", "HP:0001251")`. 

Then you generate `HPOs.disease.tsv` file - use it in `plot_cohort_within_omim.R` as `disease_cathegory_doc`.

Then simply run `plot_cohort_within_omim.R` as `Rscript plot_cohort_within_omim.R` or line by line in `R`!

Since the plot is based on all OMIM diseases, it takes several (2-4 hours) to produce the plots. It makes sense to execute the code line by line.