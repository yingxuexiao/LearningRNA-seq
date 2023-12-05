if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("sccomp")

library(sccomp)

single_cell_object |>
  sccomp_glm( 
    formula_composition = ~ type, 
    .sample =  sample, 
    .cell_group = cell_group, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  )

seurat_obj |>
  sccomp_glm( 
    formula_composition = ~ 0 + type, 
    contrasts =  c("typecancer - typehealthy", "typehealthy - typecancer"),
    .sample = sample,
    .cell_group = cell_group, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  )


seurat_obj |>
  sccomp_glm( 
    formula_composition = ~ 0 + type, 
    contrasts =  c("typecancer - typehealthy", "typehealthy - typecancer"),
    .sample = sample,
    .cell_group = cell_group, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  )


library(loo)
# Fit first model
model_with_factor_association = 
  seurat_obj |>
  sccomp_glm( 
    formula_composition = ~ type, 
    .sample =  sample, 
    .cell_group = cell_group, 
    check_outliers = FALSE, 
    bimodal_mean_variability_association = TRUE,
    cores = 1, 
    enable_loo = TRUE
  )


# Fit second model
model_without_association = 
  seurat_obj |>
  sccomp_glm( 
    formula_composition = ~ 1, 
    .sample =  sample, 
    .cell_group = cell_group, 
    check_outliers = FALSE, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 , 
    enable_loo = TRUE
  )


# Compare models
loo_compare(
  model_with_factor_association |> attr("fit") |> loo(),
  model_without_association |> attr("fit") |> loo()
)

###differential variability,binary factor
# model the cell-group variability also dependent on the type, 
# and so test differences in variability.

res = 
  seurat_obj |>
  sccomp_glm( 
    formula_composition = ~ type, 
    formula_variability = ~ type,
    .sample = sample,
    .cell_group = cell_group,
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  )

res

##Visualisation
plots = plot_summary(res) 

plots$boxplot

plots$credible_intervals_1D

##Visualisation of the MCMC chains from the posterior distribution
library(magrittr)

res %>% attr("fit") %>% rstan::traceplot("beta[2,1]")

plots = plot_summary(res)
plots$credible_intervals_1D


