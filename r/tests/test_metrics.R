# test all of the metrics to see whether they build!

devtools::load_all()

metrics <- list(
  metric_pheno_h2               = pheno_h2(),
  metric_pheno_comp_mean        = pheno_comp_mean("genetic"),
  metric_pheno_comp_var         = pheno_comp_var("genetic"),
  metric_pheno_comp_cor         = pheno_comp_cor("genetic"),
  metric_pheno_comp_xcor        = pheno_comp_xcor("genetic"),
  metric_latent_pheno_h2        = pheno_latent_h2(),
  metric_latent_pheno_comp_mean = pheno_latent_comp_mean("genetic"),
  metric_latent_pheno_comp_var  = pheno_latent_comp_var("genetic"),
  metric_latent_pheno_comp_cor  = pheno_latent_comp_cor("genetic"),
  metric_latent_pheno_comp_xcor = pheno_latent_comp_xcor("genetic")
)

for (metric in metrics) {
  build_metric(metric)
  print(paste0("Managed to build ", metric$name))
}
