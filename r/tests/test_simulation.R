devtools::load_all("r")

# Instantiate the builder
sim <- Simulation$new()

# Configure simulation parameters
sim$simulation(
  n_generations = 10,
  n_individuals = 1000,
  output_dir = "output/",
  random_seed = 123
)

# Configure genome parameters
n_loci <- 4000
sim$genome(
  n_loci = n_loci,
  locus_mafs = runif(n_loci, 0.1, 0.5),
  locus_recombination = runif(n_loci, 0, 0.1),
  locus_mutation = rep(0, n_loci)
)

# Configure phenome parameters
n_pheno <- 2
sim$phenome(
  n_pheno = n_pheno,
  names = c("height", "weight"),
  loci = c(2000, 2000),
  h2_genetic = c(0.5, 0.5),
  h2_environmental = c(0.5, 0.5),
  h2_vertical = c(0.0, 0.0),
  genetic_cor = diag(n_pheno),
  environmental_cor = diag(n_pheno)
)

# Configure mating parameters
sim$mating(
  mating_type = "assortative",
  n_iterations = 1000,
  temp_init = 1.0,
  temp_decay = 0.95,
  mate_cor = diag(n_pheno)
)

# Configure metrics
sim$metrics(
  metrics = list(pheno_h2())
)

# Display the configuration summary
sim$run()
