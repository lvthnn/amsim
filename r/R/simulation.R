simulate_replication <- function(params, rep_id = NULL) {
  if (!is.null(rep_id)) {
    rep_dir <- fs::path(
      params[["simulation"]][["output_dir"]],
      format("rep_%03d", rep_id)
    )
    random_seed <- params[["simulation"]][["random_seed"]] + rep_id
  } else {
    rep_dir <- params[["simulation"]][["output_dir"]]
    random_seed <- params[["simulation"]][["random_seed"]]
  }
  .builder <- .builder_new()
  .builder_simulation(
    .builder,
    n_generations = params[["simulation"]][["n_generations"]],
    n_individuals = params[["simulation"]][["n_individuals"]],
    output_dir = rep_dir,
    random_seed = random_seed
  )
}


#' @export
Simulation <- R6::R6Class(
  classname = "Simulation",

  public = list(

    simulation = function(
      n_generations,
      n_individuals,
      output_dir,
      random_seed = NULL
    ) {
      checkmate::assert_count(n_generations)
      checkmate::assert_count(n_individuals)
      checkmate::assert_string(output_dir)
      checkmate::assert_count(random_seed)

      private$.params[["simulation"]] = list(
        n_generations = n_generations,
        n_individuals = n_individuals,
        output_dir = output_dir,
        random_seed = random_seed
      )

      invisible(self)
    },

    genome = function(
      n_loci,
      locus_mafs,
      locus_recombination,
      locus_mutation
    ) {
      checkmate::assert_count(n_loci)
      checkmate::assert_numeric(
        locus_mafs,
        len = n_loci,
        lower = 0,
        upper = 1,
        any.missing = FALSE,
        all.missing = FALSE
)
      checkmate::assert_numeric(
        locus_recombination,
        len = n_loci,
        lower = 0,
        upper = 1,
        any.missing = FALSE,
        all.missing = FALSE
      )
      checkmate::assert_numeric(
        locus_mutation,
        len = n_loci,
        lower = 0,
        upper = 1,
        any.missing = FALSE,
        all.missing = FALSE
      )

      private$.params[["genome"]] = list(
        n_loci = n_loci,
        locus_mafs = locus_mafs,
        locus_recombination = locus_recombination,
        locus_mutation = locus_mutation
      )

      invisible(self)
    },

    phenome = function(
      n_phenotypes,
      names,
      loci,
      h2_genetic,
      h2_environmental,
      h2_vertical,
      genetic_cor,
      environmental_cor
    ) {
      checkmate::assert_count(n_phenotypes)
      checkmate::assert_character(
        names,
        min.chars = 1,
        len = n_phenotypes,
        any.missing = FALSE,
        all.missing = FALSE
      )
      checkmate::assert_integerish(
        loci,
        lower = 0,
        upper = private$.params[["genome"]][["n_loci"]],
        len = n_phenotypes,
        any.missing = FALSE,
        all.missing = FALSE
      )
      assert_probability(h2_genetic, len = n_phenotypes)
      assert_probability(h2_environmental, len = n_phenotypes)
      assert_probability(h2_vertical, len = n_phenotypes)
      assert_correlation_matrix(genetic_cor)
      assert_correlation_matrix(environmental_cor)

      private$.params[["phenome"]] <- list(
        n_phenotypes = n_phenotypes,
        names = names,
        loci = loci,
        h2_genetic = h2_genetic,
        h2_environmental = h2_environmental,
        h2_vertical = h2_vertical,
        genetic_cor = as.vector(genetic_cor),
        environmental_cor = as.vector(environmental_cor)
      )

      invisible(self)
    },

    mating = function(
      mating_type,
      n_iterations = NULL,
      temp_init = NULL,
      temp_decay = NULL,
      mate_cor = NULL
    ) {
      checkmate::assert_count(n_iterations)
      checkmate::assert_numeric(temp_init, lower = 0)
      checkmate::assert_numeric(temp_decay, lower = 0, upper = 1)
      assert_correlation_matrix(mate_cor, cross = TRUE)

      private$.params[["mating"]] <- list(
        mating_type = mating_type,
        n_iterations = n_iterations,
        temp_init = temp_init,
        temp_decay = temp_decay,
        mate_cor = as.vector(mate_cor)
      )

      invisible(self)
    },

    metrics = function(metrics) {
      private$.params[["mating"]] <- list(
        metrics = metrics
      )

      invisible(self)
    }
  ),

  private = list(
    .params = list()
  ),

  active = list(
    params = function() return(private$.params)
  )
)
