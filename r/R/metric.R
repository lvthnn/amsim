.Metric <- R6::R6Class(
  classname = "Metric",

  public = list(
    initialize = function(build_fn, params = NULL) {
      private$.build_fn <- build_fn
      private$.params   <- if (!is.null(params)) params else list()
    },

    build = function() {
      if (length(private$.params) == 0)
        return(private$.build_fn())
      else
        return(do.call(private$.params, private$.params))
    }
  ),

  private = list(
    .build_fn = NULL,
    .params   = NULL
  )
)

build_metrics <- function(metrics) {
  return(lapply(metrics, function(metric) metric$build()))
}

#' @export
pheno_h2 <- function() return(.Metric$new(build_fn = .pheno_h2))

#' @export
pheno_comp_mean <- function(component_type) {
  .Metric$new(
    build_fn = .pheno_comp_mean,
    params = list(component_type = component_type)
  )
}

#' @export
pheno_comp_var <- function(component_type) {
  .Metric$new(
    build_fn = .pheno_comp_var,
    params = list(component_type = component_type)
  )
}

#' @export
pheno_comp_cor <- function(component_type) {
  .Metric$new(
    build_fn = .pheno_comp_cor,
    params = list(component_type = component_type)
  )
}

#' @export
pheno_comp_xcor <- function(component_type) {
  .Metric$new(
    build_fn = .pheno_comp_xcor,
    params = list(component_type = component_type)
  )
}

#' @export
pheno_latent_h2 <- function() return(.Metric$new(build_fn = .pheno_h2))

#' @export
pheno_latent_comp_mean <- function(component_type) {
  .Metric$new(
    build_fn = .pheno_comp_mean,
    params = list(component_type = component_type)
  )
}

#' @export
pheno_latent_comp_var <- function(component_type) {
  .Metric$new(
    build_fn = .pheno_comp_var,
    params = list(component_type = component_type)
  )
}

#' @export
pheno_latent_comp_cor <- function(component_type) {
  .Metric$new(
    build_fn = .pheno_comp_cor,
    params = list(component_type = component_type)
  )
}

#' @export
pheno_latent_comp_xcor <- function(component_type) {
  .Metric$new(
    build_fn = .pheno_comp_xcor,
    params = list(component_type = component_type)
  )
}
