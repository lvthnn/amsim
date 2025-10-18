#include <cmath>
#include <iostream>

#include <amsim/metric.h>
#include <amsim/component_type.h>
#include <amsim/simulation_context.h>

#if defined(__APPLE__)
  #include <Accelerate/Accelerate.h>
#else
  #include <cblas.h>
#endif

namespace amsim {
  Metric::Metric(MetricFunc f, const std::string &name,
                 const std::size_t n_rows, const std::size_t n_cols,
                 std::optional<std::vector<std::string>> names)
    : name(name),
      n_rows(n_rows),
      n_cols(n_cols),
      named((names)),
      names((named) ? std::move(*names) : std::vector<std::string>{}),
      f_(std::move(f)) {
    buf_.resize(n_rows * n_cols);
  }

  std::string Metric::header() {
    std::string header = "it\t"; 
    for (std::size_t el = 0; el < buf_.size(); el++) {
      std::string el_name = named ? names[el] : std::to_string(el);
      header += name + "::" + el_name;
      if (el < buf_.size() - 1) header += "\t";
    }
    return header;
  }

  std::string Metric::stream(const SimulationContext &ctx) {
    buf_ = f_(ctx);
    std::string res;
    for (std::size_t el = 0; el < buf_.size(); el++)
      res += std::to_string(buf_[el]) + ((el < (buf_.size() - 1)) ? "\t" : "");
      
    return res;
  }

  namespace metrics {
    Metric make_metric(MetricFunc f, const std::string &name,
                       const std::size_t n_rows, const std::size_t n_cols,
                       std::optional<std::vector<std::string>> el_names) {
      return Metric{std::move(f), name, n_rows, n_cols, el_names};
    }

    namespace phenome {
      MetricFunc f_comp_cor(ComponentType type) {
        return [type](const SimulationContext &ctx) -> std::vector<double> {
          const std::size_t n_pheno = ctx.phenotypes.size();
          const std::size_t n_ind   = ctx.phenotypes[0].n_ind();
          std::vector<double> cor_buf(n_pheno * n_pheno);
          std::vector<double> std_buf(n_pheno * n_ind);
          const std::vector<double> ones(n_ind, 1.0);

          for (std::size_t pheno = 0; pheno < n_pheno; pheno++) {
            if (ctx.phenotypes[pheno].comp_var(type) == 0)
              continue;
            double* buf_cur = &std_buf[pheno * n_ind];
            double mean_cur = ctx.phenotypes[pheno].comp_mean(type);
            double prec_cur = 1.0 / std::sqrt(ctx.phenotypes[pheno].comp_var(type));

            cblas_dcopy(n_ind, ctx.phenotypes[pheno](type), 1, buf_cur, 1);
            cblas_daxpy(n_ind, -mean_cur, ones.data(), 1, buf_cur, 1);
            cblas_dscal(n_ind, prec_cur, buf_cur, 1);
          }

          for (std::size_t i = 0; i < n_pheno; i++) {
            const double* ptr_i = &std_buf[i * n_ind];
            for (std::size_t j = i; j < n_pheno; j++) {
              if (i == j) {
                cor_buf[i * n_pheno + j] = 1.0;
                continue;
              }
              const double* ptr_j = &std_buf[j * n_ind];
              cor_buf[i * n_pheno + j] = (1.0 / static_cast<double>(n_ind)) *
                                          cblas_ddot(n_ind, ptr_i, 1, ptr_j, 1);
              cor_buf[j * n_pheno + i] = cor_buf[i * n_pheno + j];
            }
          }

          return cor_buf;
        };
      }

      MetricFunc f_comp_var(ComponentType type) {
        return [type](const SimulationContext &ctx) -> std::vector<double> {
          std::vector<double> comp_vars(ctx.phenotypes.size());
          for (std::size_t pheno = 0; pheno < ctx.phenotypes.size(); pheno++)
            comp_vars[pheno] = ctx.phenotypes[pheno].comp_var(type);
          return comp_vars;
        };
      }

      // between-mate phenotype component correlation
      
      // phenotype heritabilities
      std::vector<double> f_h2(const SimulationContext &ctx) {
        const std::size_t n_pheno = ctx.phenotypes.size();
        std::vector<double> pheno_h2(n_pheno);

        for (std::size_t pheno = 0; pheno < n_pheno; pheno++) {
          const Phenotype& pheno_cur = ctx.phenotypes[pheno];
          pheno_h2[pheno] = pheno_cur.comp_var(ComponentType::GENETIC) /
                            pheno_cur.comp_var(ComponentType::TOTAL);
        }

        return pheno_h2;
      }


      // latent phenotype heritability

      // within-individual latent phenotype component correlation

      // between-mate latent phenotype component correlation
    }

    Metric comp_cor(const std::size_t n_pheno, ComponentType type) {
      std::string name;
      switch (type) {
        case ComponentType::GENETIC:
          name = "gen_cor";
          break;
        case ComponentType::ENVIRONMENTAL:
          name = "env_cor";
          break;
        case ComponentType::VERTICAL:
          name = "vert_cor";
          break;
        case ComponentType::TOTAL:
          name = "pheno_cor";
          break;
      }

      return make_metric(phenome::f_comp_cor(type), name, n_pheno, n_pheno);
    }

    Metric comp_var(const std::size_t n_pheno, ComponentType type) {
      std::string name = "i_hate_this";
      return make_metric(phenome::f_comp_var(type), name, n_pheno, 1);
    }

    Metric pheno_h2(const std::size_t n_pheno) {
      return make_metric(phenome::f_h2, "h2", n_pheno, 1);
    }
  }
}
