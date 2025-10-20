#include <cmath>

#include <amsim/stats.h>
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
      // rework this metric to use statistics header
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
            stats::standardise(n_ind, ctx.phenotypes[pheno](type), 1, buf_cur, 1,
                               mean_cur, prec_cur);
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

      // between-mate phenotype component correlation
      MetricFunc f_comp_xcor(ComponentType type) {
        return [type](const SimulationContext &ctx) -> std::vector<double> {
          const std::size_t n_pheno               = ctx.phenotypes.size();
          const std::size_t n_sex                 = ctx.phenotypes[0].n_ind() / 2;
          const std::vector<std::size_t> matching = ctx.model.state;

          const double*     buf_male   = ctx.phenotypes[0](type);
          const double*     buf_female = buf_male + n_sex;
          const int         lda_buf    = 2 * n_sex;

          std::vector<double> cor_buf(n_pheno * n_pheno);

          std::vector<double> female_reordered(n_sex * n_pheno);
          for (std::size_t pheno = 0; pheno < n_pheno; ++pheno) {
            const double* src = buf_female + pheno * lda_buf;
            double* dst = female_reordered.data() + pheno * n_sex;
            for (std::size_t male_idx = 0; male_idx < n_sex; ++male_idx) {
              dst[male_idx] = src[matching[male_idx]];
            }
          }

          stats::cor(n_sex, n_pheno, n_pheno, buf_male, lda_buf,
                     female_reordered.data(), n_sex, cor_buf.data(), n_pheno);

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

      // phenotype heritabilities
      std::vector<double> f_pheno_h2(const SimulationContext &ctx) {
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
      std::vector<double> f_latent_h2(const SimulationContext &ctx) {
        const std::size_t n_pheno = ctx.phenotypes.size();
        const std::size_t n_sex   = ctx.phenotypes[0].n_ind() / 2;

        const std::vector<double>& U  = ctx.model.cor_U;
        const std::vector<double>& VT = ctx.model.cor_VT;
        std::vector<double> latent_h2(2 * n_pheno);

        const double* buf_gen_male   = ctx.phenotypes[0](ComponentType::GENETIC);
        const double* buf_tot_male   = ctx.phenotypes[0](ComponentType::TOTAL);
        const double* buf_gen_female = buf_gen_male + n_sex;
        const double* buf_tot_female = buf_tot_male + n_sex;
        const int     lda_buf        = 2 * n_sex;

        std::vector<double> latent_gen(n_sex);
        std::vector<double> latent_tot(n_sex);

        for (std::size_t dim = 0; dim < n_pheno; ++dim) {
          // compute male latent phenotype
          const double* U_col = U.data() + dim * n_pheno;

          cblas_dgemv(CblasColMajor, CblasNoTrans, n_sex, n_pheno,
                      1.0, buf_gen_male, lda_buf, U_col, 1,
                      0.0, latent_gen.data(), 1);

          cblas_dgemv(CblasColMajor, CblasNoTrans, n_sex, n_pheno,
                      1.0, buf_tot_male, lda_buf, U_col, 1,
                      0.0, latent_tot.data(), 1);

          latent_h2[dim] = stats::var(n_sex, latent_gen.data(), 1) /
                           stats::var(n_sex, latent_tot.data(), 1);

          // compute female latent phenotype
          const double* VT_row = VT.data() + dim;

          cblas_dgemv(CblasColMajor, CblasNoTrans, n_sex, n_pheno,
                      1.0, buf_gen_female, lda_buf, VT_row, n_pheno,
                      0.0, latent_gen.data(), 1);
          cblas_dgemv(CblasColMajor, CblasNoTrans, n_sex, n_pheno,
                      1.0, buf_tot_female, lda_buf, VT_row, n_pheno,
                      0.0, latent_tot.data(), 1);

          latent_h2[n_pheno + dim] = stats::var(n_sex, latent_gen.data(), 1) /
                                     stats::var(n_sex, latent_tot.data(), 1);
        }

        return latent_h2;
      }

      // within-individual latent phenotype component correlation
      MetricFunc f_latent_comp_cor(ComponentType type) {
        return [type](const SimulationContext &ctx) -> std::vector<double> {
          const std::size_t n_pheno = ctx.phenotypes.size();
          const std::size_t n_sex   = ctx.phenotypes[0].n_ind() / 2;

          const std::vector<double>& U  = ctx.model.cor_U;
          const std::vector<double>& VT = ctx.model.cor_VT;

          const double* buf_male   = ctx.phenotypes[0](type);
          const double* buf_female = buf_male + n_sex;
          const int     lda_buf    = 2 * n_sex;

          std::vector<double> latent_male(n_sex * n_pheno);
          std::vector<double> latent_female(n_sex * n_pheno);

          // compute male latent phenotypes
          cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                      n_sex, n_pheno, n_pheno,
                      1.0, buf_male, lda_buf,
                      U.data(), n_pheno,
                      0.0, latent_male.data(), n_sex);

          // compute female latent phenotypes
          cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                      n_sex, n_pheno, n_pheno,
                      1.0, buf_female, lda_buf,
                      VT.data(), n_pheno,
                      0.0, latent_female.data(), n_sex);

          std::vector<double> cor_buf(2 * n_pheno * n_pheno);

          stats::cor(n_sex, n_pheno, latent_male.data(), n_sex,
                     cor_buf.data(), n_pheno);

          stats::cor(n_sex, n_pheno, latent_female.data(), n_sex,
                     cor_buf.data() + n_pheno * n_pheno, n_pheno);

          return cor_buf;
        };
      }

      // between-mate latent phenotype component correlation
      MetricFunc f_latent_comp_xcor(ComponentType type) {
        return [type](const SimulationContext &ctx) -> std::vector<double> {
          const std::size_t n_pheno = ctx.phenotypes.size();
          const std::size_t n_sex   = ctx.phenotypes[0].n_ind() / 2;
          const std::vector<std::size_t> matching = ctx.model.state;

          const std::vector<double>& U  = ctx.model.cor_U;
          const std::vector<double>& VT = ctx.model.cor_VT;

          const double* buf_male   = ctx.phenotypes[0](type);
          const double* buf_female = buf_male + n_sex;
          const int     lda_buf    = 2 * n_sex;

          // Compute latent phenotypes for all dimensions
          std::vector<double> latent_male(n_sex * n_pheno);
          std::vector<double> latent_female_unordered(n_sex * n_pheno);

          // Males: project onto all columns of U
          cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                      n_sex, n_pheno, n_pheno,
                      1.0, buf_male, lda_buf,
                      U.data(), n_pheno,
                      0.0, latent_male.data(), n_sex);

          // Females: project onto all columns of V (rows of VT)
          cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                      n_sex, n_pheno, n_pheno,
                      1.0, buf_female, lda_buf,
                      VT.data(), n_pheno,
                      0.0, latent_female_unordered.data(), n_sex);

          // Reorder females according to mating
          std::vector<double> latent_female(n_sex * n_pheno);
          for (std::size_t dim = 0; dim < n_pheno; ++dim) {
            const double* src = latent_female_unordered.data() + dim * n_sex;
            double* dst = latent_female.data() + dim * n_sex;
            for (std::size_t male_idx = 0; male_idx < n_sex; ++male_idx) {
              dst[male_idx] = src[matching[male_idx]];
            }
          }

          // Compute cross-correlation matrix
          std::vector<double> cor_buf(n_pheno * n_pheno);

          stats::cor(n_sex, n_pheno, n_pheno,
                     latent_male.data(), n_sex,
                     latent_female.data(), n_sex,
                     cor_buf.data(), n_pheno);

          return cor_buf;
        };
      }
    }

    // latent phenotype component variance
    MetricFunc f_latent_comp_var(ComponentType type) {
      return [type](const SimulationContext &ctx) -> std::vector<double> {
        const std::size_t n_pheno = ctx.phenotypes.size();
        const std::size_t n_sex   = ctx.phenotypes[0].n_ind() / 2;
        const std::size_t n_ind   = 2 * n_sex;

        const std::vector<double>& U  = ctx.model.cor_U;
        const std::vector<double>& VT = ctx.model.cor_VT;

        const double* buf_male   = ctx.phenotypes[0](type);
        const double* buf_female = buf_male + n_sex;
        const int     lda_buf    = 2 * n_sex;

        std::vector<double> latent_all(n_ind * n_pheno);

        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    n_sex, n_pheno, n_pheno,
                    1.0, buf_male, lda_buf,
                    U.data(), n_pheno,
                    0.0, latent_all.data(), n_ind);

        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                    n_sex, n_pheno, n_pheno,
                    1.0, buf_female, lda_buf,
                    VT.data(), n_pheno,
                    0.0, latent_all.data() + n_sex, n_ind);

        std::vector<double> var_buf(n_pheno);

        for (std::size_t dim = 0; dim < n_pheno; ++dim) {
          const double* latent_dim = latent_all.data() + dim * n_ind;
          var_buf[dim] = stats::var(n_ind, latent_dim, 1);
        }

        return var_buf;
      };
    }

    Metric comp_var(const std::size_t n_pheno, ComponentType type) {
      std::string name = "implement_name_comp_var";
      return make_metric(phenome::f_comp_var(type), name, n_pheno, 1);
    }

    Metric pheno_h2(const std::size_t n_pheno) {
      return make_metric(phenome::f_pheno_h2, "h2", n_pheno, 1);
    }

    Metric latent_h2(const std::size_t n_pheno) {
      return make_metric(phenome::f_latent_h2, "latent_h2", 2, n_pheno);
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

    Metric comp_xcor(const std::size_t n_pheno, ComponentType type) {
      std::string name;
      switch (type) {
        case ComponentType::GENETIC:
          name = "gen_xcor";
          break;
        case ComponentType::ENVIRONMENTAL:
          name = "env_xcor";
          break;
        case ComponentType::VERTICAL:
          name = "vert_xcor";
          break;
        case ComponentType::TOTAL:
          name = "pheno_xcor";
          break;
      }
      return make_metric(phenome::f_comp_xcor(type), name, n_pheno, n_pheno);
    }
  }
}
