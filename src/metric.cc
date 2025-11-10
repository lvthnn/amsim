#include <cmath>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

#include <amsim/component_type.h>
#include <amsim/metric.h>
#include <amsim/simulation_context.h>
#include <amsim/stats.h>

namespace amsim {
Metric::Metric(
    MetricFunc f,
    const std::string& name,
    const std::size_t n_rows,
    const std::size_t n_cols,
    std::vector<std::string> labels,
    const bool require_lat)
    : name(name),
      n_rows(n_rows),
      n_cols(n_cols),
      require_lat(require_lat),
      f_(std::move(f)),
      buf_(n_rows * n_cols),
      labels_(std::move(labels)) {}

std::string Metric::header() {
  std::string header = "gen\t";
  for (std::size_t el = 0; el < buf_.size(); ++el) {
    header += (!labels_.empty()) ? labels_[el] : std::to_string(el);
    if (el < buf_.size() - 1) header += "\t";
  }
  return header;
}

std::string Metric::stream(const SimulationContext& ctx) {
  buf_ = f_(ctx);
  std::string res;
  for (std::size_t el = 0; el < buf_.size(); ++el)
    res += std::to_string(buf_[el]) + ((el < (buf_.size() - 1)) ? "\t" : "");
  return res;
}

namespace metrics {
Metric make_metric(
    MetricFunc f,
    const std::string& name,
    const std::size_t n_rows,
    const std::size_t n_cols,
    const std::vector<std::string> labels,
    const bool require_lat) {
  return Metric{std::move(f), name, n_rows, n_cols, labels, require_lat};
}

namespace phenome {
std::vector<double> f_pheno_h2(const SimulationContext& ctx) {
  const std::size_t n_pheno = ctx.n_pheno;
  std::vector<double> pheno_h2(n_pheno);

  for (std::size_t pheno = 0; pheno < n_pheno; ++pheno) {
    const Phenotype& pheno_cur = ctx.phenotypes[pheno];
    pheno_h2[pheno] = pheno_cur.comp_var(ComponentType::GENETIC) /
                      pheno_cur.comp_var(ComponentType::TOTAL);
  }

  return pheno_h2;
}

MetricFunc f_comp_cor(ComponentType type) {
  return [type](const SimulationContext& ctx) -> std::vector<double> {
    const std::size_t n_pheno = ctx.n_pheno;
    const std::size_t n_ind = ctx.n_ind;
    std::vector<double> cor_buf(n_pheno * n_pheno, 0.0);
    std::vector<double> std_buf(n_pheno * n_ind, 0.0);

    for (std::size_t pheno = 0; pheno < n_pheno; ++pheno) {
      if (ctx.phenotypes[pheno].comp_var(type) == 0) continue;
      double* buf_cur = &std_buf[pheno * n_ind];
      double mean_cur = ctx.phenotypes[pheno].comp_mean(type);
      double sd_cur = std::sqrt(ctx.phenotypes[pheno].comp_var(type));
      stats::standardise(
          n_ind, ctx.phenotypes[pheno](type), 1, buf_cur, 1, mean_cur, sd_cur);
    }

    for (std::size_t i = 0; i < n_pheno; ++i) {
      const double* ptr_i = &std_buf[i * n_ind];
      for (std::size_t j = i; j < n_pheno; ++j) {
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
  return [type](const SimulationContext& ctx) -> std::vector<double> {
    const std::size_t n_pheno = ctx.n_pheno;
    const std::size_t n_sex = ctx.n_sex;

    const std::vector<std::size_t> matching = ctx.model.state;

    const double* buf_male = ctx.phenotypes[0](type);
    const double* buf_female = buf_male + n_sex;
    const int lda_buf = 2 * n_sex;

    std::vector<double> cor_buf(n_pheno * n_pheno);

    std::vector<double> female_reordered(n_sex * n_pheno);
    for (std::size_t pheno = 0; pheno < n_pheno; ++pheno) {
      const double* src = buf_female + pheno * lda_buf;
      double* dst = female_reordered.data() + pheno * n_sex;
      for (std::size_t male_idx = 0; male_idx < n_sex; ++male_idx) {
        dst[male_idx] = src[matching[male_idx]];
      }
    }

    stats::cor(
        n_sex,
        n_pheno,
        n_pheno,
        buf_male,
        lda_buf,
        female_reordered.data(),
        n_sex,
        cor_buf.data(),
        n_pheno);

    return cor_buf;
  };
}

MetricFunc f_comp_mean(ComponentType type) {
  return [type](const SimulationContext& ctx) -> std::vector<double> {
    std::vector<double> comp_means(ctx.n_pheno);
    for (std::size_t pheno = 0; pheno < ctx.n_pheno; ++pheno)
      comp_means[pheno] = ctx.phenotypes[pheno].comp_mean(type);
    return comp_means;
  };
}

// phenotype component variance
MetricFunc f_comp_var(ComponentType type) {
  return [type](const SimulationContext& ctx) -> std::vector<double> {
    std::vector<double> comp_vars(ctx.n_pheno);
    for (std::size_t pheno = 0; pheno < ctx.n_pheno; ++pheno)
      comp_vars[pheno] = ctx.phenotypes[pheno].comp_var(type);
    return comp_vars;
  };
}

// latent phenotype heritability
std::vector<double> f_latent_h2(const SimulationContext& ctx) {
  const std::size_t n_pheno = ctx.n_pheno;
  const std::size_t n_sex = ctx.n_sex;
  std::vector<double> latent_h2(2 * n_pheno);

  for (std::size_t pheno = 0; pheno < n_pheno; ++pheno) {
    // gather the required buffer pointers
    const double* lat_gen_male = ctx.buf.latent(pheno, ComponentType::GENETIC);
    const double* lat_tot_male = ctx.buf.latent(pheno, ComponentType::TOTAL);
    const double* lat_gen_female = lat_gen_male + n_sex;
    const double* lat_tot_female = lat_tot_male + n_sex;

    // compute the heritabilities
    latent_h2[pheno] =
        stats::var(n_sex, lat_gen_male, 1) / stats::var(n_sex, lat_tot_male, 1);

    latent_h2[pheno + n_pheno] = stats::var(n_sex, lat_gen_female, 1) /
                                 stats::var(n_sex, lat_tot_female, 1);
  }

  return latent_h2;
}

// latent phenotype component correlation
MetricFunc f_latent_comp_cor(ComponentType type) {
  return [type](const SimulationContext& ctx) -> std::vector<double> {
    const std::size_t n_pheno = ctx.n_pheno;
    const std::size_t n_ind = ctx.n_ind;
    double* latent_comp = ctx.buf.latent(type);

    std::vector<double> cor_buf(n_pheno * n_pheno);

    stats::cor(n_ind, n_pheno, latent_comp, n_ind, cor_buf.data(), n_pheno);

    return cor_buf;
  };
}

// between-mate latent phenotype component correlation
MetricFunc f_latent_comp_xcor(ComponentType type) {
  return [type](const SimulationContext& ctx) -> std::vector<double> {
    const std::size_t n_pheno = ctx.n_pheno;
    const std::size_t n_sex = ctx.n_sex;
    const std::vector<std::size_t> matching = ctx.model.state;

    const double* latent_comp = ctx.buf.latent(type);
    const std::size_t lda = 2 * n_sex;

    // males are in the first n_sex rows, females in the next n_sex rows
    const double* latent_male = latent_comp;
    const double* latent_female_unordered = latent_comp + n_sex;

    // reorder females according to mating
    std::vector<double> latent_female(n_sex * n_pheno);
    for (std::size_t dim = 0; dim < n_pheno; ++dim) {
      const double* src = latent_female_unordered + dim * lda;
      double* dst = latent_female.data() + dim * n_sex;
      for (std::size_t male_idx = 0; male_idx < n_sex; ++male_idx) {
        dst[male_idx] = src[matching[male_idx]];
      }
    }

    std::vector<double> cor_buf(n_pheno * n_pheno);

    stats::cor(
        n_sex,
        n_pheno,
        n_pheno,
        latent_male,
        lda,
        latent_female.data(),
        n_sex,
        cor_buf.data(),
        n_pheno);

    return cor_buf;
  };
}

MetricFunc f_latent_comp_mean(ComponentType type) {
  return [type](const SimulationContext& ctx) -> std::vector<double> {
    const std::size_t n_pheno = ctx.n_pheno;
    const std::size_t n_sex = ctx.n_sex;
    const std::size_t n_ind = 2 * n_sex;
    const double* latent_all = ctx.buf.latent(type);

    std::vector<double> mean_buf(n_pheno);

    for (std::size_t dim = 0; dim < n_pheno; ++dim) {
      const double* latent_dim = latent_all + dim * n_ind;
      mean_buf[dim] = stats::mean(n_ind, latent_dim, 1);
    }

    return mean_buf;
  };
}

// latent phenotype component variance
MetricFunc f_latent_comp_var(ComponentType type) {
  return [type](const SimulationContext& ctx) -> std::vector<double> {
    const std::size_t n_pheno = ctx.n_pheno;
    const std::size_t n_sex = ctx.n_sex;
    const std::size_t n_ind = 2 * n_sex;
    const double* latent_all = ctx.buf.latent(type);

    std::vector<double> var_buf(n_pheno);

    for (std::size_t dim = 0; dim < n_pheno; ++dim) {
      const double* latent_dim = latent_all + dim * n_ind;
      var_buf[dim] = stats::var(n_ind, latent_dim, 1);
    }

    return var_buf;
  };
}

}  // namespace phenome
}  // namespace metrics
}  // namespace amsim
