#include <cmath>

#include <amsim/metric.h>
#include <amsim/componenttype.h>

#if defined(__APPLE__)
  #include <Accelerate/Accelerate.h>
#else
  #include <cblas.h>
#endif

namespace amsim {
  Metric::Metric(MetricFunc f, const std::string& file, std::string name,
                 const std::size_t n_rows, const std::size_t n_cols)
    : f_(std::forward<MetricFunc>(f)),
      file_(file, std::ios::out),
      name_(std::move(name)),
      n_rows_(n_rows),
      n_cols_(n_cols) {
    if (!file_.is_open()) throw std::runtime_error("couldn't open metric file");
    buf_.resize(n_rows_ * n_cols_);
    it_ = 0;
  }

  void Metric::stream(std::vector<Phenotype>& phenotypes, Genome& genome) {
    buf_ = f_(phenotypes, genome);

    // print the header if this is the first iteration
    if (it_ == 0) {
      file_ << "it\t";
      for (std::size_t el = 0; el < buf_.size(); el++) {
        file_ << name_ << "::" << el;
        if (el < buf_.size() - 1) file_ << "\t";
      }
      file_ << "\n";
    }

    // stream data in buffer to file
    file_ << it_ << "\t";
    for (std::size_t el = 0; el < buf_.size(); el++) {
      file_ << buf_[el];
      if (el < buf_.size() - 1) file_ << "\t";
    }
    file_ << "\n";

    it_++;
  }

  namespace metrics {
    Metric make_metric(MetricFunc f, const std::string& name,
                       const std::string& file, const std::size_t n_rows,
                       const std::size_t n_cols) {
      return Metric{std::move(f), name, file, n_rows, n_cols};
    }

    namespace phenome {
      MetricFunc fcomp_cor(ComponentType type) {
        return [type](std::vector<Phenotype>& phenotypes, Genome& genome) -> std::vector<double> {
          const std::size_t n_pheno = phenotypes.size();
          const std::size_t n_ind   = phenotypes[0].n_ind();
          std::vector<double> cor_buf(n_pheno * n_pheno);
          std::vector<double> std_buf(n_pheno * n_ind);
          const std::vector<double> ones(n_ind, 1.0);

          for (std::size_t pheno = 0; pheno < n_pheno; pheno++) {
            if (phenotypes[pheno].comp_var(type) == 0)
              continue;
            double* buf_cur = &std_buf[pheno * n_ind];
            double mean_cur = phenotypes[pheno].comp_mean(type);
            double prec_cur = 1.0 / std::sqrt(phenotypes[pheno].comp_var(type));

            cblas_dcopy(n_ind, phenotypes[pheno](type), 1, buf_cur, 1);
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
              cor_buf[i * n_pheno + j] = (1.0 / static_cast<double>(n_ind)) * cblas_ddot(n_ind, ptr_i, 1, ptr_j, 1);
              cor_buf[j * n_pheno + i] = cor_buf[i * n_pheno + j];
            }
          }

          return cor_buf;
        };
      }
    }

    Metric comp_cor(const std::string& file, const std::size_t n_pheno, ComponentType type) {
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

      return make_metric(phenome::fcomp_cor(type), name, file, n_pheno, n_pheno);
    }
  }
}