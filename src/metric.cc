#include <cmath>
#include <filesystem>

#include <amsim/metric.h>
#include <amsim/component_type.h>

#if defined(__APPLE__)
  #include <Accelerate/Accelerate.h>
#else
  #include <cblas.h>
#endif

namespace amsim {
  Metric::Metric(MetricFunc f, const std::string &file, const std::string &name,
                 const std::size_t n_rows, const std::size_t n_cols,
                 std::optional<std::vector<std::string>> el_names)
    : f_(std::forward<MetricFunc>(f)),
      name_(name),
      n_rows_(n_rows),
      n_cols_(n_cols) {

    buf_.resize(n_rows_ * n_cols_);
    if (el_names) {
      named_ = true;
      el_names_ = std::move(*el_names);
    }
  }

  void Metric::set_dir(const std::string &dir) {
    if (!std::filesystem::exists(dir))
      std::filesystem::create_directory(dir); 

    file_path_ = (std::filesystem::path(dir) / file_path_).string();
    file_.open(file_path_); 

    if (!file_.is_open())
      throw std::runtime_error("couldn't open metric file");
  }

  void Metric::stream(PhenotypeList &phenotypes, Genome &genome) {
    buf_ = f_(phenotypes, genome);

    // print the header if this is the first iteration
    if (it_ == 0) {
      file_ << "it\t";
      for (std::size_t el = 0; el < buf_.size(); el++) {
        std::string el_name = named_ ? el_names_[el] : std::to_string(el);
        file_ << name_ << "::" << el_name;
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
    Metric make_metric(MetricFunc f, const std::string &name,
                       const std::string &file, const std::size_t n_rows,
                       const std::size_t n_cols,
                       std::optional<std::vector<std::string>> el_names) {
      return Metric{std::move(f), file, name, n_rows, n_cols, el_names};
    }

    namespace phenome {
      MetricFunc f_comp_cor(ComponentType type) {
        return [type](PhenotypeList &phenotypes, Genome &genome) -> std::vector<double> {
          const std::size_t n_pheno = phenotypes.size();
          const std::size_t n_ind   = phenotypes[0].get().n_ind();
          std::vector<double> cor_buf(n_pheno * n_pheno);
          std::vector<double> std_buf(n_pheno * n_ind);
          const std::vector<double> ones(n_ind, 1.0);

          for (std::size_t pheno = 0; pheno < n_pheno; pheno++) {
            if (phenotypes[pheno].get().comp_var(type) == 0)
              continue;
            double* buf_cur = &std_buf[pheno * n_ind];
            double mean_cur = phenotypes[pheno].get().comp_mean(type);
            double prec_cur = 1.0 / std::sqrt(phenotypes[pheno].get().comp_var(type));

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
              cor_buf[i * n_pheno + j] = (1.0 / static_cast<double>(n_ind)) *
                                          cblas_ddot(n_ind, ptr_i, 1, ptr_j, 1);
              cor_buf[j * n_pheno + i] = cor_buf[i * n_pheno + j];
            }
          }

          return cor_buf;
        };
      }

      // between-mate phenotype component correlation
      
      // phenotype heritabilities
      std::vector<double> f_h2(PhenotypeList &phenotypes, Genome &genome) {
        const std::size_t n_pheno = phenotypes.size();
        std::vector<double> pheno_h2(n_pheno);

        for (std::size_t pheno = 0; pheno < n_pheno; pheno++) {
          Phenotype &pheno_cur = phenotypes[pheno].get();
          pheno_h2[pheno] = pheno_cur.comp_var(ComponentType::GENETIC) /
                            pheno_cur.comp_var(ComponentType::TOTAL);
        }

        return pheno_h2;
      }

      // latent phenotype heritability

      // within-individual latent phenotype component correlation

      // between-mate latent phenotype component correlation
    }

    Metric comp_cor(const std::string& file, const std::size_t n_pheno,
                    ComponentType type) {
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

      return make_metric(phenome::f_comp_cor(type), file, name, n_pheno, n_pheno);
    }

    Metric pheno_h2(const std::string& file, const std::size_t n_pheno) {
      return make_metric(phenome::f_h2, file, "h2", n_pheno, 1);
    }
  }
}
