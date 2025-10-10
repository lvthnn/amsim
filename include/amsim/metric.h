#ifndef AMSIMCPP_METRIC_H
#define AMSIMCPP_METRIC_H

#pragma once
#include <string>
#include <fstream>

#include <amsim/genome.h>
#include <amsim/phenotype.h>

namespace amsim {

  template<class F, typename T>
  class Metric {
  public:
    Metric(F &&f, const std::string& name, const std::string& file,
           const std::size_t n_rows, const std::size_t n_cols = 1)
      : f_(std::forward<F>(f)),
        name_(name),
        n_rows_(n_rows),
        n_cols_(n_cols),
        file_(file) {
      buf_.resize(n_rows_ * n_cols_);
      n_it_ = 0;
    }

    // the vector of phenotypes should be managed by a simulation class, which
    // we've yet to implement
    void log(std::vector<Phenotype> &phenotypes, Genome &genome) {
      buf_ = f_(phenotypes, genome);

      // print the header first
      file_ << "it\t";
      for (std::size_t el = 0; el < buf_.size(); el++) {
        file_ << name_ << "[" << el << "]";
        if (el < buf_.size() - 1) file_ << "\t";
      }
      file_ << "\n";

      // stream data in buffer to file
      file_ << n_it_ << "\t";
      for (const T& val : buf_)
        file_ << "\t" << val << "\t";
      file_ << "\n";

      n_it_++;
    }

  private:
    F f_;
    const std::string name_;
    const std::size_t n_rows_;
    const std::size_t n_cols_;
    std::size_t n_it_;
    std::ofstream file_;
    std::vector<T> buf_;
  };

  namespace metrics {
    template<typename F>
    inline auto make_metric(F f, const std::string& name,
                            const std::string& file, std::size_t n_rows,
                            std::size_t n_cols = 1) {
      using T = typename std::invoke_result_t<F, std::vector<Phenotype>&, Genome&>::value_type;
      return Metric<F, T>(std::move(f), name, file, n_rows, n_cols);
    }

    namespace fgenome {
      std::vector<double> loc_maf(std::vector<Phenotype> &phenotype, Genome &genome) { return genome.v_lmaf(); }
      std::vector<double> loc_mean(std::vector<Phenotype> &phenotype, Genome &genome) { return genome.v_lmean(); }
      std::vector<double> loc_var(std::vector<Phenotype> &phenotype, Genome &genome) { return genome.v_lvar(); }

      std::vector<double> cor_gen(std::vector<Phenotype> &phenotype, Genome &genome) {

      }
    }

    inline auto loc_maf(const std::string& file, std::size_t n_loc) { return make_metric(fgenome::loc_maf, "loc_mafs", file, n_loc); }
    inline auto loc_mean(const std::string& file, std::size_t n_loc) { return make_metric(fgenome::loc_mean, "loc_means", file, n_loc); }
    inline auto loc_var(const std::string& file, std::size_t n_loc) { return make_metric(fgenome::loc_var, "loc_vars", file, n_loc); }
  }
}

#endif //AMSIMCPP_METRIC_H