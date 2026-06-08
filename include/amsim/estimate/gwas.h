#pragma once

#include <amsim/core/params.h>
#include <amsim/core/utils.h>
#include <amsim/estimate/sample.h>
#include <amsim/io/writer.h>
#include <amsim/sample/proband.h>

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <utility>

namespace amsim {

template <Proband P>
class GWASEstimator : public SampleEstimatorStrategy<P> {
 public:
  GWASEstimator(
      const Params& params,
      std::string sample_name,
      std::filesystem::path sample_dir,
      std::string name = "gwas",
      std::size_t n_pcs = 0,
      double pval_threshold = 5e-8)
      : SampleEstimatorStrategy<P>(
            std::move(name),
            std::move(sample_name),
            sample_dir,
            params.pheno.names,
            {"l2_effect", "fpr", "tpr", "pgs_r2", "pgs_rmse"},
            params.pheno.n_pheno,
            5),
        n_pheno_(params.pheno.n_pheno),
        n_loc_(params.geno.n_loc),
        pval_threshold_(pval_threshold),
        n_pcs_(n_pcs),
        pheno_names_(params.pheno.names),
        beta_true_(n_pheno_, Eigen::VectorXd::Zero(n_loc_)),
        causal_mask_(n_pheno_, Eigen::VectorXd::Zero(n_loc_)) {
    utils::check_plink2();
    for (std::size_t p = 0; p < n_pheno_; ++p)
      for (std::size_t i = 0; i < params.pheno.pheno_loc[p].size(); ++i) {
        std::size_t loc = params.pheno.pheno_loc[p][i];
        beta_true_[p](loc) = params.pheno.pheno_effects[p](i);
        causal_mask_[p](loc) = 1.0;
      }
  }

  void compute() override {
    runGWAS();

    for (std::size_t p = 0; p < n_pheno_; ++p) {
      auto [betas, pvals] = parseGWAS(p);

      auto sig = (pvals.array() < pval_threshold_).template cast<double>();
      auto& mask = causal_mask_[p];
      double n_causal = mask.sum();
      double n_null = n_loc_ - n_causal;

      double l2 = (betas - beta_true_[p]).norm();
      double tpr = n_causal > 0 ? (sig * mask.array()).sum() / n_causal : 0.0;
      double fpr =
          n_null > 0 ? (sig * (1.0 - mask.array())).sum() / n_null : 0.0;
      auto [r2, rmse] = pgsMetrics(p);

      this->data_.row(p) << l2, fpr, tpr, r2, rmse;
    }
  }

 private:
  std::size_t n_pheno_;
  std::size_t n_loc_;
  double pval_threshold_;
  std::size_t n_pcs_;
  std::vector<std::string> pheno_names_;

  std::vector<Eigen::VectorXd> beta_true_;
  std::vector<Eigen::VectorXd> causal_mask_;

  void runGWAS() {
    if (n_pcs_ > 0) {
      utils::system_throttled(std::format(
          "plink2 --bfile {} --pca approx {} --out {}",
          (this->sample_dir_ / "data").string(),
          n_pcs_,
          (this->sample_dir_ / (this->name_ + "_pca")).string()));

      utils::system_throttled(std::format(
          "plink2 --bfile {} --pheno {} --covar {} "
          "--glm --variance-standardize --no-psam-pheno "
          "--threads 1 --out {}",
          (this->sample_dir_ / "data").string(),
          (this->sample_dir_ / "data.pheno").string(),
          (this->sample_dir_ / (this->name_ + "_pca.eigenvec")).string(),
          (this->sample_dir_ / this->name_).string()));
    } else {
      utils::system_throttled(std::format(
          "plink2 --bfile {} --pheno {} --glm allow-no-covars "
          "--variance-standardize --no-psam-pheno "
          "--threads 1 --out {}",
          (this->sample_dir_ / "data").string(),
          (this->sample_dir_ / "data.pheno").string(),
          (this->sample_dir_ / this->name_).string()));
    }
  }

  std::pair<Eigen::VectorXd, Eigen::VectorXd> parseGWAS(std::size_t p) {
    Eigen::VectorXd betas = Eigen::VectorXd::Zero(n_loc_);
    Eigen::VectorXd pvals = Eigen::VectorXd::Ones(n_loc_);

    std::ifstream f(
        this->sample_dir_ / std::format("{}.{}.glm.linear", this->name_, pheno_names_[p]));

    if (!f) return {betas, pvals};

    std::string line;
    std::getline(f, line);  // header
    while (std::getline(f, line)) {
      std::istringstream ss(line);
      std::string chrom;
      std::string pos;
      std::string id;
      std::string ref;
      std::string alt;
      std::string prov_ref;
      std::string a1;
      std::string omitted;
      double a1_freq;
      std::string test;
      std::size_t obs;
      double beta;
      double se;
      double t;
      double pval;

      ss >> chrom >> pos >> id >> ref >> alt >> prov_ref >> a1 >> omitted >>
          a1_freq >> test >> obs >> beta >> se >> t >> pval;

      if (test != "ADD") continue;

      std::size_t loc = std::stoul(id.substr(3));
      betas(loc) = beta;
      pvals(loc) = pval;
    }
    return {betas, pvals};
  }

  std::pair<double, double> pgsMetrics(std::size_t p) {
    constexpr auto NaN = std::numeric_limits<double>::quiet_NaN();
    auto glm = this->sample_dir_ /
               std::format("{}.{}.glm.linear", this->name_, pheno_names_[p]);
    auto score_file =
        this->sample_dir_ / std::format("{}_score_{}.txt", this->name_, p);
    auto pfix = this->sample_dir_ / std::format("{}_pgs_{}", this->name_, p);

    {
      std::ifstream glm_in(glm);
      std::ofstream score_out(score_file);
      if (!glm_in) return {NaN, NaN};

      std::string line;
      std::getline(glm_in, line);
      score_out << line << "\n";

      bool any_hits = false;
      while (std::getline(glm_in, line)) {
        std::istringstream ss(line);
        std::string chrom;
        std::string pos;
        std::string id;
        std::string ref;
        std::string alt;
        std::string prov_ref;
        std::string a1;
        std::string omitted;
        std::string test;
        double a1_freq;
        double beta;
        double se;
        double t;
        double pval;
        std::size_t obs;
        ss >> chrom >> pos >> id >> ref >> alt >> prov_ref >> a1 >> omitted >>
            a1_freq >> test >> obs >> beta >> se >> t >> pval;
        if (test == "ADD" && pval < pval_threshold_) {
          score_out << line << "\n";
          any_hits = true;
        }
      }
      if (!any_hits) return {0.0, NaN};
    }

    utils::system_throttled(std::format(
        "plink2 --bfile {} --score {} 3 7 12 header "
        "--threads 1 --out {}",
        (this->sample_dir_ / "data").string(),
        score_file.string(),
        pfix.string()));

    std::ifstream pgs_file(pfix.string() + ".sscore");
    std::ifstream gen_file(this->sample_dir_ / "data.genetic.pheno");
    if (!pgs_file || !gen_file) return {NaN, NaN};

    std::string pgs_line;
    std::string gen_line;
    std::getline(pgs_file, pgs_line);  // header
    std::getline(gen_file, gen_line);  // header

    std::vector<double> pgs_vals;
    std::vector<double> gen_vals;
    while (std::getline(pgs_file, pgs_line) &&
           std::getline(gen_file, gen_line)) {
      std::istringstream pgs_ss(pgs_line);
      std::string fid;
      std::string iid;
      double allele_ct;
      double dosage_sum;
      double score;
      pgs_ss >> fid >> iid >> allele_ct >> dosage_sum >> score;
      pgs_vals.push_back(score);

      std::istringstream gen_ss(gen_line);
      std::string gen_fid;
      std::string gen_iid;
      gen_ss >> gen_fid >> gen_iid;
      double gen_val = 0.0;
      for (std::size_t col = 0; col <= p; ++col) gen_ss >> gen_val;
      gen_vals.push_back(gen_val);
    }

    if (pgs_vals.empty() || pgs_vals.size() != gen_vals.size())
      return {NaN, NaN};

    Eigen::VectorXd pgs =
        Eigen::Map<Eigen::VectorXd>(pgs_vals.data(), pgs_vals.size());
    Eigen::VectorXd gen =
        Eigen::Map<Eigen::VectorXd>(gen_vals.data(), gen_vals.size());

    auto pc = pgs.array() - pgs.mean();
    auto gc = gen.array() - gen.mean();
    double num = (pc * gc).sum();
    double r2 = (num * num) / (pc.square().sum() * gc.square().sum());
    double rmse = std::sqrt((pgs - gen).squaredNorm() / pgs.size());
    return {r2, rmse};
  }
};

template <Proband P>
inline SampleEstimator<P> SampleGWASEstimator(
    std::string name = "gwas",
    std::size_t n_pcs = 0,
    double pval_threshold = 5e-8) {
  return SampleEstimator<P>{
      .name = name,
      .fn = [name = std::move(name), n_pcs, pval_threshold](
                const Params& params,
                std::size_t /*n_probands*/,
                const std::filesystem::path& sample_dir) {
        return std::make_unique<GWASEstimator<P>>(
            params,
            sample_dir.filename().string(),
            sample_dir,
            name,
            n_pcs,
            pval_threshold);
      }};
}

}  // namespace amsim
