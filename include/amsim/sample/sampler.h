#pragma once

#include <amsim/core.h>
#include <amsim/estimate.h>
#include <amsim/sample.h>

#include <filesystem>
#include <variant>

namespace amsim {

namespace details {

template <Proband P>
auto get_prefix = [](auto member_enum) -> std::string_view {
  for (std::size_t i = 0; i < ProbandData<P>::ProbandSize; ++i) {
    if (ProbandData<P>::ProbandMembers[i].member == member_enum)
      return ProbandData<P>::ProbandMembers[i].prefix;
  }
  return "";
};

}  // namespace details

inline WeightFunction Uniform() {
  return [](const Eigen::MatrixXd& /*agg*/, Eigen::VectorXd& res) {
    res.setConstant(1.0);
  };
}

inline WeightFunction Logistic(const Eigen::VectorXd& effects) {
  return [effects](const Eigen::MatrixXd& agg, Eigen::VectorXd& res) {
    res = 1.0 / (1.0 + (-agg * effects).array().exp());
  };
}

struct SampleDescription {
  std::string name;
  std::string proband_type;
  std::size_t n_probands;
  std::vector<std::string> on;
  std::vector<std::string> of;
  std::string agg = "identity";
  std::string weight_function = "uniform()";
  std::vector<std::string> estimators;
};

// user-facing specification struct
template <Proband P>
struct Sample {
  using ProbandEnum = ProbandData<P>::ProbandEnum;

  // selection parameters — which data, from whom, and what to do
  std::string name;
  std::size_t n_probands;
  std::optional<std::vector<std::string>> on;
  std::optional<std::vector<Component>> on_components;
  ProbandEnum of = ProbandData<P>::ProbandDefault;
  Aggregator agg = ProbandData<P>::AggDefault;
  bool decompress_genotypes = false;

  // weighting method — probability of selecting based on proband aggregate
  WeightFunction weighting;

  // estimators are declared here
  SampleEstimators<P> estimators;
};

using SampleSpec = std::variant<
    Sample<Proband::Individual>,
    Sample<Proband::Mate>,
    Sample<Proband::Family>>;

class Sampler {
 public:
  template <Proband P>
  explicit Sampler(Sample<P> sample, const Params& params)
      : self_(std::make_unique<Model<P>>(std::move(sample), params)){};

  void operator()(const State& state) { (*self_)(state); }

 private:
  // type-erased computing interface
  struct Concept {
    virtual ~Concept() = default;
    virtual void draw(const State& state) = 0;
    virtual void estimate(const State& state) = 0;
    virtual void operator()(const State& state) = 0;
  };

  // implementation strategy
  template <Proband P>
  struct Model : Concept {
    explicit Model(Sample<P> sample, const Params& params)
        : name(sample.name),
          n_probands(sample.n_probands),
          sample_dir(params.sim.out_dir / name),
          of(sample.of),
          agg(std::move(sample.agg)),
          weighting(std::move(sample.weighting)),
          decompress_genotypes(sample.decompress_genotypes),
          n_probands_total(
              P == Proband::Individual ? params.geno.n_ind
                                       : params.geno.n_ind / 2),
          n_sex(params.geno.n_ind / 2),
          n_members(__builtin_popcountll(static_cast<uint8_t>(sample.of))),
          n_pheno(params.pheno.n_pheno),
          n_loc(params.geno.n_loc),
          n_on(sample.on.has_value() ? sample.on.value().size() : 1),
          on_components(
              sample.on_components.has_value()
                  ? std::move(sample.on_components.value())
                  : std::vector<Component>(n_on, Component::Total)),
          members(n_probands_total * n_members, n_on),
          aggregates(n_probands_total, n_on),
          unif(n_probands_total),
          keys(n_probands_total),
          names(params.pheno.names),
          phenotypes(sample.n_probands * ProbandData<P>::ProbandSize, n_pheno),
          selected(n_probands_total) {
      if (n_probands > n_probands_total)
        throw std::runtime_error(
            "n_probands (" + std::to_string(n_probands) + ") exceeds " +
            "available sampling units (" + std::to_string(n_probands_total) +
            ") for sample '" + name + "'");

      if (sample.on.has_value())
        for (const auto& pheno_name : sample.on.value())
          on_indices.push_back(params.pheno.pheno_ids.at(pheno_name));

      std::filesystem::create_directories(sample_dir);
      for (const auto& estimator : sample.estimators)
        estimators.emplace_back(estimator(params, n_probands, sample_dir));
    };

    // fields obtained from sampler specification
    std::string name;
    std::size_t n_probands;
    std::filesystem::path sample_dir;
    ProbandData<P>::ProbandEnum of;
    Aggregator agg;
    WeightFunction weighting;
    std::vector<std::unique_ptr<SampleEstimatorStrategy<P>>> estimators;
    bool decompress_genotypes;

    // derived dimensions
    std::size_t n_probands_total;
    std::size_t n_sex;
    std::size_t n_members;
    std::size_t n_pheno;
    std::size_t n_loc;
    std::size_t n_on;

    // phenotype indices and components for sampling
    std::vector<std::size_t> on_indices;
    std::vector<Component> on_components;

    // buffers to store sampling keys and data
    Eigen::MatrixXd members;
    Eigen::MatrixXd aggregates;
    Eigen::VectorXd unif;
    Eigen::VectorXd keys;

    // extracted proband phenotypes
    std::vector<std::string> names;
    Eigen::MatrixXd phenotypes;
    Eigen::MatrixXd genotypes;
    std::vector<std::size_t> selected;

    void fillAggregates(const State& state);
    void extractProbands(const State& state);

    void addGenotype(
        const State& state,
        std::size_t proband_id,
        std::size_t locus,
        const ProbandMember<typename ProbandData<P>::ProbandEnum>& member,
        std::uint8_t& byte,
        std::size_t& bit_pos,
        std::fstream& bed) const;

    void writeBIM() const;
    void writeBED(const State& state) const;
    void writeFAM(const State& state) const;
    void writePHENO(
        const State& state, Component type = Component::Total) const;

    void writePLINK(const State& state) const {
      writeBIM();
      writeBED(state);
      writeFAM(state);
      for (auto comp :
           {Component::Total,
            Component::Genetic,
            Component::Environmental,
            Component::Vertical})
        writePHENO(state, comp);
    }

    void draw(const State& state) override;
    void estimate(const State& state) override;
    void operator()(const State& state) override;
  };

  std::unique_ptr<Concept> self_;
};

template <Proband P>
inline void Sampler::Model<P>::fillAggregates(const State& state) {
  if (on_indices.empty()) return;

  for (std::size_t on = 0; on < n_on; ++on) {
    std::size_t pheno_id = on_indices[on];
    Component comp = on_components[on];

    for (std::size_t prob = 0; prob < n_probands_total; ++prob) {
      std::size_t row = prob * n_members;
      std::size_t member_idx = 0;

      for (std::size_t m = 0; m < ProbandData<P>::ProbandSize; ++m) {
        const auto& member = ProbandData<P>::ProbandMembers[m];
        if (of & member.self)
          members(row + member_idx++, on) =
              member_pheno(state, member, prob, pheno_id, comp);
      }
    }
  }

  // aggregate each of the probands
  for (std::size_t prob = 0; prob < n_probands_total; ++prob) {
    auto proband = members.middleRows(prob * n_members, n_members);

    if (agg == Aggregator::Mean)
      aggregates.row(prob) = proband.colwise().mean();

    else if (agg == Aggregator::Max)
      aggregates.row(prob) = proband.colwise().maxCoeff();

    else if (agg == Aggregator::Min)
      aggregates.row(prob) = proband.colwise().minCoeff();

    else if (agg == Aggregator::Identity)
      aggregates.row(prob) = proband.row(0);
  }
}

template <Proband P>
inline void Sampler::Model<P>::extractProbands(const State& state) {
  for (std::size_t pheno_id = 0; pheno_id < n_pheno; ++pheno_id) {
    for (std::size_t prob = 0; prob < n_probands; ++prob) {
      std::size_t prob_id = selected[prob];
      std::size_t row = prob * ProbandData<P>::ProbandSize;

      for (std::size_t m = 0; m < ProbandData<P>::ProbandSize; ++m) {
        const auto& member_data = ProbandData<P>::ProbandMembers[m];
        phenotypes(row + m, pheno_id) =
            member_pheno(state, member_data, prob_id, pheno_id);
      }
    }
  }
}

template <Proband P>
inline void Sampler::Model<P>::addGenotype(
    const State& state,
    std::size_t proband_id,
    std::size_t locus,
    const ProbandMember<typename ProbandData<P>::ProbandEnum>& member,
    std::uint8_t& byte,
    std::size_t& bit_pos,
    std::fstream& bed) const {
  byte |= (member_geno_plink(state, member, proband_id, locus) << bit_pos);
  bit_pos += 2;
  if (bit_pos == 8) {
    bed.put(byte);
    byte = 0;
    bit_pos = 0;
  }
}

template <Proband P>
inline void Sampler::Model<P>::writeBIM() const {
  std::fstream bim_file(sample_dir / "data.bim", std::ios::out);

  if (!bim_file.is_open())
    throw std::runtime_error(
        "Could not open BIM file output stream " +
        (sample_dir / "data.bim").string());

  for (std::size_t loc = 0; loc < n_loc; ++loc)
    bim_file << std::format("1\tSNP{}\t0\t{}\tA\tG\n", loc, loc + 1);
}

template <Proband P>
inline void Sampler::Model<P>::writeBED(const State& state) const {
  if (state.geno().view() != HaploView::LocusMajor)
    throw std::runtime_error(
        "Sampler::Model<P>::writeBED: require locus-major layout");

  // open the file
  std::fstream bed_file(
      sample_dir / "data.bed", std::ios::out | std::ios::binary);

  if (!bed_file.is_open())
    throw std::runtime_error(
        "Could not open BED file stream " + (sample_dir / "data.bed").string());

  // magic numbers
  bed_file.put(0x6c).put(0x1b).put(0x01);

  // start looping
  for (std::size_t loc = 0; loc < n_loc; ++loc) {
    std::uint8_t byte = 0;
    std::size_t bit_pos = 0;

    for (std::size_t prob = 0; prob < n_probands; ++prob) {
      std::size_t id = selected[prob];

      // add the genotypes of all members in the proband
      for (std::size_t m = 0; m < ProbandData<P>::ProbandSize; ++m)
        addGenotype(
            state,
            id,
            loc,
            ProbandData<P>::ProbandMembers[m],
            byte,
            bit_pos,
            bed_file);
    }

    // flush partial buffer at end of locus
    if (bit_pos != 0) bed_file.put(byte);
  }
}

template <Proband P>
inline void Sampler::Model<P>::writeFAM(const State& state) const {
  std::fstream fam_file(sample_dir / "data.fam", std::ios::out);

  if (!fam_file.is_open())
    throw std::runtime_error(
        "Could not open FAM file output stream " +
        (sample_dir / "data.bim").string());

  for (std::size_t prob = 0; prob < n_probands; ++prob) {
    std::size_t id = selected[prob];

    for (std::size_t m = 0; m < ProbandData<P>::ProbandSize; ++m) {
      auto member = ProbandData<P>::ProbandMembers[m];

      std::size_t fid = prob;
      std::size_t iid = member_index(state, member, id);

      int sex = (member.sex == Sex::Unknown)
                    ? (static_cast<int>(iid >= n_sex) + 1)
                    : static_cast<int>(member.sex);

      std::string istr = member_id<P>(state, member.self, id);
      std::string pstr = member_id<P>(state, member.father, id);
      std::string mstr = member_id<P>(state, member.mother, id);

      fam_file << std::format(
          "FAM{}\t{}\t{}\t{}\t{}\t-9\n", fid, istr, pstr, mstr, sex);
    }
  }
}

template <Proband P>
inline void Sampler::Model<P>::writePHENO(
    const State& state, Component type) const {
  auto path = (type == Component::Total)
                  ? sample_dir / "data.pheno"
                  : sample_dir / std::format("data.{}.pheno", to_string(type));

  std::fstream pheno_file(path, std::ios::out);

  if (!pheno_file.is_open())
    throw std::runtime_error(
        "Could not open PHENO file output stream " + path.string());

  std::string header = "FID\tIID";
  for (const std::string& pheno : names) header += "\t" + pheno;
  header += "\n";

  pheno_file << header;

  for (std::size_t prob = 0; prob < n_probands; ++prob) {
    std::size_t id = selected[prob];

    for (std::size_t m = 0; m < ProbandData<P>::ProbandSize; ++m) {
      auto member = ProbandData<P>::ProbandMembers[m];

      std::size_t fid = prob;
      std::string istr = member_id<P>(state, member.self, id);

      pheno_file << std::format("FAM{}\t{}", fid, istr);

      for (std::size_t pheno = 0; pheno < n_pheno; ++pheno)
        pheno_file << std::format(
            "\t{}", member_pheno(state, member, id, pheno, type));

      pheno_file << "\n";
    }
  }
}

template <Proband P>
inline void Sampler::Model<P>::draw(const State& state) {
  // fill the aggregation buffer with phenotypes of the respective members
  fillAggregates(state);

  // generate uniform weights for Efraimidis-Spirakis
  rng::UniformRange::fill(unif.data(), n_probands_total);

  // run the aggregation buffer
  weighting(aggregates, keys);
  keys = unif.array().pow(1.0 / keys.array());

  std::iota(selected.begin(), selected.end(), 0);
  std::nth_element(
      selected.begin(),
      selected.begin() + n_probands,
      selected.end(),
      [&](auto a, auto b) { return keys(a) > keys(b); });

  // extract selected proband data
  extractProbands(state);

  // write to PLINK if necessary
  writePLINK(state);
}

template <Proband P>
inline void Sampler::Model<P>::estimate(const State& state) {
  for (const auto& estimator : estimators)
    (*estimator)(state.gen, state.rep, phenotypes, genotypes);
}

template <Proband P>
inline void Sampler::Model<P>::operator()(const State& state) {
  std::vector<std::filesystem::path> to_remove;
  for (const auto& entry : std::filesystem::directory_iterator(sample_dir))
    to_remove.push_back(entry.path());
  for (const auto& p : to_remove) std::filesystem::remove(p);
  draw(state);
  estimate(state);
}

template struct Sampler::Model<Proband::Individual>;
template struct Sampler::Model<Proband::Mate>;
template struct Sampler::Model<Proband::Family>;

class ComputeSampleEstimates {
 public:
  ComputeSampleEstimates(
      const Params& params, const std::vector<SampleSpec>& samples) {
    for (const auto& sample : samples)
      estimators_.emplace_back(std::visit(
          [&params](auto&& spec) { return Sampler(spec, params); }, sample));
  }

  void operator()(const State& state);

 private:
  std::vector<Sampler> estimators_;
};

inline void ComputeSampleEstimates::operator()(const State& state) {
  for (auto& estimator : estimators_) estimator(state);
}

}  // namespace amsim
