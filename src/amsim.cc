// This file is part of amsim, copyright (C) 2025-2026 Kári Hlynsson.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <amsim/core.h>
#include <amsim/estimate.h>
#include <amsim/init.h>
#include <amsim/io/config.h>
#include <amsim/io/parse.h>
#include <amsim/simulation.h>
#include <getopt.h>

#include <Eigen/Dense>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

constexpr std::string OS_NAME() {
  std::string os = AMSIM_OS;
  if (os == "Darwin") return "macOS";
  return os;
}

void display_version() {
  std::string version_str = std::format(
      "amsim v{} {} {} ({} build, {})",
      AMSIM_VERSION,
      OS_NAME(),
      AMSIM_ARCH,
      AMSIM_BUILD_TYPE,
      __DATE__);

  std::cout << version_str << std::endl;
}

void display_header() {
  std::string header_str = std::format(
      "{:<40}{:>40}\n{:<40}{:>40}",
      std::format(
          "amsim v{} {} {} ({} build)",
          AMSIM_VERSION,
          OS_NAME(),
          AMSIM_ARCH,
          AMSIM_BUILD_TYPE),
      "https://github.com/lvthnn/amsim",
      "(C) 2025-2026 Kári Hlynsson",
      "GNU General Public License v3");

  std::cout << header_str << std::endl;
}

void display_help_short() {
  std::string help_str = R"(
  amsim <option(s)>
  amsim --config <config_file> <option(s)>
  amsim --template <config_path>
  amsim --template-cmd

To see all options, run "amsim --help".)";

  std::cout << help_str << std::endl;
}

void display_help() {
  std::string help_str = R"(
For full documentation and usage examples, see https://github.com/lvthnn/amsim.

program options:
  -h, --help
  -v, --version

global options:
  --n-individuals <n_ind>
  --n-loci <n_loc>
  --n-generations <n_gen>
  --n-threads <n_thr>
  --n-replicates <n_rep>
  --random-seed <rng_seed>
  --output-dir <out_dir>
  --output-name <out_name>
  --log-level {debug | info | warning | error | none}
  --log-to-output
  --save-config
  --share-init-state
  --no-run

genome options:
  --locus-maf [--val <value> | --file <path> | --dist <distribution>]
  --locus-rec [--val <value> | --file <path> | --dist <distribution>]
  --locus-mut [--val <value> | --file <path> | --dist <distribution>]

phenome options:
  --phenotype <pheno_name>
    --var-genetic <var_gen>
    --var-environmental <var_env>
    --var-vertical <var_vert>
    --n-loci <n_loc>
    --loci <loc>
    --effects <eff>
  --pheno-gen-cor [--value <vals> | --file <path> | --singular-values <svals>]
  --pheno-env-cor [--value <vals> | --file <path> | --singular-values <svals>]

mating options:
  --mating {assortative | random}
    --mate-cor [--value <vals> | --file <path> | --singular-values <svals>]
    --tol-inf <tol-inf>
    --max-itr <max-itr>
    --temp-init <tmp_init>
    --temp-decay <tmp_decay>

sampling and estimation options:
  --estimator {genotype-mean |
               genotype-var |
               genotype-maf |
               genotype-cov |
               genotype-cor |
               heritability |
               pheno-mean({genetic | environmental | vertical | total}) |
               pheno-var({genetic | environmental | vertical | total}) |
               pheno-cov({genetic | environmental | vertical | total}) |
               mate-cor({genetic | environmental | vertical | total})}
  --sample-estimator <name>
    --type {gwas(<n_pcs>, <pval_thresh>) | haseman-elston | greml | external}
    --exec <exec_cmd>
    --n-rows <n_rows>
    --n-cols <n_cols>
    --row-names <row_names>
    --col-names <col_names>
  --sample <sample_name>
    --proband {individual | mate | family}
    --n_probands <n-prob>
    --weight {logistic(<pheno_coefs>) | uniform()}
    --on <weight_on>
    --of <weight_of>
    --agg {max | min | mean | identity}
    --estimators <est_names>
  )";

  std::cout << help_str << std::endl;
}

enum Context : uint32_t {
  // Base domain
  Global = 0,
  Genome = 1 << 0,
  Phenotype = 1 << 1,
  Mating = 1 << 2,
  Sample = 1 << 3,
  SampleEstimator = 1 << 4,
  PopulationEstimator = 1 << 5,

  // Subcomponent
  InitMAFs = 1 << 6,
  RecombinationProbs = 1 << 7,
  MutationProbs = 1 << 8,
  CausalLoci = 1 << 9,
  GeneticComponent = 1 << 10,
  EnvironmentalComponent = 1 << 11,
  EffectAssignment = 1 << 12,
  CorrelationMatrix = 1 << 13,

  // Derived contexts
  GenomeInitMAFs = Genome | InitMAFs,
  GenomeRecombinationProbs = Genome | RecombinationProbs,
  GenomeMutationProbs = Genome | MutationProbs,
  GenomeProbabilities = InitMAFs | RecombinationProbs | MutationProbs,
  PhenotypeCausalLoci = Phenotype | CausalLoci,
  PhenotypeGeneticCorrelation = Phenotype | GeneticComponent,
  PhenotypeEnvironmentalCorrelation = Phenotype | EnvironmentalComponent,
  MatingCorrelation = Mating | CorrelationMatrix,

  // Base domain
  ContextDomain = Global | Genome | Phenotype | Mating | Sample |
                  SampleEstimator | PopulationEstimator,

  // To allow help messages for flags
  Help = 1 << 14
};

Context context_domain(Context context) {
  return static_cast<Context>(context & ContextDomain);
}

enum Option {
  GlobalNumIndividuals,
  GlobalNumLoci,
  GlobalNumGenerations,
  GlobalNumThreads,
  GlobalNumReplicates,
  GlobalRandomSeed,
  GlobalOutputDirectory,
  GlobalOutputName,
  GlobalLogLevel,
  GlobalLogNoFile,
  GlobalSaveConfig,
  GlobalShareInitState,
  GlobalNoRun,
  GenomeLocusInitMAFs,
  GenomeLocusRecombinationProbs,
  GenomeLocusMutationProbs,
  PhenotypeDecl,
  PhenotypeVarGenetic,
  PhenotypeVarEnvironmental,
  PhenotypeVarVertical,
  PhenotypeLocusEffects,
  PhenotypeLocusIndices,
  PhenotypeGeneticCor,
  PhenotypeEnvironmentalCor,
  MatingDecl,
  MatingCor,
  MatingErrorTolerance,
  MatingMaxIterations,
  MatingAnnealingTempInit,
  MatingAnnealingTempDecay,
  SampleDecl,
  SampleProbandType,
  SampleNumProbands,
  SampleWeightFunction,
  SampleWeightOnPhenotypes,
  SampleWeightOfMembers,
  SampleWeightAggregation,
  SampleEstimators,
  SampleEstimatorDecl,
  SampleEstimatorType,
  SampleEstimatorExec,
  SampleEstimatorNumRows,
  SampleEstimatorNumCols,
  SampleEstimatorRowNames,
  SampleEstimatorColNames,
  PopulationEstimatorDecl,
  VirtualNumLoci,  // these options are context-sensitive
  VirtualFile,
  VirtualValue,
  VirtualDistribution,
  VirtualSingularValues
};

std::unordered_map<
    std::string,
    std::function<amsim::PopulationEstimator(const std::vector<std::string>&)>>
    estimators = {
        {"genotype-mean",
         [](const std::vector<std::string>& /*s*/) {
           return amsim::PopulationGenotypeMean();
         }},
        {"genotype-var",
         [](const std::vector<std::string>& /*s*/) {
           return amsim::PopulationGenotypeVar();
         }},
        {"genotype-maf",
         [](const std::vector<std::string>& /*s*/) {
           return amsim::PopulationGenotypeMAF();
         }},
        {"genotype-cov",
         [](const std::vector<std::string>& /*s*/) {
           return amsim::PopulationGenotypeCov();
         }},
        {"genotype-cor",
         [](const std::vector<std::string>& /*s*/) {
           return amsim::PopulationGenotypeCor();
         }},
        {"heritability",
         [](const std::vector<std::string>& /*s*/) {
           return amsim::PopulationHeritability();
         }},
        {"pheno-mean",
         [](const std::vector<std::string>& s) {
           amsim::Component component = amsim::Component::Total;
           if (!s.empty()) component = amsim::Component_from_string(s[0]);
           return PopulationComponentMean(component);
         }},
        {"pheno-var",
         [](const std::vector<std::string>& s) {
           amsim::Component component = amsim::Component::Total;
           if (!s.empty()) component = amsim::Component_from_string(s[0]);
           return PopulationComponentMean(component);
         }},
        {"pheno-cor",
         [](const std::vector<std::string>& s) {
           amsim::Component component = amsim::Component::Total;
           if (!s.empty()) component = amsim::Component_from_string(s[0]);
           return PopulationComponentCor(component);
         }},
        {"pheno-cov",
         [](const std::vector<std::string>& s) {
           amsim::Component component = amsim::Component::Total;
           if (!s.empty()) component = amsim::Component_from_string(s[0]);
           return PopulationComponentCov(component);
         }},
        {"mate-cor", [](const std::vector<std::string>& s) {
           amsim::Component component = amsim::Component::Total;
           if (!s.empty()) component = amsim::Component_from_string(s[0]);
           return PopulationMateCor(component);
         }}};

int main(int argc, char* argv[]) {
  try {
    if (argc == 1) {
      display_header();
      display_help_short();
      exit(EXIT_SUCCESS);
    }

    Context context = Context::Global;
    amsim::Simulation simulation;

    std::size_t n_replicates = 1;
    std::size_t n_threads = 1;

    bool run = true;
    bool save_config = false;
    std::filesystem::path config_path;

    std::vector<amsim::SampleDecl> sample_decl;
    std::vector<amsim::SampleEstimatorDecl> sample_estimator_decl;

    int opt;

    // NOLINTBEGIN(modernize-use-designated-initializers)
    struct option long_opts[] = {
        // Basic options
        {"help", no_argument, nullptr, 'h'},
        {"version", no_argument, nullptr, 'v'},
        // global options
        {"n-individuals", required_argument, nullptr, GlobalNumIndividuals},
        {"n-generations", required_argument, nullptr, GlobalNumGenerations},
        {"n-replicates", required_argument, nullptr, GlobalNumReplicates},
        {"n-threads", required_argument, nullptr, GlobalNumThreads},
        {"random-seed", required_argument, nullptr, GlobalRandomSeed},
        {"out-dir", required_argument, nullptr, GlobalOutputDirectory},
        {"out-name", required_argument, nullptr, GlobalOutputName},
        {"output-dir", required_argument, nullptr, GlobalOutputDirectory},
        {"output-name", required_argument, nullptr, GlobalOutputName},
        {"log-to-output", no_argument, nullptr, GlobalLogNoFile},
        {"log-level", required_argument, nullptr, GlobalLogLevel},
        {"save-config", no_argument, nullptr, GlobalSaveConfig},
        {"share-init-state", no_argument, nullptr, GlobalShareInitState},
        {"no-run", no_argument, nullptr, GlobalNoRun},
        // Genome options
        {"loc-maf", no_argument, nullptr, GenomeLocusInitMAFs},
        {"loc-rec", no_argument, nullptr, GenomeLocusRecombinationProbs},
        {"loc-mut", no_argument, nullptr, GenomeLocusMutationProbs},
        {"locus-maf", no_argument, nullptr, GenomeLocusInitMAFs},
        {"locus-rec", no_argument, nullptr, GenomeLocusRecombinationProbs},
        {"locus-mut", no_argument, nullptr, GenomeLocusMutationProbs},
        // Phenotype options
        {"phenotype", required_argument, nullptr, PhenotypeDecl},
        {"var-genetic", required_argument, nullptr, PhenotypeVarGenetic},
        {"var-environmental",
         required_argument,
         nullptr,
         PhenotypeVarEnvironmental},
        {"var-vertical", required_argument, nullptr, PhenotypeVarVertical},
        {"effects", no_argument, nullptr, PhenotypeLocusEffects},
        {"loci", no_argument, nullptr, PhenotypeLocusIndices},
        {"pheno-gen-cor", no_argument, nullptr, PhenotypeGeneticCor},
        {"pheno-env-cor", no_argument, nullptr, PhenotypeEnvironmentalCor},
        {"pheno-genetic-cor", no_argument, nullptr, PhenotypeGeneticCor},
        {"pheno-environmental-cor",
         no_argument,
         nullptr,
         PhenotypeEnvironmentalCor},
        // Mating options
        {"mating", required_argument, nullptr, MatingDecl},
        {"mate-cor", no_argument, nullptr, MatingCor},
        {"tol-inf", required_argument, nullptr, MatingErrorTolerance},
        {"max-itr", required_argument, nullptr, MatingMaxIterations},
        {"temp-init", required_argument, nullptr, MatingAnnealingTempInit},
        {"temp-decay", required_argument, nullptr, MatingAnnealingTempDecay},
        // Population estimator options
        {"estimator", required_argument, nullptr, PopulationEstimatorDecl},
        // Sample options
        {"sample", required_argument, nullptr, SampleDecl},
        {"proband", required_argument, nullptr, SampleProbandType},
        {"n-probands", required_argument, nullptr, SampleNumProbands},
        {"weight", required_argument, nullptr, SampleWeightFunction},
        {"on", required_argument, nullptr, SampleWeightOnPhenotypes},
        {"of", required_argument, nullptr, SampleWeightOfMembers},
        {"agg", required_argument, nullptr, SampleWeightAggregation},
        {"estimators", required_argument, nullptr, SampleEstimators},
        // Sample estimator options
        {"sample-estimator", required_argument, nullptr, SampleEstimatorDecl},
        {"type", required_argument, nullptr, SampleEstimatorType},
        {"exec", required_argument, nullptr, SampleEstimatorExec},
        {"n-rows", required_argument, nullptr, SampleEstimatorNumRows},
        {"n-cols", required_argument, nullptr, SampleEstimatorNumCols},
        {"row-names", required_argument, nullptr, SampleEstimatorRowNames},
        {"col-names", required_argument, nullptr, SampleEstimatorColNames},
        // Virtual (context-sensitive) options
        {"n-loci", required_argument, nullptr, VirtualNumLoci},
        {"file", required_argument, nullptr, VirtualFile},
        {"dist", required_argument, nullptr, VirtualDistribution},
        {"value", required_argument, nullptr, VirtualValue},
        {"singular-values", required_argument, nullptr, VirtualSingularValues},
        {nullptr, 0, nullptr, 0}};
    // NOLINTEND(modernize-use-designated-initializers)

    auto check_context = [&](Context expected, std::string_view flag) {
      if (context_domain(context) != expected)
        throw std::runtime_error(
            std::format("Unexpected flag '{}' in current context", flag));
    };

    // NOLINTBEGIN(bugprone-switch-missing-default-case)
    while ((opt = getopt_long(argc, argv, "hv", long_opts, nullptr)) != -1) {
      // Default options
      switch (opt) {
        case 'v':
          display_version();
          exit(EXIT_SUCCESS);
        case '?':
          throw std::runtime_error(
              std::format("Unrecognised option '{}'", argv[optind - 1]));
        case 'h':
          context = Context::Help;
          display_header();
          display_help();
          exit(EXIT_SUCCESS);
      }

      // Global options
      switch (opt) {
        case Option::GlobalNumIndividuals:
          context = Context::Global;
          simulation.n_individuals = std::stoull(optarg);
          continue;
        case Option::GlobalNumGenerations:
          context = Context::Global;
          simulation.n_generations = std::stoull(optarg);
          continue;
        case Option::GlobalNumThreads:
          context = Context::Global;
          n_threads = std::stoull(optarg);
          continue;
        case Option::GlobalNumReplicates:
          context = Context::Global;
          n_replicates = std::stoull(optarg);
          continue;
        case Option::GlobalRandomSeed:
          context = Context::Global;
          simulation.random_seed = std::stoull(optarg);
          continue;
        case Option::GlobalOutputDirectory:
          context = Context::Global;
          simulation.output_dir = optarg;
          continue;
        case Option::GlobalOutputName:
          context = Context::Global;
          simulation.output_name = optarg;
          continue;
        case Option::GlobalLogLevel:
          context = Context::Global;
          simulation.log_level = amsim::LogLevel_from_string(optarg);
          continue;
        case Option::GlobalLogNoFile:
          context = Context::Global;
          simulation.log_to_file = false;
          continue;
        case Option::GlobalSaveConfig:
          context = Context::Global;
          save_config = true;
          continue;
        case Option::GlobalShareInitState:
          context = Context::Global;
          simulation.share_init_state = true;
          continue;
        case Option::GlobalNoRun:
          context = Context::Global;
          run = false;
          continue;
      }

      // Genome configuration
      switch (opt) {
        case Option::GenomeLocusInitMAFs:
          context = Context::GenomeInitMAFs;
          continue;
        case Option::GenomeLocusRecombinationProbs:
          context = Context::GenomeRecombinationProbs;
          continue;
        case Option::GenomeLocusMutationProbs:
          context = Context::GenomeMutationProbs;
          continue;
      }

      // Phenome configuration
      switch (opt) {
        case Option::PhenotypeDecl: {
          context = Context::Phenotype;
          amsim::Phenotype pheno_new{.name = optarg};
          simulation.phenotypes.push_back(pheno_new);
          continue;
        }
        case Option::PhenotypeVarGenetic:
          check_context(Context::Phenotype, "--var-genetic");
          simulation.phenotypes.back().var_genetic = std::stod(optarg);
          continue;
        case Option::PhenotypeVarEnvironmental:
          check_context(Context::Phenotype, "--var-environmental");
          simulation.phenotypes.back().var_environmental = std::stod(optarg);
          continue;
        case Option::PhenotypeVarVertical:
          check_context(Context::Phenotype, "--var-vertical");
          simulation.phenotypes.back().var_vertical = std::stod(optarg);
          continue;
        case Option::PhenotypeGeneticCor:
          context = Context::PhenotypeGeneticCorrelation;
          continue;
        case Option::PhenotypeEnvironmentalCor:
          context = Context::PhenotypeEnvironmentalCorrelation;
          continue;
        case Option::PhenotypeLocusEffects:
          check_context(Context::Phenotype, "--effects");
          context = Context::PhenotypeCausalLoci;
          continue;
        case Option::PhenotypeLocusIndices:
          check_context(Context::Phenotype, "--loci");
          context = Context::PhenotypeCausalLoci;
          continue;
      }

      // Mating configuration
      switch (opt) {
        case Option::MatingDecl:
          if (std::string_view(optarg) == "random")
            simulation.mating.type = "random";
          if (std::string_view(optarg) == "assortative") {
            context = Context::Mating;
            simulation.mating.type = "assortative";
          }
          continue;
        case Option::MatingCor:
          check_context(Context::Mating, "--mate-cor");
          context = Context::MatingCorrelation;
          continue;
        case Option::MatingErrorTolerance:
          check_context(Context::Mating, "--tol-inf");
          simulation.mating.tolerance = std::stod(optarg);
          continue;
        case Option::MatingMaxIterations:
          check_context(Context::Mating, "--max-itr");
          simulation.mating.max_iterations = std::stoull(optarg);
          continue;
        case Option::MatingAnnealingTempInit:
          check_context(Context::Mating, "--temp-init");
          simulation.mating.initial_temperature = std::stod(optarg);
          continue;
        case Option::MatingAnnealingTempDecay:
          check_context(Context::Mating, "--temp-decay");
          simulation.mating.temperature_decay = std::stod(optarg);
          continue;
      }

      // Estimators and sampling
      switch (opt) {
        case Option::PopulationEstimatorDecl: {
          context = Context::PopulationEstimator;

          auto [name, params] = amsim::parse_function(optarg);
          auto it = estimators.find(name);

          if (it == estimators.end())
            throw std::runtime_error(
                "Could not find population estimator " + name);

          simulation.estimators.push_back(estimators[name](params));

          continue;
        }
        case Option::SampleEstimatorDecl: {
          context = Context::SampleEstimator;
          amsim::SampleEstimatorDecl sample_est_desc{.name = optarg};
          sample_estimator_decl.push_back(sample_est_desc);
          continue;
        }
        case Option::SampleEstimatorType: {
          check_context(Context::SampleEstimator, "--type");
          auto [name, params] = amsim::parse_function(optarg);
          sample_estimator_decl.back().type = name;
          if (!params.empty()) sample_estimator_decl.back().params = params;
          continue;
        }
        case Option::SampleEstimatorExec:
          check_context(Context::SampleEstimator, "--exec");
          sample_estimator_decl.back().exec = optarg;
          continue;
        case Option::SampleEstimatorNumRows:
          check_context(Context::SampleEstimator, "--n-rows");
          sample_estimator_decl.back().n_rows = std::stoull(optarg);
          continue;
        case Option::SampleEstimatorNumCols:
          check_context(Context::SampleEstimator, "--n-cols");
          sample_estimator_decl.back().n_cols = std::stoull(optarg);
          continue;
        case Option::SampleEstimatorRowNames:
          check_context(Context::SampleEstimator, "--row-names");
          sample_estimator_decl.back().row_names =
              amsim::utils::split_string(optarg);
          continue;
        case Option::SampleEstimatorColNames:
          check_context(Context::SampleEstimator, "--col-names");
          sample_estimator_decl.back().col_names =
              amsim::utils::split_string(optarg);
          continue;
        case Option::SampleDecl: {
          context = Context::Sample;
          amsim::SampleDecl sample_desc;
          sample_desc.name = optarg;
          sample_decl.push_back(sample_desc);
          continue;
        }
        case Option::SampleProbandType:
          sample_decl.back().proband_type = optarg;
          continue;
        case Option::SampleNumProbands:
          check_context(Context::Sample, "--n-probands");
          sample_decl.back().n_probands = std::stoull(optarg);
          continue;
        case Option::SampleWeightOnPhenotypes: {
          check_context(Context::Sample, "--on");
          sample_decl.back().on = amsim::utils::split_string(optarg);
          continue;
        }
        case Option::SampleWeightOfMembers:
          check_context(Context::Sample, "--of");
          sample_decl.back().of = amsim::utils::split_string(optarg);
          continue;
        case Option::SampleWeightAggregation:
          check_context(Context::Sample, "--agg");
          sample_decl.back().agg = optarg;
          continue;
        case Option::SampleWeightFunction:
          check_context(Context::Sample, "--weight");
          sample_decl.back().weight_function = optarg;
          continue;
        case Option::SampleEstimators:
          check_context(Context::Sample, "--estimators");
          sample_decl.back().estimators = amsim::utils::split_string(optarg);
          continue;
      }

      // Context-sensitive operation flags
      Context domain = context_domain(context);

      switch (opt) {
        case Option::VirtualNumLoci: {
          std::size_t n_loc = std::stoull(optarg);
          if (domain == Context::Phenotype)
            simulation.phenotypes.back().n_causal_loci = n_loc;
          else
            simulation.genome.n_loci = n_loc;
          continue;
        }
        case Option::VirtualValue: {
          if (context == PhenotypeCausalLoci) {
            std::vector<std::size_t> value = amsim::parse_indices(optarg);
            simulation.phenotypes.back().causal_loci = value;
            continue;
          }

          Eigen::MatrixXd matrix = amsim::parse_matrix_value(optarg);
          if (context == Context::GenomeInitMAFs)
            simulation.genome.v_maf = matrix(0);
          if (context == Context::GenomeRecombinationProbs)
            simulation.genome.v_rec = matrix(0);
          if (context == Context::GenomeMutationProbs)
            simulation.genome.v_mut = matrix(0);
          if (context == Context::PhenotypeGeneticCorrelation)
            simulation.genetic_component_cor = matrix;
          if (context == Context::PhenotypeEnvironmentalCorrelation)
            simulation.environmental_component_cor = matrix;
          if (context == Context::MatingCorrelation)
            simulation.mating.mate_cor = matrix;
          continue;
        }
        case Option::VirtualFile: {
          if (context == PhenotypeCausalLoci) {
            auto file = amsim::File<std::vector<std::size_t>>{.path = optarg};
            simulation.phenotypes.back().causal_loci = file;
            continue;
          }

          auto file = amsim::File<Eigen::MatrixXd>{.path = optarg};
          if (context == Context::GenomeInitMAFs)
            simulation.genome.v_maf = file;
          if (context == Context::GenomeRecombinationProbs)
            simulation.genome.v_rec = file;
          if (context == Context::GenomeMutationProbs)
            simulation.genome.v_mut = file;
          if (context == Context::PhenotypeGeneticCorrelation)
            simulation.genetic_component_cor = file;
          if (context == Context::PhenotypeEnvironmentalCorrelation)
            simulation.environmental_component_cor = file;
          if (context == Context::MatingCorrelation)
            simulation.mating.mate_cor = file;
          continue;
        }
        case Option::VirtualSingularValues: {
          if (context == Context::PhenotypeGeneticCorrelation) {
            Eigen::MatrixXd matrix =
                amsim::utils::matrix_from_singular_values(optarg, true);
            simulation.genetic_component_cor = matrix;
          }
          if (context == Context::PhenotypeEnvironmentalCorrelation) {
            Eigen::MatrixXd matrix =
                amsim::utils::matrix_from_singular_values(optarg, true);
            simulation.environmental_component_cor = matrix;
          }
          if (context == Context::MatingCorrelation) {
            Eigen::MatrixXd matrix =
                amsim::utils::matrix_from_singular_values(optarg);
            simulation.mating.mate_cor = matrix;
          }
          continue;
        }
        case Option::VirtualDistribution: {
          amsim::Distribution dist = amsim::parse_distribution(
              optarg, (context & Context::GenomeProbabilities) != 0U);

          if (context == Context::GenomeInitMAFs)
            simulation.genome.v_maf = dist;
          if (context == Context::GenomeRecombinationProbs)
            simulation.genome.v_rec = dist;
          if (context == Context::GenomeMutationProbs)
            simulation.genome.v_mut = dist;
          if (context_domain(context) == Context::Phenotype)
            simulation.phenotypes.back().effects = dist;
          continue;
        }
      }
      // NOLINTEND(bugprone-switch-missing-default-case)
    }

    // validate sample estimator declarations
    for (const auto& se : sample_estimator_decl) {
      if (se.type == "external") {
        if (!se.exec.has_value())
          throw std::runtime_error(
              std::format(
                  "sample estimator '{}': type 'external' requires --exec",
                  se.name));
        if (!se.n_rows.has_value() || !se.n_cols.has_value())
          throw std::runtime_error(
              std::format(
                  "sample estimator '{}': type 'external' requires --n-rows "
                  "and --n-cols",
                  se.name));
      } else if (se.exec.has_value()) {
        throw std::runtime_error(
            std::format(
                "sample estimator '{}': --exec is only valid for type "
                "'external'",
                se.name));
      }
    }

    // error if any sample has no estimators
    for (const auto& s : sample_decl)
      if (s.estimators.empty())
        throw std::runtime_error(
            std::format("sample '{}' has no estimators declared", s.name));

    std::vector<amsim::SampleSpec> samples =
        amsim::build_samples(sample_decl, sample_estimator_decl);

    simulation.samples = std::move(samples);

    if (save_config) {
      if (simulation.output_name.has_value())
        config_path = simulation.output_dir /
                      ("config_" + simulation.output_name.value() + ".toml");
      else
        config_path = simulation.output_dir / "config.toml";

      amsim::write_config(
          simulation,
          n_replicates,
          n_threads,
          sample_decl,
          sample_estimator_decl,
          config_path);
    }

    if (run) amsim::run_simulation(simulation, n_replicates, n_threads);
    exit(EXIT_SUCCESS);
  } catch (const std::exception& e) {
    amsim::Log::error(e.what());
    exit(EXIT_FAILURE);
  }
}
