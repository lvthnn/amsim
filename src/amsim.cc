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
#include <amsim/io.h>
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
  MatrixSpecification = 1 << 13,

  // Derived contexts
  GenomeInitMAFs = Genome | InitMAFs,
  GenomeRecombinationProbs = Genome | RecombinationProbs,
  GenomeMutationProbs = Genome | MutationProbs,
  GenomeProbabilities = InitMAFs | RecombinationProbs | MutationProbs,
  PhenotypeCausalLoci = Phenotype | CausalLoci,
  PhenotypeGeneticCorrelation =
      Phenotype | GeneticComponent | MatrixSpecification,
  PhenotypeEnvironmentalCorrelation =
      Phenotype | EnvironmentalComponent | MatrixSpecification,
  MatingCorrelation = Mating | MatrixSpecification,

  // Base domain
  ContextDomain = Global | Genome | Phenotype | Mating | Sample |
                  SampleEstimator | PopulationEstimator,
};

Context context_domain(Context context) {
  return static_cast<Context>(context & ContextDomain);
}

enum Option {
  GlobalNumIndividuals,
  GlobalNumGenerations,
  GlobalNumThreads,
  GlobalNumReplicates,
  GlobalRandomSeed,
  GlobalOutputDirectory,
  GlobalOutputName,
  GlobalLogLevel,
  GlobalLogNoFile,
  GlobalPedigreeMaxDepth,
  GlobalPedigreeWarmup,
  GlobalSaveConfig,
  GlobalLoadConfig,
  GlobalShareInitState,
  GlobalNoRun,
  GenomeLocusInitMAFs,
  GenomeLocusRecombinationProbs,
  GenomeLocusMutationProbs,
  PhenotypeSpec,
  PhenotypeVarGenetic,
  PhenotypeVarEnvironmental,
  PhenotypeVarVertical,
  PhenotypeLocusEffects,
  PhenotypeLocusIndices,
  PhenotypeGeneticCor,
  PhenotypeEnvironmentalCor,
  MatingSpec,
  MatingCor,
  MatingErrorTolerance,
  MatingMaxIterations,
  MatingAnnealingTempInit,
  MatingAnnealingTempDecay,
  SampleSpec,
  SampleProbandType,
  SampleNumProbands,
  SampleWeightFunction,
  SampleWeightOnPhenotypes,
  SampleWeightOfMembers,
  SampleWeightAggregation,
  SampleEstimators,
  SampleEstimatorSpec,
  SampleEstimatorType,
  SampleEstimatorExec,
  SampleEstimatorNumRows,
  SampleEstimatorNumCols,
  SampleEstimatorRowNames,
  SampleEstimatorColNames,
  PopulationEstimatorSpec,
  VirtualNumLoci,  // these options are context-sensitive
  VirtualFile,
  VirtualValue,
  VirtualDistribution,
  VirtualSingularValues
};

std::string phenotype_docs() {
  return R"(
--phenotype <name> declares a new phenotype and opens its context. Every
flag below, up to the next --phenotype (or a flag from a different domain),
configures that phenotype. Flags --pheno-gen-cor and --pheno-env-cor assign
correlation structure for the genetic and environmental components of traits,
respectively.

  --phenotype <name>
    --n-loci <n_loc>
    --loci [--value <idx> | --file <path>]
    --effects [--value <val> | --file <path> | --dist <distribution>]
    --var-genetic <var_gen>
    --var-environmental <var_env>
    --var-vertical <var_vert>
  --pheno-gen-cor [--value <vals> | --file <path> | --singular-values <sv>]
  --pheno-env-cor [--value <vals> | --file <path> | --singular-values <sv>]

A phenotype is the sum of up to three independent components:

    phenotype = genetic + environmental + vertical

--var-genetic, --var-environmental and --var-vertical set the share of
variance contributed by each; they need not sum to 1, since they are
renormalised internally so that the resulting heritabilities sum to one.
The genetic component is a linear combination of allele dosages at <n_loc>
causal loci, weighted by an effect-size vector; that vector is rescaled
internally so the genetic component alone carries variance h2_genetic. The
environmental component is drawn from an iid normal distribution for each
individual.

  --loci
    Which of the genome's <n_loci> loci are causal for this phenotype, given
    as indices, either directly via --value, or loaded from a file via --file.
    If omitted, <n_loc> loci are chosen at random. Mutually exclusive with
    --pheno-gen-cor / --pheno-env-cor below.

  --effects
    Per-locus effect sizes for the genetic component, over the causal loci above.
    Either a constant (--value), a vector loaded from file (--file), or drawn from
    a distribution via --dist, e.g. normal(0,1). Defaults to a Rademacher
    distribution (+-1 with equal probability) if omitted.

With more than one phenotype declared, --pheno-gen-cor and --pheno-env-cor set
the cross-phenotype correlation of the genetic and environmental components
respectively -- e.g. two phenotypes sharing part of their genetic basis (pleiotropy)
without either being causal for the other's loci. Each is a P x P correlation matrix
over the P declared phenotypes, read symmetrically (component correlation, not
directional):

             height weight
    height [ 1.00   0.50 ]
    weight [ 0.50   1.00 ]

Supply it directly (--value / --file), or generate a random one with a prescribe
eigenspectrum via --singular-values (P values in, a random symmetric P x P matrix
with exactly those eigenvalues out). Defaults to the identity matrix (uncorrelated
components) if not specified.

For more on the --value / --file / --dist / --singular-values flags themselves,
run --help on each.
  )";
}

std::string mating_docs() {
  return R"(
--mating {random | assortative} selects how mates are paired each generation.

  random
    Mates are drawn uniformly at random from the opposite sex. No sub-flags
    apply. This is the default if --mating is never specified.

  assortative
    Mates are paired so that the phenotypic correlation among mated pairs
    matches a target correlation matrix, --mate-cor.

    --mate-cor [--value <vals> | --file <path> | --singular-values <sv>]
      The target correlation matrix among mates, over the P declared
      phenotypes. Entry (i, j) is the correlation between the husband's
      phenotype i and the wife's phenotype j; unlike --pheno-gen-cor /
      --pheno-env-cor, this matrix is not forced symmetric, so cross-trait
      and directional assortment are both expressible:

                          wife:
                          height weight
          husband height [ 0.40  0.10 ]
                  weight [ 0.05  0.20 ]

      A single phenotype reduces this to the familiar scalar spousal correlation.

    --tol-inf <tol_inf>
      Convergence tolerance on the empirical-vs-target mate correlation
      (matrix norm of the difference). Default 1e-7.

    --max-itr <max_itr>
      Maximum number of annealing proposals attempted before giving up and
      returning the best pairing found so far. Default 2,000,000.

    --temp-init <temp_init>
    --temp-decay <temp_decay>
      Initial temperature and per-iteration decay of the mate matching algorithm.
      A higher --temp-init (default: 1.0) accepts more correlation-worsening swaps
      early on, helping escape local optima at the cost of slower convergence;
      --temp-decay (default 0.99) controls how quickly that exploration cools into
      pure hill-climbing. Each iteration of the algorithm, the temperature is
      updated via temp *= temp_decay.
  )";
}

std::string sample_docs() {
  return R"(
amsim models ascertainment / participation bias by drawing a sub-sample of
the population with inclusion probability governed by a weighting function of
(possibly aggregated) phenotype data, rather than uniformly. Two related
declarations make this up:

  --sample <name>
    --proband {individual | mate | family}
    --n-probands <n_prob>
    --on <phenotype(s)>
    --of <member(s)>
    --agg {mean | max | min | identity}
    --weight {uniform() | logistic(<coef>, ...)}
    --estimators <est_name(s)>

  --sample-estimator <name>
    --type {gwas(<n_pcs>, <pval_thresh>) | haseman-elston | greml | external |
            sample-mean | sample-var | sample-cov | sample-mate-cor}
    --exec <exec_cmd>                                      (type external only)
    --n-rows <n_rows>                                      (type external only)
    --n-cols <n_cols>                                      (type external only)
    --row-names <name(s)>                                  (type external only)
    --col-names <name(s)>                                  (type external only)

--proband fixes the unit of sampling and, with it, which family members are
addressable through --of:

  individual   a single person (--of is not applicable)
  mate         a mated couple: husband, wife, and their parents-in-law. Valid
                options are {'husband', 'wife', 'husbandfather', 'husbandmother',
                'wifefather', 'wifemother', 'couple', 'husbandinlaws', 'wifeinlaws',
                'parents', 'husbandfamily', 'wifefamily', 'males', 'females', 'all'}.
  family       a nuclear family: father, mother, their children, and whoever their
                children married. Valid options are {'father', 'mother', 'son',
                'daughter', 'sonwife', 'daughterhusband', 'parents', 'siblings',
                'males', 'females', 'all'}.

For each proband, --on selects which phenotype(s) enter the weighting
calculation and --of selects whose values are used; --agg then reduces the
values across the selected members down to one number per phenotype per
proband (mean/max/min, or identity for a single member). --weight turns
that aggregated vector into an inclusion probability:

  uniform()             every proband is equally likely to be sampled -- i.e.,
                         no ascertainment bias, a random sample.
  logistic(b1, b2, ...) inclusion probability = 1 / (1 + exp(-agg * b)),
                         a logistic model in the aggregated phenotype
                         values; one coefficient per value selected via
                         --on, so ascertainment can depend on, and bias
                         estimates with respect to, one or more phenotypes.

amsim then draws --n-probands individuals under that probability and attaches
every sample estimator named in --estimators (each declared separately via its
own --sample-estimator block) to the resulting sub-sample, so the same estimator
machinery used at the population level can be recomputed under whatever
ascertainment or sampling scheme is being studied.

--sample-estimator works like --estimator (see --help --estimator) but is declared
by name so it can be referenced from one or more --sample blocks via --estimators;
--type additionally supports sample-only estimators (sample-mean, sample-var,
sample-cov, sample-mate-cor) and an external escape hatch that shells out to
--exec with the sample matrix (--n-rows x --n-cols) written to a temporary file.
  )";
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

void display_help() {
  std::string help_str = R"(
For full documentation and usage examples, see https://github.com/lvthnn/amsim.
Running amsim --help <flag> prints usage guides on the various command flags
listed below. Example: amsim --help --n-individuals.

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
  --pedigree-max-depth <max_depth>
  --pedigree-warmup
  --log-level <log_level>
  --log-to-output
  --save-config
  --load-config <config_path>
  --share-init-state
  --no-run

genome options:
  --locus-maf [--val <value> | --file <path> | --dist <dist>]
  --locus-rec [--val <value> | --file <path> | --dist <dist>]
  --locus-mut [--val <value> | --file <path> | --dist <dist>]

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
  --estimator {genotype-mean | genotype-var | genotype-maf | genotype-cov |
               genotype-cor | heritability | pheno-mean(<component>) |
               pheno-var(<component>) | pheno-cov(<component>) |
               mate-cor(<component>) | cousin-cov(<component>) |
               ancestor-cov(<component>)}
  --sample-estimator <name>
    --type {gwas(<n_pcs>, <pval_thresh>) | haseman-elston | greml | external |
            sample-mean | sample-var | sample-cov | sample-mate-cor}
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
    --agg <agg_fn>
    --estimators <est_names>
  )";

  std::cout << help_str << std::endl;
}

void display_help_short() {
  std::string help_str = R"(
  amsim <option(s)>
  amsim --load-config <config_file> <option(s)>

To see all options, run "amsim --help".)";

  std::cout << help_str << std::endl;
}

void display_help_option(Option option, const std::string& option_flag) {
  std::unordered_map<Option, std::string> docs;

  docs[Option::GlobalNumIndividuals] = R"(
Number of individuals to simulate in each generation. Positive, non-zero integer.
  )";

  docs[Option::GlobalNumGenerations] = R"(
Number of generations to simulate. Positive, non-zero integer.
  )";

  docs[Option::GlobalNumThreads] = R"(
Number of threads to use for multithreaded replicate simulation. Positive, non-zero
integer.
  )";

  docs[Option::GlobalNumReplicates] = R"(
Number of replicate simulations to perform. Can be started from either the same
initial population or different ones based on whether the --share-init-state flag
is specified. Positive, non-zero integer.
  )";

  docs[Option::GlobalRandomSeed] = R"(
Random seed for reproducible results. Positive, non-zero integer.
  )";

  docs[Option::GlobalOutputDirectory] = R"(
Path to output directory of generated results.
  )";

  docs[Option::GlobalOutputName] = R"(
Basename used for files produced during simulation.
  )";

  docs[Option::GlobalLogLevel] = R"(
Specifies the verbosity of amsim's internal logger system. Options are 'none',
'error', 'warning', 'info', 'debug'.
  )";

  docs[Option::GlobalLogNoFile] = R"(
If enabled, amsim logs directly to the shell instead of writing to file.
  )";

  docs[Option::GlobalPedigreeMaxDepth] = R"(
Maximum depth of pedigree structure.
  )";

  docs[Option::GlobalPedigreeWarmup] = R"(
If enabled, simulates panmictic generations before generating the founder sibling
pairs so that the pedigree is at full capacity immediately in the first generation.
This option is useful when employing kinship-based estimators so estimates are
available immediately from the initial generation. Otherwise, the estimators
return NA until the pedigree is filled to its desired depth.
  )";

  docs[Option::GlobalSaveConfig] = R"(
If enabled, saves configuration to a TOML file for reuse. Can be paired with the
--no-run flag to generate a configuration file without performing simulation.
  )";

  docs[Option::GlobalLoadConfig] = R"(
Attempts to load configuration file from the provided path. Parameters can be
overloaded by supplying them in conjunction with this option. In the case that a
declarative option is specified, such as adding a phenotype, then this is added
to the existing configuration.
  )";

  docs[Option::GlobalShareInitState] = R"(
If enabled, starts simulations from the same founder population and state.
  )";

  docs[Option::GlobalNoRun] = R"(
If enabled, the software exits without performing simulation. This is useful when
specifying configurations to be saved to file via --save-config.
  )";

  docs[Option::GenomeLocusInitMAFs] = R"(
Specify the minor allele frequencies (MAFs) of the genetic loci in the simulation.
In a population numbering n individuals, the MAF of a locus is the total dosage of
the minor allele divided by the total number of alleles for the locus.

The MAFs can be specified through directly through the --value flag, in which case
<n_loci> double values in the range [0.0, 1.0] must be specified. On the other hand,
these data may be supplied through the --file option, which takes as argument a path
to a file. Lastly, the --dist flag can be used to specify a probability distribution
from which these frequencies are drawn.

For more information on virtual flags and their use, run --help on flags --value,
--file, --dist, and --singular-values.
  )";

  docs[Option::GenomeLocusRecombinationProbs] = R"(
Specify the recombination probabilities of the genetic loci in the simulation.
Each entry gives the probability that a gamete formed at meiosis recombines at
that position rather than continuing to inherit from the same parental strand. A
value of 0.5 corresponds to free recombination or independence between adjoining
loci, whereas a value near 0.0 corresponds to tight linkage.

As with --loc-maf, these can be specified directly through the --value flag
(<n_loci> values in [0.0, 1.0]), loaded from a file via --file, or drawn from a
distribution via --dist. Since these are probabilities, --dist is restricted to
[0,1]-supported distributions (uniform, beta).

For more information on virtual flags and their use, run --help on flags --value,
--file, and --dist.
  )";

  docs[Option::GenomeLocusMutationProbs] = R"(
Specify the per-locus mutation probabilities, i.e., the chance that a given allele
is flipped to its complement when passed from parent to offspring. One value per
locus, in [0.0, 1.0], with default 0.0 (no mutation).

As with --loc-maf and --loc-rec, these can be specified directly through the
--value flag, loaded from a file via --file, or drawn from a distribution via
--dist, restricted to [0,1]-supported distributions (uniform, beta) since these
are probabilities.

For more information on virtual flags and their use, run --help on flags --value,
--file, and --dist.
  )";

  docs[Option::PhenotypeSpec] = phenotype_docs();
  docs[Option::PhenotypeVarGenetic] = phenotype_docs();
  docs[Option::PhenotypeVarEnvironmental] = phenotype_docs();
  docs[Option::PhenotypeVarVertical] = phenotype_docs();
  docs[Option::PhenotypeLocusEffects] = phenotype_docs();
  docs[Option::PhenotypeLocusIndices] = phenotype_docs();
  docs[Option::PhenotypeGeneticCor] = phenotype_docs();
  docs[Option::PhenotypeEnvironmentalCor] = phenotype_docs();

  docs[Option::MatingSpec] = mating_docs();
  docs[Option::MatingCor] = mating_docs();
  docs[Option::MatingErrorTolerance] = mating_docs();
  docs[Option::MatingMaxIterations] = mating_docs();
  docs[Option::MatingAnnealingTempInit] = mating_docs();
  docs[Option::MatingAnnealingTempDecay] = mating_docs();

  docs[Option::SampleSpec] = sample_docs();
  docs[Option::SampleProbandType] = sample_docs();
  docs[Option::SampleNumProbands] = sample_docs();
  docs[Option::SampleWeightFunction] = sample_docs();
  docs[Option::SampleWeightOnPhenotypes] = sample_docs();
  docs[Option::SampleWeightOfMembers] = sample_docs();
  docs[Option::SampleWeightAggregation] = sample_docs();

  docs[Option::SampleEstimators] = sample_docs();
  docs[Option::SampleEstimatorSpec] = sample_docs();
  docs[Option::SampleEstimatorType] = sample_docs();
  docs[Option::SampleEstimatorNumRows] = sample_docs();
  docs[Option::SampleEstimatorNumCols] = sample_docs();
  docs[Option::SampleEstimatorRowNames] = sample_docs();
  docs[Option::SampleEstimatorColNames] = sample_docs();

  docs[Option::PopulationEstimatorSpec] = R"(
--estimator <name(<args>)> declares a population-wide estimator, or a statistic
computed once per generation over the full population. Repeat the flag to add more
than one.

  --estimator {genotype-mean | genotype-var | genotype-maf | genotype-cov |
               genotype-cor | heritability | pheno-mean(<component>) |
               pheno-var(<component>) | pheno-cor(<component>) |
               pheno-cov(<component>) | mate-cor(<component>) |
               cousin-cov(<component>, <degree>) |
               ancestor-cov(<component>, <degree>)}

genotype-* and heritability take no arguments and operate on the genome buffer /
true heritability directly. The rest operate on a phenotype component and take an
optional <component>, one of genetic, environmental, vertical, or total. The
default component selected is total.

cousin-cov and ancestor-cov additionally take an optional <degree> (default: 1):
the pedigree distance to compute the covariance at degree 1 for first cousins /
grandparent-grandchild, degree 2 for second cousins / great-grandparents, and so
on. Resolving degree d requires a pedigree deep enough to reach it; see
--pedigree-max-depth and --pedigree-warmup.

Examples:
  --estimator genotype-maf
  --estimator heritability
  --estimator 'pheno-cor(genetic)'
  --estimator 'cousin-cov(total,2)'

Quoting matters here: unescaped, the shell will try to interpret the parentheses
itself rather than passing them through to amsim, so wrap any functional-syntax
argument in quotes.
  )";

  docs[Option::VirtualValue] = R"(
--value <val(s)> supplies data inline on the command line. What is expected
depends on which flag opened the current context:

  --loc-maf / --loc-rec / --loc-mut   a single value (broadcast to every locus)
                                       or one value per locus
  --loci                              whitespace- or comma-separated causal
                                       locus indices
  --effects                           a single value (broadcast to every causal
                                       locus) or one value per locus
  --pheno-gen-cor / --pheno-env-cor   a symmetric P x P matrix, entered as rows
  / --mate-cor                         separated by ';' or a newline, values
                                       within a row separated by a comma or
                                       space; a single row collapses to a vector

Example: --loc-maf --value 0.25                       (one MAF for every locus)
         --pheno-gen-cor --value "1,0.5;0.5,1"           (an inline 2x2 matrix)
  )";

  docs[Option::VirtualFile] = R"(
--file <path> loads data from a file. The format is plain text: values separated
by a comma or space, rows separated by a semicolon or newline; a file with a
single row is read as a vector rather than a 1xN matrix. Tabs are not a recognised
separator.

Which flag opens the context determines what the loaded data is used for, exactly
as with --value: locus MAFs / recombination / mutation probabilities, causal locus
indices, effect sizes, or a phenotype gen-cor / env-cor / mate-cor matrix.
  )";

  docs[Option::VirtualDistribution] = R"(
--dist <distribution> draws the requested quantity from a probability distribution
instead of supplying it directly, e.g. --dist 'normal(0,1)'. Note that the argument
must be quoted. Supported distributions and their parameters:

  uniform(lo, hi)
  beta(a, b)
  exponential(lambda)
  gamma(shape, scale)
  normal(mean, sd)
  laplace(location, scale)
  students_t(df)
  )";

  docs[Option::VirtualNumLoci] = R"(
--n-loci <n> sets a count that depends on where it appears. At global / genome
scope it sets the number of loci simulated across the whole genome. Inside a
--phenotype block, it instead sets that phenotype's number of causal loci
(defaulting to genome-loci / n-phenotypes if never given).
  )";

  docs[Option::VirtualSingularValues] = R"(
--singular-values <sv(s)> builds a correlation or cross-correlation matrix from
a prescribed spectrum instead of specifying every entry by hand: give P values
and get back a random P x P matrix with exactly that spectrum, which is
otherwise awkward to construct by hand while keeping the matrix a valid
correlation matrix.

For --pheno-gen-cor / --pheno-env-cor the result is forced symmetric (a
random eigendecomposition U*diag(sv)*U^T), since those are correlations
between components of the same kind. For --mate-cor the result is not
forced symmetric (a random singular value decomposition U*diag(sv)*V^T),
since husband-phenotype-i-to-wife-phenotype-j correlations need not equal
their transpose.
  )";

  std::string subtitle = "Help for option " + option_flag;
  std::cout << std::endl << subtitle << std::endl;
  std::cout << std::string(subtitle.length(), '-');
  std::cout << docs[option] << std::endl;
}

std::vector<option> get_options() {
  // NOLINTBEGIN(modernize-use-designated-initializers)
  return {
      // other options
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
      {"pedigree-max-depth",
       required_argument,
       nullptr,
       GlobalPedigreeMaxDepth},
      {"pedigree-warmup", no_argument, nullptr, GlobalPedigreeWarmup},
      {"log-to-output", no_argument, nullptr, GlobalLogNoFile},
      {"log-level", required_argument, nullptr, GlobalLogLevel},
      {"save-config", no_argument, nullptr, GlobalSaveConfig},
      {"load-config", required_argument, nullptr, GlobalLoadConfig},
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
      {"phenotype", required_argument, nullptr, PhenotypeSpec},
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
      {"mating", required_argument, nullptr, MatingSpec},
      {"mate-cor", no_argument, nullptr, MatingCor},
      {"tol-inf", required_argument, nullptr, MatingErrorTolerance},
      {"max-itr", required_argument, nullptr, MatingMaxIterations},
      {"temp-init", required_argument, nullptr, MatingAnnealingTempInit},
      {"temp-decay", required_argument, nullptr, MatingAnnealingTempDecay},

      // Population estimator options
      {"estimator", required_argument, nullptr, PopulationEstimatorSpec},

      // Sample options
      {"sample", required_argument, nullptr, SampleSpec},
      {"proband", required_argument, nullptr, SampleProbandType},
      {"n-probands", required_argument, nullptr, SampleNumProbands},
      {"weight", required_argument, nullptr, SampleWeightFunction},
      {"on", required_argument, nullptr, SampleWeightOnPhenotypes},
      {"of", required_argument, nullptr, SampleWeightOfMembers},
      {"agg", required_argument, nullptr, SampleWeightAggregation},
      {"estimators", required_argument, nullptr, SampleEstimators},

      // Sample estimator options
      {"sample-estimator", required_argument, nullptr, SampleEstimatorSpec},
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
}

int main(int argc, char* argv[]) {
  // get options
  std::vector<option> opts = get_options();

  // context helps us resolve what options are legal
  Context context = Context::Global;

  auto check_context = [&](Context expected, std::string_view flag) {
    if (context_domain(context) != expected)
      throw std::runtime_error(
          std::format("Unexpected flag '{}' in current context", flag));
  };

  // auxiliary boolean variables to manage CLI control flow
  bool run = true;
  bool save_config = false;

  // for options processing
  int opt;

  try {
    // invocation without arguments — display help and exit
    if (argc == 1) {
      display_header();
      display_help_short();
      exit(EXIT_SUCCESS);
    }

    // check for the --help flag; if specified, either display the full message
    // if no other arguments, or show option documentation
    for (int i = 1; i < argc; ++i) {
      if (std::string_view(argv[i]) == "--help" ||
          std::string_view(argv[i]) == "-h") {
        if (argc == 2) {
          display_header();
          display_help();
          exit(EXIT_SUCCESS);
        }
        if (argc > 3) {
          display_header();
          display_help_short();
          exit(EXIT_FAILURE);
        }

        Option opt;
        std::string opt_flag = argv[(i == 1) ? 2 : 1];
        auto pos = opt_flag.find_first_not_of('-');
        std::string opt_name = opt_flag.substr(pos);

        for (const auto& option : opts) {
          if (option.name == nullptr) break;
          if (option.name == opt_name) {
            opt = static_cast<Option>(option.val);
            display_header();
            display_help_option(opt, opt_flag);
            exit(EXIT_SUCCESS);
          }
        }

        // we've not found the option — show a short help message and exit with
        // failure
        display_header();
        display_help_short();
        exit(EXIT_FAILURE);
      }
    }

    // set up a configuration object; read from file if specified by user
    amsim::SimulationSpec spec;

    for (int i = 1; i < argc; ++i) {
      if (std::string_view(argv[i]) == "--load-config") {
        if ((i + 1) == argc) break;
        std::filesystem::path config_path = argv[i + 1];
        amsim::ConfigReader config(config_path);
        spec = config.result();
      }
    }

    // NOLINTBEGIN(bugprone-switch-missing-default-case)
    while ((opt = getopt_long(argc, argv, "hv", opts.data(), nullptr)) != -1) {
      // Default options
      switch (opt) {
        case 'v':
          display_version();
          exit(EXIT_SUCCESS);
        case 'h':
          // no need to do anything here since the help flag is resolved
          continue;
        case '?':
          throw std::runtime_error(
              std::format("Unrecognised option '{}'", argv[optind - 1]));
      }

      // Global options
      switch (opt) {
        case Option::GlobalNumIndividuals:
          context = Context::Global;
          spec.n_individuals = std::stoull(optarg);
          continue;
        case Option::GlobalNumGenerations:
          context = Context::Global;
          spec.n_generations = std::stoull(optarg);
          continue;
        case Option::GlobalNumThreads:
          context = Context::Global;
          spec.n_threads = std::stoull(optarg);
          continue;
        case Option::GlobalNumReplicates:
          context = Context::Global;
          spec.n_replicates = std::stoull(optarg);
          continue;
        case Option::GlobalRandomSeed:
          context = Context::Global;
          spec.random_seed = std::stoull(optarg);
          continue;
        case Option::GlobalOutputDirectory:
          context = Context::Global;
          spec.output_dir = optarg;
          continue;
        case Option::GlobalOutputName:
          context = Context::Global;
          spec.output_name = optarg;
          continue;
        case Option::GlobalPedigreeMaxDepth:
          context = Context::Global;
          spec.pedigree_max_depth = std::stoull(optarg);
          continue;
        case Option::GlobalPedigreeWarmup:
          context = Context::Global;
          spec.pedigree_warmup = true;
          continue;
        case Option::GlobalLogLevel:
          context = Context::Global;
          spec.log_level = amsim::LogLevel_from_string(optarg);
          continue;
        case Option::GlobalLogNoFile:
          context = Context::Global;
          spec.log_to_file = false;
          continue;
        case Option::GlobalSaveConfig:
          context = Context::Global;
          save_config = true;
          continue;
        case Option::GlobalLoadConfig:
          // load config is loaded in before any processing begins, so this
          // so this is left empty
          continue;
        case Option::GlobalShareInitState:
          context = Context::Global;
          spec.share_init_state = true;
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
        case Option::PhenotypeSpec: {
          context = Context::Phenotype;
          amsim::Phenotype pheno_new{.name = optarg};
          spec.phenotypes.push_back(pheno_new);
          continue;
        }
        case Option::PhenotypeVarGenetic:
          check_context(Context::Phenotype, "--var-genetic");
          spec.phenotypes.back().var_genetic = std::stod(optarg);
          continue;
        case Option::PhenotypeVarEnvironmental:
          check_context(Context::Phenotype, "--var-environmental");
          spec.phenotypes.back().var_environmental = std::stod(optarg);
          continue;
        case Option::PhenotypeVarVertical:
          check_context(Context::Phenotype, "--var-vertical");
          spec.phenotypes.back().var_vertical = std::stod(optarg);
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
        case Option::MatingSpec:
          if (std::string_view(optarg) == "random") spec.mating.type = "random";
          if (std::string_view(optarg) == "assortative") {
            context = Context::Mating;
            spec.mating.type = "assortative";
          }
          continue;
        case Option::MatingCor:
          check_context(Context::Mating, "--mate-cor");
          context = Context::MatingCorrelation;
          continue;
        case Option::MatingErrorTolerance:
          check_context(Context::Mating, "--tol-inf");
          spec.mating.tolerance = std::stod(optarg);
          continue;
        case Option::MatingMaxIterations:
          check_context(Context::Mating, "--max-itr");
          spec.mating.max_iterations = std::stoull(optarg);
          continue;
        case Option::MatingAnnealingTempInit:
          check_context(Context::Mating, "--temp-init");
          spec.mating.initial_temperature = std::stod(optarg);
          continue;
        case Option::MatingAnnealingTempDecay:
          check_context(Context::Mating, "--temp-decay");
          spec.mating.temperature_decay = std::stod(optarg);
          continue;
      }

      // Estimators and sampling
      switch (opt) {
        case Option::PopulationEstimatorSpec: {
          context = Context::PopulationEstimator;

          auto [name, params] = amsim::parse_function(optarg);
          auto estimator = amsim::build_population_estimator(name, params);
          estimator.name = optarg;
          spec.estimators.push_back(estimator);

          continue;
        }
        case Option::SampleEstimatorSpec: {
          context = Context::SampleEstimator;
          amsim::SampleEstimatorSpec sample_est_desc{.name = optarg};
          spec.sample_estimator_spec.push_back(sample_est_desc);
          continue;
        }
        case Option::SampleEstimatorType: {
          check_context(Context::SampleEstimator, "--type");
          auto [name, params] = amsim::parse_function(optarg);
          spec.sample_estimator_spec.back().type = name;
          if (!params.empty())
            spec.sample_estimator_spec.back().params = params;
          continue;
        }
        case Option::SampleEstimatorExec:
          check_context(Context::SampleEstimator, "--exec");
          spec.sample_estimator_spec.back().exec = optarg;
          continue;
        case Option::SampleEstimatorNumRows:
          check_context(Context::SampleEstimator, "--n-rows");
          spec.sample_estimator_spec.back().n_rows = std::stoull(optarg);
          continue;
        case Option::SampleEstimatorNumCols:
          check_context(Context::SampleEstimator, "--n-cols");
          spec.sample_estimator_spec.back().n_cols = std::stoull(optarg);
          continue;
        case Option::SampleEstimatorRowNames:
          check_context(Context::SampleEstimator, "--row-names");
          spec.sample_estimator_spec.back().row_names =
              amsim::utils::split_string(optarg);
          continue;
        case Option::SampleEstimatorColNames:
          check_context(Context::SampleEstimator, "--col-names");
          spec.sample_estimator_spec.back().col_names =
              amsim::utils::split_string(optarg);
          continue;
        case Option::SampleSpec: {
          context = Context::Sample;
          amsim::SampleSpec sample_desc;
          sample_desc.name = optarg;
          spec.sample_spec.push_back(sample_desc);
          continue;
        }
        case Option::SampleProbandType:
          spec.sample_spec.back().proband_type = optarg;
          continue;
        case Option::SampleNumProbands:
          check_context(Context::Sample, "--n-probands");
          spec.sample_spec.back().n_probands = std::stoull(optarg);
          continue;
        case Option::SampleWeightOnPhenotypes: {
          check_context(Context::Sample, "--on");
          spec.sample_spec.back().on = amsim::utils::split_string(optarg);
          continue;
        }
        case Option::SampleWeightOfMembers:
          check_context(Context::Sample, "--of");
          spec.sample_spec.back().of = amsim::utils::split_string(optarg);
          continue;
        case Option::SampleWeightAggregation:
          check_context(Context::Sample, "--agg");
          spec.sample_spec.back().agg = optarg;
          continue;
        case Option::SampleWeightFunction:
          check_context(Context::Sample, "--weight");
          spec.sample_spec.back().weight_function = optarg;
          continue;
        case Option::SampleEstimators:
          check_context(Context::Sample, "--estimators");
          spec.sample_spec.back().estimators =
              amsim::utils::split_string(optarg);
          continue;
      }

      // Context-sensitive (virtual) operation flags
      Context domain = context_domain(context);
      switch (opt) {
        case Option::VirtualNumLoci: {
          std::size_t n_loc = std::stoull(optarg);
          if (domain == Context::Phenotype)
            spec.phenotypes.back().n_causal_loci = n_loc;
          else
            spec.genome.n_loci = n_loc;
          continue;
        }
        case Option::VirtualValue: {
          if (context == PhenotypeCausalLoci) {
            std::vector<std::size_t> value = amsim::parse_indices(optarg);
            spec.phenotypes.back().causal_loci = value;
            continue;
          }

          Eigen::MatrixXd matrix = amsim::parse_matrix_value(optarg);
          if (context == Context::GenomeInitMAFs) spec.genome.v_maf = matrix(0);
          if (context == Context::GenomeRecombinationProbs)
            spec.genome.v_rec = matrix(0);
          if (context == Context::GenomeMutationProbs)
            spec.genome.v_mut = matrix(0);
          if (context == Context::PhenotypeGeneticCorrelation)
            spec.genetic_component_cor = matrix;
          if (context == Context::PhenotypeEnvironmentalCorrelation)
            spec.environmental_component_cor = matrix;
          if (context == Context::MatingCorrelation)
            spec.mating.mate_cor = matrix;
          continue;
        }
        case Option::VirtualFile: {
          if (context == PhenotypeCausalLoci) {
            auto file = amsim::File<std::vector<std::size_t>>{.path = optarg};
            spec.phenotypes.back().causal_loci = file;
            continue;
          }

          auto file = amsim::File<Eigen::MatrixXd>{.path = optarg};
          if (context == Context::GenomeInitMAFs) spec.genome.v_maf = file;
          if (context == Context::GenomeRecombinationProbs)
            spec.genome.v_rec = file;
          if (context == Context::GenomeMutationProbs) spec.genome.v_mut = file;
          if (context == Context::PhenotypeGeneticCorrelation)
            spec.genetic_component_cor = file;
          if (context == Context::PhenotypeEnvironmentalCorrelation)
            spec.environmental_component_cor = file;
          if (context == Context::MatingCorrelation)
            spec.mating.mate_cor = file;
          continue;
        }
        case Option::VirtualSingularValues: {
          std::vector<double> values = amsim::parse_doubles(optarg);
          if (context == Context::PhenotypeGeneticCorrelation) {
            spec.genetic_component_cor =
                amsim::utils::matrix_from_singular_values(values, true);
          }
          if (context == Context::PhenotypeEnvironmentalCorrelation) {
            spec.environmental_component_cor =
                amsim::utils::matrix_from_singular_values(values, true);
          }
          if (context == Context::MatingCorrelation) {
            spec.mating.mate_cor =
                amsim::utils::matrix_from_singular_values(values);
          }
          continue;
        }
        case Option::VirtualDistribution: {
          amsim::Distribution dist = amsim::parse_distribution(
              optarg, (context & Context::GenomeProbabilities) != 0U);

          if (context == Context::GenomeInitMAFs) spec.genome.v_maf = dist;
          if (context == Context::GenomeRecombinationProbs)
            spec.genome.v_rec = dist;
          if (context == Context::GenomeMutationProbs) spec.genome.v_mut = dist;
          if (context_domain(context) == Context::Phenotype)
            spec.phenotypes.back().effects = dist;
          continue;
        }
      }
    }
    // NOLINTEND(bugprone-switch-missing-default-case)
    if (save_config) amsim::ConfigWriter config(spec);
    if (run) amsim::run_simulation(spec);
    exit(EXIT_SUCCESS);
  } catch (const std::exception& e) {
    amsim::Log::error(e.what());
    exit(EXIT_FAILURE);
  }
}
