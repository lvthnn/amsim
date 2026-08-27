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

#include <amsim/core/generation.h>
#include <amsim/core/params.h>
#include <amsim/data/mating.h>

#include <algorithm>
#include <cstddef>
#include <deque>

namespace amsim {

struct PedigreeNode {
  std::size_t index;
  std::size_t depth;
  Sex sex;  // this is immediate from index but convenient
};

struct PedigreePath {
  std::uint64_t seq;
  std::size_t seq_len;
  std::optional<std::size_t> pivot;
};

class Pedigree {
 public:
  explicit Pedigree(const Params& params)
      : max_depth_(params.sim.pedigree_max_depth),
        n_ind_(params.sim.n_ind),
        n_sex_(params.sim.n_ind / 2) {}

  void push(const Matching& matching, const Matching& inv_matching);

  std::size_t depth() const;

  std::size_t max_depth() const;

  std::vector<std::vector<PedigreeNode>> find_cousins(
      const PedigreeNode& self, std::size_t degree = 1) const;

  std::vector<std::vector<PedigreeNode>> find_cousins(
      std::size_t degree = 1) const;

  std::vector<std::vector<PedigreeNode>> find_ancestors(
      const PedigreeNode& self, std::size_t degree = 1) const;

  std::vector<std::vector<PedigreeNode>> find_ancestors(
      std::size_t degree = 1) const;

 private:
  std::size_t max_depth_;
  std::size_t n_ind_;
  std::size_t n_sex_;
  std::deque<std::pair<Matching, Matching>> mate_history_;

  PedigreeNode node(std::size_t index, std::size_t depth) const;

  PedigreeNode father(const PedigreeNode& self) const;
  PedigreeNode mother(const PedigreeNode& self) const;
  PedigreeNode sibling(const PedigreeNode& self) const;
  PedigreeNode son(const PedigreeNode& self) const;
  PedigreeNode daughter(const PedigreeNode& self) const;
  PedigreeNode parent(const PedigreeNode& self, Sex sex) const;
  PedigreeNode offspring(const PedigreeNode& self, Sex sex) const;

  std::vector<PedigreeNode> traversal(
      const PedigreeNode& init, const PedigreePath& path) const;
};

inline void Pedigree::push(
    const Matching& matching, const Matching& inv_matching) {
  if (mate_history_.size() == max_depth_) mate_history_.pop_back();

  mate_history_.emplace_front(matching, inv_matching);
}

inline std::size_t Pedigree::depth() const {
  return mate_history_.size();
}

inline std::size_t Pedigree::max_depth() const {
  return max_depth_;
}

inline std::vector<std::vector<PedigreeNode>> Pedigree::find_cousins(
    const PedigreeNode& self, std::size_t degree) const {
  if (degree == 0)
    throw std::invalid_argument(
        "Pedigree::find_cousins: degree must be non-zero.");

  std::size_t seq_len = 2 * degree;
  std::vector<std::vector<PedigreeNode>> cousins;

  // 2^(2d) bits in total
  for (std::uint64_t seq = 0; seq < (1ULL << seq_len); ++seq) {
    PedigreePath path{.seq = seq, .seq_len = seq_len, .pivot = degree};
    cousins.push_back(traversal(self, path));
  }

  return cousins;
}

inline std::vector<std::vector<PedigreeNode>> Pedigree::find_cousins(
    std::size_t degree) const {
  if (degree == 0)
    throw std::invalid_argument(
        "Pedigree::find_cousins: degree must be non-zero.");

  auto cousins_per_ind = static_cast<std::size_t>(1ULL << (2 * degree));
  std::vector<std::vector<PedigreeNode>> cousins(n_ind_ * cousins_per_ind);

  for (std::size_t ind = 0; ind < n_ind_; ++ind) {
    auto ind_cousins = find_cousins(node(ind, 0), degree);
    std::ranges::move(ind_cousins, (cousins.begin() + (ind * cousins_per_ind)));
  }

  return cousins;
}

inline std::vector<std::vector<PedigreeNode>> Pedigree::find_ancestors(
    const PedigreeNode& self, std::size_t degree) const {
  if (degree == 0)
    throw std::invalid_argument(
        "Pedigree::find_ancestors: degree must be non-zero");

  std::size_t seq_len = degree;
  std::vector<std::vector<PedigreeNode>> ancestors;

  for (std::uint64_t seq = 0; seq < (1ULL << seq_len); ++seq) {
    PedigreePath path{.seq = seq, .seq_len = seq_len};
    ancestors.push_back(traversal(self, path));
  }

  return ancestors;
}

inline std::vector<std::vector<PedigreeNode>> Pedigree::find_ancestors(
    std::size_t degree) const {
  if (degree == 0)
    throw std::invalid_argument(
        "Pedigree::find_ancestors: degree must be non-zero");

  auto ancestors_per_ind = static_cast<std::size_t>(1ULL << degree);
  std::vector<std::vector<PedigreeNode>> ancestors(n_ind_ * ancestors_per_ind);

  for (std::size_t ind = 0; ind < n_ind_; ++ind) {
    auto ind_ancestors = find_ancestors(node(ind, 0), degree);
    std::ranges::move(
        ind_ancestors, ancestors.begin() + (ancestors_per_ind * ind));
  }

  return ancestors;
}

inline PedigreeNode Pedigree::node(std::size_t index, std::size_t depth) const {
  if (depth >= max_depth_)
    throw std::invalid_argument(
        "Pedigree::node depth argument exceeds maximum pedigree depth");

  return PedigreeNode{
      .index = index,
      .depth = depth,
      .sex = (index < n_sex_) ? Sex::Male : Sex::Female};
}

inline PedigreeNode Pedigree::father(const PedigreeNode& self) const {
  if (self.depth == max_depth_)
    throw std::range_error(
        "Pedigree::father invoked at maximum pedigree depth");

  return PedigreeNode{
      .index = (self.sex == Sex::Male)
                   ? self.index
                   : mate_history_[self.depth + 1].second[self.index - n_sex_],
      .depth = self.depth + 1,
      .sex = Sex::Male};
}

inline PedigreeNode Pedigree::mother(const PedigreeNode& self) const {
  if (self.depth == max_depth_)
    throw std::range_error(
        "Pedigree::mother invoked at maximum pedigree depth");

  return PedigreeNode{
      .index = (self.sex == Sex::Male)
                   ? n_sex_ + mate_history_[self.depth + 1].first[self.index]
                   : self.index,
      .depth = self.depth + 1,
      .sex = Sex::Female};
}

inline PedigreeNode Pedigree::parent(const PedigreeNode& self, Sex sex) const {
  if (sex == Sex::Male) return father(self);
  if (sex == Sex::Female) return mother(self);
  throw std::invalid_argument(
      "Pedigree::parent invoked with parental sex Sex::Unknown");
}

inline PedigreeNode Pedigree::sibling(const PedigreeNode& self) const {
  if (self.depth == max_depth_)
    throw std::range_error(
        "Pedigree::sibling invoked at maximum pedigree depth");

  return PedigreeNode{
      .index = (self.sex == Sex::Male)
                   ? n_sex_ + mate_history_[self.depth + 1].first[self.index]
                   : mate_history_[self.depth + 1].second[self.index - n_sex_],
      .depth = self.depth,
      .sex = (self.sex == Sex::Male) ? Sex::Female : Sex::Male};
}

inline PedigreeNode Pedigree::son(const PedigreeNode& self) const {
  if (self.depth == 0)
    throw std::range_error(
        "Pedigree::son invoked at current generation of pedigree");

  return PedigreeNode{
      .index = (self.sex == Sex::Male)
                   ? self.index
                   : mate_history_[self.depth].second[self.index - n_sex_],
      .depth = self.depth - 1,
      .sex = Sex::Male};
}

inline PedigreeNode Pedigree::daughter(const PedigreeNode& self) const {
  if (self.depth == 0)
    throw std::range_error(
        "Pedigree::son invoked at current generation of pedigree");

  return PedigreeNode{
      .index = (self.sex == Sex::Male)
                   ? n_sex_ + mate_history_[self.depth].first[self.index]
                   : self.index,
      .depth = self.depth - 1,
      .sex = Sex::Female};
}

inline PedigreeNode Pedigree::offspring(
    const PedigreeNode& self, Sex sex) const {
  if (sex == Sex::Male) return son(self);
  if (sex == Sex::Female) return daughter(self);
  throw std::invalid_argument("Pedigree::offspring invoked with Sex::Unknown");
}

inline std::vector<PedigreeNode> Pedigree::traversal(
    const PedigreeNode& init, const PedigreePath& path) const {
  std::vector<PedigreeNode> result(path.seq_len + 1);
  result[0] = init;

  for (std::size_t l = 0; l < path.seq_len; ++l) {
    Sex sex = bit_to_sex((path.seq >> l) & 1ULL);
    bool up = l < path.pivot.value_or(path.seq_len);

    // If a path is pivoted — that is, can be decomposed into an ascent sequence
    // followed by a descent sequence — then the "apex" node is irrelevant in
    // the bit sequence seeing as couples are strictly monogamous. Thus, the
    // node preceding the apex node simply steps over to its sibling which
    // always yields a valid cousin-cousin sequence.
    if (path.pivot.has_value() && (l == path.pivot.value()))
      result[l + 1] = offspring(sibling(result[l]), sex);
    else
      result[l + 1] = up ? parent(result[l], sex) : offspring(result[l], sex);
  }

  return result;
}

}  // namespace amsim
