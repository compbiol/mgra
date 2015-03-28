//
// Created by Nikita Kartashov on 17/03/2015.
//

#ifndef MGRA_STATISTICS_PRODUCER_HPP_
#define MGRA_STATISTICS_PRODUCER_HPP_

#include "structures/branch.hpp"

namespace algo {

  template <class graph_pack_t>
  struct StatisticsProducer {
    using mcolor_t = typename graph_pack_t::mcolor_type;
    using BranchHelper = structure::Branch<mcolor_t>;
    using tree_vector = typename BranchHelper::tree_vector;
    using branch_t = typename BranchHelper::branch_t;
    using statistic_t = std::pair<branch_t, size_t>;
    using statistic_vector = std::vector<statistic_t>;

    static statistic_vector make_statistics(graph_pack_t& graph_pack,
        tree_vector const& known_subtrees = cfg::get().phylotrees) {
      auto& edges_count = graph_pack.stats.simple_multiedges_count;

      // No statistics - no trees
      assert(!edges_count.empty());

      auto complete_color = cfg::get().complete_color();

      statistic_vector edges_statistics;

      std::transform(std::begin(edges_count), std::end(edges_count), std::back_inserter(edges_statistics),
          [&complete_color](std::pair<mcolor_t, size_t> const& multiedge) {
            auto compliment = mcolor_t(complete_color, multiedge.first, mcolor_t::Difference);
            return statistic_t(branch_t(multiedge.first, compliment), multiedge.second);
          });

      statistic_vector color_edges_pairs;

      // Filter complete colors
      std::copy_if(std::begin(edges_statistics),
          std::end(edges_statistics),
          std::back_inserter(color_edges_pairs),
          [](statistic_t const& statistic) {
            return !statistic.first.first.empty() && !statistic.first.second.empty();
          });

      // Disregard irregular edges
      for (auto& statistic: edges_statistics) {
        statistic.second -= graph_pack.stats.irrer_multiedges_count[statistic.first.first];
        statistic.second -= graph_pack.stats.irrer_multiedges_count[statistic.first.second];
      }

      size_t cumulative_score = 0;
      std::for_each(std::begin(color_edges_pairs), std::end(color_edges_pairs),
          [&cumulative_score](statistic_t const& statistic) {
            cumulative_score += statistic.second;
          });

      for(auto& branch: BranchHelper::break_trees_into_branches(known_subtrees)) {
        color_edges_pairs.push_back(statistic_t(branch, cumulative_score));
      }

      // Flip edges, so the smaller color is on the left
      std::for_each(std::begin(color_edges_pairs), std::end(color_edges_pairs), [](statistic_t& statistic) {
        auto& branch = statistic.first;
        if (statistic.first.first.size() > statistic.first.second.size()) {
          statistic = statistic_t(BranchHelper::flip(branch), statistic.second);
        }
      });

      std::sort(std::begin(color_edges_pairs), std::end(color_edges_pairs),
          // Descending by number of edges sort
          [](statistic_t left, statistic_t right) {
            return left.second > right.second;
          });

      return color_edges_pairs;
    }
  };
}

#endif //_MGRA_STATISTICS_PRODUCER_HPP_
