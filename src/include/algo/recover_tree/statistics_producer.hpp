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
    using branch_t = typename BranchHelper::branch_t;
    using statistic_t = std::pair<branch_t, size_t>;
    using statistic_vector = std::vector<statistic_t>;

    static statistic_vector make_statistics(graph_pack_t& graph_pack) {
      std::map<branch_t, size_t> edges_statistics = graph_pack.stats.simple_paths_count;
      // No statistics - no trees
      assert(!edges_statistics.empty());
      statistic_vector color_edges_pairs;

      // Filter complete colors
      std::copy_if(std::begin(edges_statistics),
          std::end(edges_statistics),
          std::back_inserter(color_edges_pairs),
          [](statistic_t statistic) {
            return !statistic.first.first.empty() && !statistic.first.second.empty();
          });

      // Disregard irregular edges
      for (auto& statistic: edges_statistics) {
        statistic.second -= graph_pack.stats.irrer_multiedges_count[statistic.first.first];
        statistic.second -= graph_pack.stats.irrer_multiedges_count[statistic.first.second];
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
