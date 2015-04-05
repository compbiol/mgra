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

    StatisticsProducer(graph_pack_t& graph_pack,
        tree_vector const& known_subtrees = cfg::get().phylotrees):
        m_graph_pack(graph_pack), m_known_subtrees(known_subtrees){
    }

    statistic_vector make_statistics() {
      auto& edges_count = m_graph_pack.stats.simple_multiedges_count;

      // No statistics - no trees
      assert(!edges_count.empty());

      auto complete_color = cfg::get().complete_color();

      statistic_vector edges_statistics;

      std::transform(std::begin(edges_count), std::end(edges_count), std::back_inserter(edges_statistics),
          [&complete_color](std::pair<mcolor_t, size_t> const& multiedge) {
            auto compliment = mcolor_t(complete_color, multiedge.first, mcolor_t::Difference);
            return statistic_t(branch_t(multiedge.first, compliment), multiedge.second);
          });

//      FIXME: fix the way indels are traversed
      for (auto& indel_statistic: m_graph_pack.stats.complement_indel_stats) {
        for (auto& statistic: edges_statistics) {
          if (statistic.first == indel_statistic.first) {
            statistic.second += indel_statistic.second;
            break;
          }
        }
      }

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
        statistic.second -= m_graph_pack.stats.irrer_multiedges_count[statistic.first.first];
        statistic.second -= m_graph_pack.stats.irrer_multiedges_count[statistic.first.second];
      }

      size_t cumulative_score = 0;
      std::for_each(std::begin(color_edges_pairs), std::end(color_edges_pairs),
          [&cumulative_score](statistic_t const& statistic) {
            cumulative_score += statistic.second;
          });

      for(auto& branch: BranchHelper::break_trees_into_branches(m_known_subtrees)) {
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

      if (cfg::get().is_debug) {
        for (auto& statistic: color_edges_pairs) {
          std::clog << cfg::get().mcolor_to_name(statistic.first.first) << " + "
              << cfg::get().mcolor_to_name(statistic.first.second) << std::endl;
        }
      }

      return color_edges_pairs;
    }

  protected:
    graph_pack_t& m_graph_pack;
    tree_vector const& m_known_subtrees;
  };
}

#endif //_MGRA_STATISTICS_PRODUCER_HPP_
