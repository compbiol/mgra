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
        tree_vector const& known_subtrees = cfg::get().phylotrees) :
        m_graph_pack(graph_pack), m_known_subtrees(known_subtrees) {
    }

    statistic_vector make_statistics() {
      // Overloaded in children, gets the statistic data from graph_pack
      populate_result();

      // No statistics - no trees
      assert(!m_result_statistics.empty());

//      FIXME: fix the way indels are traversed
      for (auto& indel_statistic: m_graph_pack.stats.complement_indel_stats) {
        for (auto& statistic: m_result_statistics) {
          if (statistic.first == indel_statistic.first) {
            statistic.second += indel_statistic.second;
            break;
          }
        }
      }

      // Filter complete colors
      m_result_statistics.erase(std::remove_if(std::begin(m_result_statistics),
          std::end(m_result_statistics),
          [](statistic_t const& statistic) {
            return statistic.first.first.empty() || statistic.first.second.empty();
          }), std::end(m_result_statistics));

      // Disregard irregular edges
      for (auto& statistic: m_result_statistics) {
        statistic.second -= m_graph_pack.stats.irrer_multiedges_count[statistic.first.first];
        statistic.second -= m_graph_pack.stats.irrer_multiedges_count[statistic.first.second];
      }

      size_t cumulative_score = 0;
      std::for_each(std::begin(m_result_statistics), std::end(m_result_statistics),
          [&cumulative_score](statistic_t const& statistic) {
            cumulative_score += statistic.second;
          });

      for (auto& branch: BranchHelper::break_trees_into_branches(m_known_subtrees)) {
        m_result_statistics.push_back(statistic_t(branch, cumulative_score));
      }

      // Flip edges, so the smaller color is on the left
      std::for_each(std::begin(m_result_statistics), std::end(m_result_statistics), [](statistic_t& statistic) {
        auto& branch = statistic.first;
        if (statistic.first.first.size() > statistic.first.second.size()) {
          statistic = statistic_t(BranchHelper::flip(branch), statistic.second);
        }
      });

      std::sort(std::begin(m_result_statistics), std::end(m_result_statistics),
          // Descending by number of edges sort
          [](statistic_t left, statistic_t right) {
            return left.second > right.second;
          });

      if (cfg::get().is_debug) {
        for (auto& statistic: m_result_statistics) {
          std::clog << cfg::get().mcolor_to_name(statistic.first.first) << " + "
              << cfg::get().mcolor_to_name(statistic.first.second) << std::endl;
        }
      }

      return m_result_statistics;
    }

  protected:
    virtual void populate_result() = 0;

    graph_pack_t& m_graph_pack;
    tree_vector const& m_known_subtrees;
    statistic_vector m_result_statistics;
  };
}

#endif //_MGRA_STATISTICS_PRODUCER_HPP_
