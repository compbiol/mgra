//
// Created by Nikita Kartashov on 05/04/2015.
//

#ifndef MGRA_SIMPLE_PATH_STATISTICS_PRODUCER_HPP_
#define MGRA_SIMPLE_PATH_STATISTICS_PRODUCER_HPP_

#include "statistics_producer.hpp"

namespace algo {
  template <class graph_pack_t>
  struct SimplePathStatisticsProducer : StatisticsProducer<graph_pack_t> {
    using mcolor_t = typename StatisticsProducer<graph_pack_t>::mcolor_t;
    using BranchHelper = typename StatisticsProducer<graph_pack_t>::BranchHelper;
    using tree_vector = typename StatisticsProducer<graph_pack_t>::tree_vector;
    using branch_t = typename StatisticsProducer<graph_pack_t>::branch_t;
    using statistic_t = typename StatisticsProducer<graph_pack_t>::statistic_t;
    using statistic_vector = typename StatisticsProducer<graph_pack_t>::statistic_vector;

    SimplePathStatisticsProducer(graph_pack_t& graph_pack,
                                 tree_vector const& known_subtrees = cfg::get().phylotrees) :
        StatisticsProducer<graph_pack_t>(graph_pack, known_subtrees) {
    }

  protected:
    void populate_result() {
      auto& simple_edges_counts = StatisticsProducer<graph_pack_t>::m_graph_pack.stats.simple_multiedges_count;
      auto complete_color = cfg::get().complete_color();
      StatisticsProducer<graph_pack_t>::m_result_statistics.reserve(simple_edges_counts.size());

      std::transform(std::begin(simple_edges_counts),
                     std::end(simple_edges_counts),
                     std::back_inserter(StatisticsProducer<graph_pack_t>::m_result_statistics),
                     [&complete_color](std::pair<mcolor_t, size_t> const& multiedge) {
                       auto compliment = mcolor_t(complete_color, multiedge.first, mcolor_t::Difference);
                       return statistic_t(branch_t(multiedge.first, compliment), multiedge.second);
                     });
    }
  };
}

#endif //MGRA_SIMPLE_PATH_STATISTICS_PRODUCER_HPP_
