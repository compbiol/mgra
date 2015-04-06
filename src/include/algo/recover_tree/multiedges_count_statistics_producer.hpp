//
// Created by Nikita Kartashov on 05/04/2015.
//

#ifndef MGRA_MULTIEDGES_COUNT_STATISTICS_PRODUCER_HPP_
#define MGRA_MULTIEDGES_COUNT_STATISTICS_PRODUCER_HPP_

#include "statistics_producer.hpp"

namespace algo {
  template <class graph_pack_t>
  struct MultiEdgesCountStatisticsProducer : StatisticsProducer<graph_pack_t> {
    using mcolor_t = typename StatisticsProducer<graph_pack_t>::mcolor_t;
    using BranchHelper = typename StatisticsProducer<graph_pack_t>::BranchHelper;
    using tree_vector = typename StatisticsProducer<graph_pack_t>::tree_vector;
    using branch_t = typename StatisticsProducer<graph_pack_t>::branch_t;
    using statistic_t = typename StatisticsProducer<graph_pack_t>::statistic_t;
    using statistic_vector = typename StatisticsProducer<graph_pack_t>::statistic_vector;

    MultiEdgesCountStatisticsProducer(graph_pack_t& graph_pack,
        tree_vector const& known_subtrees = cfg::get().phylotrees) :
        StatisticsProducer<graph_pack_t>(graph_pack, known_subtrees) {
    }

  private:
    void populate_result() {
      auto& multiedges_count = StatisticsProducer<graph_pack_t>::m_graph_pack.stats.multiedges_count;
      std::copy(multiedges_count.begin(),
          multiedges_count.end(),
          std::back_inserter(StatisticsProducer<graph_pack_t>::m_result_statistics));
    }
  };
}

#endif //MGRA_MULTIEDGES_COUNT_STATISTICS_PRODUCER_HPP_
