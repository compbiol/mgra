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
      std::map<branch_t, size_t> intermediate_results;

      for (statistic_t path_statistic: StatisticsProducer<graph_pack_t>::m_graph_pack.stats.simple_paths) {
        intermediate_results[path_statistic.first] = path_statistic.second / 2;
      }

      for (statistic_t cycle_statistic: StatisticsProducer<graph_pack_t>::m_graph_pack.stats.simple_cycles) {
        intermediate_results[cycle_statistic.first] = cycle_statistic.second / 2 - 1;
      }

      auto complete_color = cfg::get().complete_color();

      for (auto const& collection: {StatisticsProducer<graph_pack_t>::m_graph_pack.stats.bag_count,
                                    StatisticsProducer<graph_pack_t>::m_graph_pack.stats.cylinder_count}) {
        std::for_each(std::begin(collection),
                      std::end(collection),
                      [&complete_color, &intermediate_results](std::pair<mcolor_t, size_t> const& multiedge) {
                        auto branch = multiedge.first.packed_compliment(complete_color);
                        intermediate_results[branch] += multiedge.second;
                      });
      }

      std::copy(intermediate_results.begin(),
                intermediate_results.end(),
                std::back_inserter(StatisticsProducer<graph_pack_t>::m_result_statistics));
    }
  };
}

#endif //MGRA_SIMPLE_PATH_STATISTICS_PRODUCER_HPP_
