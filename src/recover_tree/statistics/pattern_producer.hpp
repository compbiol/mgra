//
// Created by Nikita Kartashov on 05/04/2015.
//

#ifndef MGRA_SIMPLE_PATH_STATISTICS_PRODUCER_HPP_
#define MGRA_SIMPLE_PATH_STATISTICS_PRODUCER_HPP_

#include "statistics_producer.hpp"

namespace algo {
  template <class graph_pack_t>
struct PatternProducer : StatisticsProducer<graph_pack_t> {
    using mcolor_t = typename StatisticsProducer<graph_pack_t>::mcolor_t;
    using tree_t = typename StatisticsProducer<graph_pack_t>::tree_t;
    using branch_t = typename StatisticsProducer<graph_pack_t>::branch_t;
    using statistic_t = typename StatisticsProducer<graph_pack_t>::statistic_t;

    PatternProducer(graph_pack_t const & graph_pack,
                                 std::vector<tree_t> const &known_subtrees = cfg::get().phylotrees) :
            StatisticsProducer<graph_pack_t>(graph_pack, known_subtrees) {
    }

protected:
    void populate_result() {
        std::map<branch_t, size_t> intermediate_results;

        for (auto const & path_statistic: StatisticsProducer<graph_pack_t>::m_graph_pack.stats.simple_paths) {
            intermediate_results[path_statistic.first] += path_statistic.second / 2;
        }

        for (auto const & cycle_statistic: StatisticsProducer<graph_pack_t>::m_graph_pack.stats.simple_cycles) {
            intermediate_results[cycle_statistic.first] += cycle_statistic.second / 2 - 1;
        }

        for (auto const &collection: {StatisticsProducer<graph_pack_t>::m_graph_pack.stats.bag_count,
                                      StatisticsProducer<graph_pack_t>::m_graph_pack.stats.cylinder_count}) {
            std::for_each(std::begin(collection), std::end(collection),
                          [&intermediate_results](std::pair<branch_t, size_t> const & digest) {
                              intermediate_results[digest.first] += digest.second;
                          });
        }

        std::copy(intermediate_results.begin(),
                  intermediate_results.end(),
                  std::back_inserter(StatisticsProducer<graph_pack_t>::m_result_statistics));
    }
};

}

#endif //MGRA_SIMPLE_PATH_STATISTICS_PRODUCER_HPP_
