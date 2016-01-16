//
// Created by Nikita Kartashov on 05/04/2015.
//

#ifndef MGRA_MULTIEDGES_COUNT_STATISTICS_PRODUCER_HPP_
#define MGRA_MULTIEDGES_COUNT_STATISTICS_PRODUCER_HPP_

#include "statistics_producer.hpp"

namespace algo {

template<class graph_pack_t>
struct DistributionProducer : StatisticsProducer<graph_pack_t> {
    using mcolor_t = typename StatisticsProducer<graph_pack_t>::mcolor_t;
    using tree_t = typename StatisticsProducer<graph_pack_t>::tree_t;
    using branch_t = typename StatisticsProducer<graph_pack_t>::branch_t;

    DistributionProducer(graph_pack_t const & graph_pack, std::vector<tree_t> const &known_subtrees = cfg::get().phylotrees) :
            StatisticsProducer<graph_pack_t>(graph_pack, known_subtrees) {
    }

private:
    void populate_result() {
        auto & multiedges_count = StatisticsProducer<graph_pack_t>::m_graph_pack.stats.multiedges_count;
        std::copy(multiedges_count.begin(), multiedges_count.end(),
                  std::back_inserter(StatisticsProducer<graph_pack_t>::m_result_statistics));
    }
};

}

#endif //MGRA_MULTIEDGES_COUNT_STATISTICS_PRODUCER_HPP_
