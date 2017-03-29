//
// Created by Nikita Kartashov on 17/03/2015.
//

#ifndef MGRA_STATISTICS_PRODUCER_HPP_
#define MGRA_STATISTICS_PRODUCER_HPP_

namespace algo {

template<class graph_pack_t>
struct StatisticsProducer {
    using mcolor_t = typename graph_pack_t::mcolor_type;
    using tree_t = typename tree_config<mcolor_t>::phylogeny_tree_t;
    using branch_t = typename structure::phyl_tree::Branch<mcolor_t>;
    using statistic_t = std::pair<branch_t, size_t>;

    StatisticsProducer(graph_pack_t const &graph_pack, std::vector<tree_t> const &known_subtrees = cfg::get().phylotrees)
        : m_graph_pack(graph_pack)
        , m_known_subtrees(known_subtrees)
    {}

    std::vector<statistic_t> make_statistics() {
        // Overloaded in children, gets the statistic data from graph_pack
        populate_result();

        // No statistics - no trees
        assert(!m_result_statistics.empty());

        // Filter complete colors
        m_result_statistics.erase(std::remove_if(std::begin(m_result_statistics),
                                                 std::end(m_result_statistics),
                                                 [](statistic_t const &statistic) {
                                                     return statistic.first.left.empty() ||
                                                            statistic.first.right.empty();
                                                 }), std::end(m_result_statistics));

        // Disregard irregular edges
        for (auto & statistic: m_result_statistics) {
            auto irregular_edges_count = m_graph_pack.stats.irrer_multiedges_count.find(statistic.first);
            if (irregular_edges_count != m_graph_pack.stats.irrer_multiedges_count.end()) {
                statistic.second -= std::min(irregular_edges_count->second, statistic.second);
            }
        }

        if (!m_known_subtrees.empty()) {
            size_t cumulative_score = 0;

            for (auto const & statistic: m_result_statistics) {
                cumulative_score += statistic.second;
            }

            /*for (auto const & tree : m_known_subtrees) {
                for (auto const & branch: structure::phyl_tree::break_tree_into_branches<mcolor_t>(m_graph_pack.multicolors.get_complete_color(), tree)) {
                    m_result_statistics.push_back(statistic_t(branch, cumulative_score));
                }
            }*/
        }

        // Flip edges, so the smaller color is on the left
        /*for (auto &statistic: m_result_statistics) {
            if (statistic.first.left.size() > statistic.first.right.size()) {
                statistic.first.canonize();
            }
        }*/

        // Descending by number of edges sort
        std::sort(std::begin(m_result_statistics), std::end(m_result_statistics),
                  [](statistic_t left, statistic_t right) {
                      return left.second > right.second;
                  });

        return m_result_statistics;
    }

protected:
    virtual void populate_result() = 0;

protected:
    graph_pack_t const & m_graph_pack;
    std::vector<tree_t> const & m_known_subtrees;
    std::vector<statistic_t> m_result_statistics;
};

}

#endif //_MGRA_STATISTICS_PRODUCER_HPP_
