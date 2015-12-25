//
// Created by Nikita Kartashov on 12/04/2015.
//

#ifndef MGRA_STATISTICS_BASED_RECOVER_TREE_ALGORITHM_HPP
#define MGRA_STATISTICS_BASED_RECOVER_TREE_ALGORITHM_HPP

#include "recover_tree_algorithm.hpp"
#include "statistics/statistics_producer.hpp"

namespace algo {

template<class graph_pack_t>
struct StatisticsBasedRecoverTreeAlgorithm : RecoverTreeAlgorithm<graph_pack_t> {
    using tree_t = typename RecoverTreeAlgorithm<graph_pack_t>::tree_t;
    using tree_ptr = typename RecoverTreeAlgorithm<graph_pack_t>::tree_ptr;
    using node_t = typename tree_t::colored_node_t;
    using node_ptr = std::shared_ptr<node_t>;
    using mcolor_t = typename graph_pack_t::mcolor_type;
    using statistic_producer_ptr = std::shared_ptr<StatisticsProducer<graph_pack_t>> const&;

    StatisticsBasedRecoverTreeAlgorithm(graph_pack_t &graph_pack,
                                        statistic_producer_ptr statistic_producer) : m_graph_pack(graph_pack),
                                                                                     m_statistic_producer(statistic_producer)
    { }

    virtual tree_ptr finalize_tree(tree_ptr tree) const {
        std::queue<node_ptr> node_queue;
        node_queue.push(tree->get_root());
        while (!node_queue.empty()) {
            auto current_node = node_queue.front();
            node_queue.pop();
            if (current_node->is_leaf()) {
                // Check nodes such as {E, F}, which can be broken without any more info
                if (current_node->get_data().size() == 2) {
                    auto broken_color = current_node->get_data().break_into_parts();
                    auto left = broken_color[0];
                    auto right = broken_color[1];
                    current_node->set_left_child(std::make_shared<node_t>(left));
                    current_node->set_right_child(std::make_shared<node_t>(right));
                }
            } else {
                node_queue.push(current_node->get_left_child());
                node_queue.push(current_node->get_right_child());
            }
        }
        return tree;
    }

protected:
    graph_pack_t &m_graph_pack;
    statistic_producer_ptr m_statistic_producer;
};

}

#endif //MGRA_STATISTICS_BASED_RECOVER_TREE_ALGORITHM_HPP
