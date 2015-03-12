#ifndef GREEDY_RECOVER_TREE_ALGORITHM_HPP__
#define GREEDY_RECOVER_TREE_ALGORITHM_HPP__

#include "recover_tree_algorithm.hpp"
#include "graph/graph_pack.hpp"
#include "structures/branch.hpp"

namespace algo {

  template <class graph_pack_t>
  struct GreedyRecoverTreeAlgorithm : RecoverTreeAlgorithm<graph_pack_t> {
    using mcolor_t = typename RecoverTreeAlgorithm<graph_pack_t>::mcolor_t;
    using tree_t = typename RecoverTreeAlgorithm<graph_pack_t>::tree_t;
    using node_t = typename tree_t::colored_node_t;
    using node_ptr = std::shared_ptr<node_t>;
    using node_vector = std::vector<node_ptr>;
    using tree_ptr = typename RecoverTreeAlgorithm<graph_pack_t>::tree_ptr;
    using tree_vector = typename RecoverTreeAlgorithm<graph_pack_t>::tree_vector;
    using branch_t = std::pair<mcolor_t, mcolor_t>;
    using statistic_t = std::pair<branch_t, size_t>;
    using statistic_vector = std::vector<statistic_t>;
    using BranchHelper = structure::Branch<mcolor_t>;

    GreedyRecoverTreeAlgorithm(graph_pack_t& graph_pack) : m_graph_pack(graph_pack) {
    }

    tree_vector recover_trees() {
      std::map<branch_t, size_t> edges_statistics = m_graph_pack.stats.multiedges_count;
      statistic_vector color_edges_pairs(std::begin(edges_statistics), std::end(edges_statistics));

      std::sort(std::begin(color_edges_pairs), std::end(color_edges_pairs),
          // Descending by number of edges sort
          [](statistic_t left, statistic_t right) {
            return left.second > right.second;
          });

      return build_trees(color_edges_pairs);
    }

    std::vector<tree_ptr> build_trees(statistic_vector& branch_statistics) {
      // Vector should be sorted in descending format
      validate_statistics(branch_statistics);

      if (branch_statistics.empty()) {
        ERROR("No branch statistics supplied")
        exit(4);
      }

      // Take the best statistic as a root node
      auto iter = std::begin(branch_statistics);
      auto root_node = BranchHelper::node_from_branch(iter->first);
      iter++;

      for (; iter != std::end(branch_statistics); ++iter) {
        BranchHelper::merge_branch_into_node(root_node, iter->first);
      }
      return {std::make_shared<tree_t>(root_node)};
    }

  private:
    /**
    * Checks if statistics are in descending order
    */
    void validate_statistics(statistic_vector const& branch_statistics) {
      // Maximum size_t value
      size_t previous_value = static_cast<size_t>(-1);
      for (auto const& statistic: branch_statistics) {
        assert(previous_value >= statistic.second &&
            "Works only on sorted statistics");
        previous_value = statistic.second;
      }
    }

    graph_pack_t& m_graph_pack;
  };
}

#endif