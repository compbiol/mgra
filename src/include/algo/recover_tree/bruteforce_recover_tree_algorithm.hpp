#ifndef BRUTEFORCE_RECOVER_TREE_ALGORITHM_HPP__
#define BRUTEFORCE_RECOVER_TREE_ALGORITHM_HPP__

#include "recover_tree_algorithm.hpp"
#include "graph/graph_pack.hpp"
#include "greedy_tree_builder.hpp"
#include "structures/print_node_visitor.hpp"

namespace algo {

  template <class graph_pack_t>
  struct BruteforceRecoverTreeAlgorithm : RecoverTreeAlgorithm<graph_pack_t> {
    using mcolor_t = typename RecoverTreeAlgorithm<graph_pack_t>::mcolor_t;
    using tree_t = typename RecoverTreeAlgorithm<graph_pack_t>::tree_t;
    using node_t = typename tree_t::colored_node_t;
    using tree_ptr = typename RecoverTreeAlgorithm<graph_pack_t>::tree_ptr;
    using branch_t = std::pair<mcolor_t, mcolor_t>;
    using statistic_t = std::pair<branch_t, size_t>;

    BruteforceRecoverTreeAlgorithm(graph_pack_t& graph_pack) : m_graph_pack(graph_pack) {}

    tree_ptr recover_tree() {
      std::map<branch_t, size_t> edges_statistics = m_graph_pack.stats.multiedges_count;
      std::vector<statistic_t> color_edges_pairs(edges_statistics.begin(), edges_statistics.end());

      std::sort(color_edges_pairs.begin(), color_edges_pairs.end(),
          // Descending by number of edges sort
          [](statistic_t left, statistic_t right) {
            return left.second > right.second;
          });

      GreedyTreeBuilder<tree_t> builder(color_edges_pairs);
      structure::PrintNodeVisitor<node_t> visitor(std::cout);
      visitor.visit(builder.get_result());

      return tree_ptr(new tree_t(builder.get_result()));
    }

  private:
    graph_pack_t& m_graph_pack;
  };
}

#endif