#ifndef BRUTEFORCE_RECOVER_TREE_ALGORITHM_HPP__
#define BRUTEFORCE_RECOVER_TREE_ALGORITHM_HPP__

#include "recover_tree_algorithm.hpp"
#include "../../graph/graph_pack.hpp"

namespace algo {

  template <class graph_pack_t>
  struct BruteforceRecoverTreeAlgorithm : RecoverTreeAlgorithm<graph_pack_t> {
    using mcolor_t = typename RecoverTreeAlgorithm<graph_pack_t>::mcolor_t;
    using tree_t = typename RecoverTreeAlgorithm<graph_pack_t>::tree_t;
    using tree_ptr = typename RecoverTreeAlgorithm<graph_pack_t>::tree_ptr;

    BruteforceRecoverTreeAlgorithm(graph_pack_t& graph_pack) : m_graph(graph_pack) {
    }

    tree_ptr recover_tree() {
      return std::make_shared<tree_t>(std::unique_ptr<typename tree_t::Node>(nullptr));
    }

  private:
    graph_pack_t& m_graph;
  };
}

#endif