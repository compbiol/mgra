#ifndef BRUTEFORCE_RECOVER_TREE_ALGORITHM_HPP__
#define BRUTEFORCE_RECOVER_TREE_ALGORITHM_HPP__

#include "recover_tree_algorithm.hpp"
#include "../../graph/graph_pack.hpp"

namespace algo {

  template <class graph_pack_t>
  struct BruteforceRecoverTreeAlgorithm : RecoverTreeAlgorithm<graph_pack_t> {
    using mcolor_t = typename RecoverTreeAlgorithm<graph_pack_t>::mcolor_t;

    BruteforceRecoverTreeAlgorithm(graph_pack_t& graph_pack) : m_graph(graph_pack) {
    }

    typename RecoverTreeAlgorithm<graph_pack_t>::tree_ptr recover_tree() {
      return std::make_shared<typename structure::BinaryTree<mcolor_t> >(nullptr);
    }

  private:
    graph_pack_t& m_graph;
  };
}

#endif