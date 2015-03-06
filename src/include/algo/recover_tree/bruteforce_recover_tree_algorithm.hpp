#ifndef BRUTEFORCE_RECOVER_TREE_ALGORITHM_HPP__
#define BRUTEFORCE_RECOVER_TREE_ALGORITHM_HPP__

#include "recover_tree_algorithm.hpp"
#include "../../graph/graph_pack.hpp"

namespace algo {

  template<class mcolor_t>
  class BruteforceRecoverTreeAlgorithm : public RecoverTreeAlgorithm<mcolor_t> {
  public:
    using graph_pack_t = GraphPack<mcolor_t>;

    BruteforceRecoverTreeAlgorithm(graph_pack_t &breakpoint_graph) : m_graph(breakpoint_graph) {}

    typename RecoverTreeAlgorithm<mcolor_t>::tree_ptr recover_tree() {
      return std::make_shared<typename structure::BinaryTree<mcolor_t> >(nullptr);
    }

  private:
    graph_pack_t& m_graph;
  };
}

#endif