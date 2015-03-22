#ifndef RECOVER_TREE_ALGORITHM_HPP_
#define RECOVER_TREE_ALGORITHM_HPP_

#include <memory>
#include <vector>
#include "structures/tree.hpp"

namespace algo {

  template <class graph_pack_t>
  struct RecoverTreeAlgorithm {
    using algo_ptr = std::shared_ptr<RecoverTreeAlgorithm>;
    using mcolor_t = typename graph_pack_t::mcolor_type;
    using tree_t = structure::BinaryTree<mcolor_t>;
    using tree_ptr = std::shared_ptr<tree_t>;
    using tree_vector = std::vector<tree_ptr>;

    virtual tree_vector recover_trees() = 0;

    virtual ~RecoverTreeAlgorithm() {
    }
  };
}

#endif
