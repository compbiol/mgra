#ifndef RECOVER_TREE_ALGORITHM_HPP__
#define RECOVER_TREE_ALGORITHM_HPP__

#include <memory>
#include "../../structures/tree.hpp"

namespace algo {

  template <class graph_pack_t>
  struct RecoverTreeAlgorithm {
    using algo_ptr = std::shared_ptr<RecoverTreeAlgorithm>;
    using mcolor_t = typename graph_pack_t::mcolor_type;
    using tree_t = structure::BinaryTree<mcolor_t>;
    using tree_ptr = std::shared_ptr<tree_t>;

    virtual tree_ptr recover_tree() = 0;

    virtual ~RecoverTreeAlgorithm() {
    }
  };
}

#endif
