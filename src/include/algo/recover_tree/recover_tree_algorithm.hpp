#ifndef RECOVER_TREE_ALGORITHM_HPP__
#define RECOVER_TREE_ALGORITHM_HPP__

#include <memory>
#include "../../structures/tree.hpp"

namespace algo {

  template<class mcolor_t>
  class RecoverTreeAlgorithm {
  public:
    using algo_ptr = std::shared_ptr<RecoverTreeAlgorithm>;
    using tree_t = structure::BinaryTree<mcolor_t>;
    using tree_ptr = std::shared_ptr<tree_t>;

    virtual tree_ptr recover_tree() = 0;

    virtual ~RecoverTreeAlgorithm() {
    }
  };
}

#endif
