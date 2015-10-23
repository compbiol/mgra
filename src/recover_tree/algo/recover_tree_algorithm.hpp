#ifndef RECOVER_TREE_ALGORITHM_HPP
#define RECOVER_TREE_ALGORITHM_HPP

#include "structures/tree.hpp"

namespace algo {

template<class graph_pack_t>
struct RecoverTreeAlgorithm {
    using algo_ptr = std::shared_ptr<RecoverTreeAlgorithm>;
    using mcolor_t = typename graph_pack_t::mcolor_type;
    using tree_t = structure::BinaryTree<mcolor_t>;
    using tree_ptr = std::shared_ptr<tree_t>;
    using tree_vector = std::vector<tree_ptr>;

    virtual tree_vector recover_trees() = 0;

    tree_vector get_result() {
        auto recovered_trees = recover_trees();
        std::for_each(std::begin(recovered_trees), std::end(recovered_trees), [this](tree_ptr const &tree) {
            return finalize_tree(tree);
        });
        return recovered_trees;
    }

    virtual tree_ptr finalize_tree(tree_ptr tree) const {
        return tree;
    }

    virtual ~RecoverTreeAlgorithm() {
    }
};

}

#endif
