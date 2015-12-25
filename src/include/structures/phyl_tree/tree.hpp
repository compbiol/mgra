#ifndef TREE_HPP
#define TREE_HPP

#include "node.hpp"
#include <structures/branch.hpp>

namespace structure {

namespace phyl_tree {

template<class mcolor_t, template<typename> class node_t = Node>
struct BinaryTree {
    using colored_node_t = node_t<mcolor_t>;
    using tree_t = BinaryTree<mcolor_t, node_t>;
    using tree_ptr = std::shared_ptr<tree_t>;
    using branch_t = std::pair<mcolor_t, mcolor_t>;
    using node_ptr = typename colored_node_t::node_ptr;
    using node_const_ptr = typename colored_node_t::node_const_ptr;

    BinaryTree() = default;

    explicit BinaryTree(node_ptr root_node)
            : root(root_node), phylogentic_root_tree(true) { }

    void set_phylogenetic_root(bool value) {
        phylogentic_root_tree = value;
    }

    node_const_ptr get_root() const {
        return root;
    }

    bool is_phylogenetic_root() const {
        return phylogentic_root_tree;
    }

    friend std::ostream &operator<<(std::ostream &out, tree_ptr const &tree) {
        out << tree->get_root();
        return out;
    }

private:
    node_ptr root;
    bool phylogentic_root_tree;

    DECL_LOGGER("BinaryTree")
};

}

}

#include "structures/phyl_tree/tree_builder.hpp"
#include "structures/phyl_tree/algorithms.hpp"

#endif
