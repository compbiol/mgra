#ifndef TREE_HPP
#define TREE_HPP

#include "structures/phyl_tree/node.hpp"

namespace structure {

namespace phyl_tree {

template<class info_t, template<typename> class data_node_t = Node>
struct BinaryTree {
    using data_t = info_t;

    using node_t = data_node_t<info_t>;
    using node_ptr = typename node_t::node_ptr;
    using node_const_ptr = typename node_t::node_const_ptr;

    using tree_t = BinaryTree<info_t, data_node_t>;
    using tree_ptr = std::shared_ptr<tree_t>;

    BinaryTree()
    : root(nullptr)
    {}

    explicit BinaryTree(node_ptr root_node)
    : root(root_node)
    { }

    node_const_ptr get_root() const {
        return root;
    }

private:
    node_ptr root;

    DECL_LOGGER("PhylogeneticTree")
};

}

}

#include "structures/phyl_tree/tree_builder.hpp"
#include "structures/phyl_tree/algorithms.hpp"
#include "structures/phyl_tree/newick_tree_printer.hpp"
#include "structures/phyl_tree/branch.hpp"

#endif
