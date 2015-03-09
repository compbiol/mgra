#ifndef TREE_HPP
#define TREE_HPP

#include "node.hpp"

namespace structure {

  template <class mcolor_t, template <typename> class node_t = Node>
  struct BinaryTree {
    using colored_node_t = node_t<mcolor_t>;
    using node_unique_ptr = typename colored_node_t::node_unique_ptr;


    BinaryTree(node_unique_ptr root_node): root(std::move(root_node)), phylogentic_root_tree(true) {}

    BinaryTree(
        std::string const& st,
        std::unordered_map<std::string, size_t> const& genome_number,
        std::vector<std::string> const& priority_name)
        : root(new colored_node_t(nullptr, st, genome_number, priority_name)), phylogentic_root_tree(false) {
      if (mcolor_t(root->get_left_child()->get_data(), root->get_right_child()->get_data(), mcolor_t::Union).size() ==
          priority_name.size()) {
        phylogentic_root_tree = true;
      }
    }

    node_unique_ptr const& get_root() const {
      return root;
    }

    bool is_phylogenetic_root() const {
      return phylogentic_root_tree;
    }

  private:
    node_unique_ptr root;
    bool phylogentic_root_tree;

    DECL_LOGGER("BinaryTree")
  };

}

#endif
