#ifndef TXT_NEWICK_TREE_HPP_
#define TXT_NEWICK_TREE_HPP_

namespace writer {

template <class tree_t>
struct TXT_NewickTree {
  using tree_ptr = typename tree_t::tree_ptr;
  using node_ptr = typename tree_t::node_ptr;
  using tree_vector = std::vector<tree_ptr>;

  TXT_NewickTree(std::ostream& out) 
  : m_out(out) 
  {
  }

  void print_tree(tree_ptr const& tree) {
    print_node(tree->get_root());
    end_tree();
    newline();
  }

  void print_trees(tree_vector const& trees) {
    for (auto& tree: trees) {
      print_tree(tree);
    }
  }

  void print_node(node_ptr node) {
    if (node->is_leaf()) {
      m_out << node->get_name();
    } else {
      start_node();
      print_node(node->get_left_child());
      comma();
      space();
      print_node(node->get_right_child());
      end_node();
    }
  }

private:
  void start_node() {
    m_out << "(";
  }

  void end_node() {
    m_out << ")";
  }

  void end_tree() {
    m_out << ";";
  }

  void comma() {
    m_out << ",";
  }

  void space() {
    m_out << " ";
  }

  void newline() {
    m_out << "\n";
  }

  std::ostream& m_out;
};

}

#endif