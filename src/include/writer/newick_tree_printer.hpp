#ifndef NEWICK_TREE_PRINTER_HPP__
#define NEWICK_TREE_PRINTER_HPP__

#include <iostream>
#include <bits/stl_bvector.h>

namespace writer {

  template <class tree_t>
  struct NewickTreePrinter {
    using tree_ptr = typename tree_t::tree_ptr;
    using node_ptr = typename tree_t::node_ptr;
    using name_vector = std::vector<std::string>;

    NewickTreePrinter(std::ostream& out): m_out(out), m_names(name_vector()) {
    }

    NewickTreePrinter(std::ostream& out, name_vector& names): m_out(out), m_names(names) {
    }

    void print_tree(tree_ptr tree) {
      print_node(tree->get_root());
      end_tree();
      newline();
    }

    void print_node(node_ptr node) {
      if (node->is_leaf()) {
        auto color_iter = node->get_data().cbegin();
        if (node->is_complete()) {
          // Get index of singleton color
          m_out << get_name(color_iter->first);
        } else {
          start_unknown_subtree();
          m_out << get_name(color_iter->first);
          for (; color_iter != node->get_data().cend(); ++color_iter) {
            comma();
            space();
            m_out << get_name(color_iter->first);
          }
          end_unknown_subtree();
        }
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
    std::string get_name(size_t color_index) {
      if (m_names.size() > color_index) {
        return m_names[color_index];
      }
      return std::to_string(color_index);
    }

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

    void start_unknown_subtree() {
      m_out << "{";
    }

    void end_unknown_subtree() {
      m_out << "}";
    }

    std::ostream& m_out;
    name_vector m_names;
  };
}

#endif