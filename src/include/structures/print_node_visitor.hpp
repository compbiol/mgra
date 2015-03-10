#ifndef PRINT_NODE_VISITOR_HPP__
#define PRINT_NODE_VISITOR_HPP__

#include <iostream>

#include "node_visitor.hpp"

namespace structure {

  template <class node_t>
  struct PrintNodeVisitor : NodeVisitor<node_t> {
    using node_ptr = typename NodeVisitor<node_t>::node_ptr;

    PrintNodeVisitor(std::ostream& out) : m_indent(0), m_out(out) {
    }

    void visit(node_ptr const& node) {
      // Check only left, tree is binary
      indented_output();
      m_out << node->get_data();
      m_out << "\n";
      if (node->has_left_child()) {
        increase_indent();
        visit(node->get_left_child());
        visit(node->get_right_child());
        decrease_indent();
      }
    }

    std::ostream& get_stream() {
      return m_out;
    }

    const size_t indent_size = 2;
  private:
    void increase_indent() {
      m_indent += indent_size;
    }

    void decrease_indent() {
      m_indent -= indent_size;
    }

    void indented_output() {
      for (size_t i = 0; i < indent_size; ++i) {
        m_out << "\t";
      }
    }

    size_t m_indent;
    std::ostream& m_out;
  };

}

#endif