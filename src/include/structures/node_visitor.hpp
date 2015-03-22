#ifndef NODE_VISITOR_HPP_
#define NODE_VISITOR_HPP_

namespace structure {

  template <class node_t>
  struct NodeVisitor {
    using node_ptr = typename node_t::node_ptr;

    NodeVisitor() {
    }

    virtual void visit(node_ptr const& node) = 0;

    virtual ~NodeVisitor() {
    }

  private:
  };
}

#endif