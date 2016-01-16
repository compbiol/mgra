#ifndef NEWICK_TREE_PRINTER_HPP_
#define NEWICK_TREE_PRINTER_HPP_

namespace structure {

namespace phyl_tree {

template<class tree_t>
struct NewickTreePrinter {
    using node_t = typename tree_t::node_t;
    using node_ptr = typename tree_t::node_ptr;

    NewickTreePrinter(std::ostream &out)
    : m_out(out) {
    }

    void print_tree(tree_t const &tree) {
        print_node(tree.get_root());
        end_tree();
    }

    void print_node(node_ptr const & node) {
        if (node->is_leaf()) {
            m_out << node->get_data();
        } else if (node->get_type() == node_t::whole_duplication) {
            assert(node->get_most_left_child() == node->get_most_right_child());
            print_node(node->get_most_left_child());
        } else {
            assert(node->has_childrens());

            start_node();
            auto child = node->cbegin();
            print_node(*child);
            for (++child; child != node->cend(); ++child) {
                comma();
                space();
                print_node(*child);
            }
            end_node();
            m_out << node->get_data();
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

    std::ostream &m_out;
};

}

}

#endif