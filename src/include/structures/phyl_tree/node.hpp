#ifndef NODE_HPP_
#define NODE_HPP_

namespace structure {

namespace phyl_tree {

template<class mcolor_t>
struct Node {
    using multicolor_t = mcolor_t;
    using node_t = Node<mcolor_t>;
    using node_ptr = std::shared_ptr<node_t>;
    using node_ptr_ref = node_ptr &;
    using node_const_ptr = node_ptr const &;

    enum tag_t { classic, whole_duplication };

    Node() = default;

    explicit Node(mcolor_t const & color)
            : Node(color, classic, nullptr, nullptr)
    { }

    Node(mcolor_t const & color, node_ptr const & lc, node_ptr const & rc)
            : Node(color, classic, lc, rc)
    { }

    Node(mcolor_t const & color, node_ptr const & lc)
        : Node(color, whole_duplication, lc, nullptr)
    { }

    void set_parent(Node *new_parent) {
        parent = new_parent;
    }

    void set_left_child(node_ptr const & node) {
        assert(node_type != whole_duplication);
        left_child = std::move(node);
        left_child->set_parent(this);
    }

    void set_right_child(node_ptr const & node) {
        assert(node_type != whole_duplication);
        right_child = std::move(node);
        right_child->set_parent(this);
    }

    mcolor_t const &get_data() const {
        return data;
    }

    tag_t get_type() const {
        return node_type;
    }

    Node const *get_parent() const {
        return parent;
    }

    node_const_ptr get_left_child() const {
        return left_child;
    }

    node_const_ptr get_right_child() const {
        return right_child;
    }

    node_ptr_ref get_left_child() {
        return left_child;
    }

    node_ptr_ref get_right_child() {
        return right_child;
    }

    bool has_left_child() const {
        return static_cast<bool>(left_child);
    }

    bool has_right_child() const {
        return static_cast<bool>(right_child);
    }

    bool is_leaf() const {
        return !has_left_child() && !has_right_child();
    }

    bool is_basic_leaf() const {
        return !has_left_child() && !has_right_child() && (data.size() <= 1);
    }

    bool is_subtree_leaf() const {
        return !has_left_child() && !has_right_child() && (data.size() > 1);
    }

private:
    Node(mcolor_t const & color, tag_t tag, node_ptr const & lc, node_ptr const & rc)
            : data (color)
            , node_type(tag)
            , parent(nullptr)
            , left_child(lc)
            , right_child(rc)
    {
        if (lc) {
            lc->set_parent(this);
        }

        if (rc) {
            rc->set_parent(this);
        }
    }

private:
    mcolor_t data;
    tag_t node_type;

    /**
    * If the tree downwards has all the leaf nodes with only one color
    */
    Node *parent;
    node_ptr left_child;
    node_ptr right_child;
};

}

}

#endif