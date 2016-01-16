#ifndef NODE_HPP_
#define NODE_HPP_

namespace structure {

namespace phyl_tree {

template<class info_t>
struct Node {
    using data_t = info_t;
    using node_t = Node<data_t>;
    using node_ptr = std::shared_ptr<node_t>;
    using node_ptr_ref = node_ptr &;
    using node_const_ptr = node_ptr const &;

    using iter = typename std::vector<node_ptr>::iterator;
    using citer = typename std::vector<node_ptr>::const_iterator;

    enum tag_t { classic, subtree, whole_duplication };

    Node() = default; //TODO think about it

    /*
     * Constructor for leaf;
     */
    explicit Node(data_t const & data)
    : Node(data, classic, std::vector<node_ptr>())
    { }

    /*
     * Constructor for whole duplication node
     */
    Node(data_t const & data, node_ptr lc)
    : Node(data, whole_duplication, std::vector<node_ptr>({lc}))
    {
        assert(static_cast<bool>(lc));
    }

    /*
     * Constructor for classic binary ancestor node
     */
    Node(data_t const & data, node_ptr lc, node_ptr rc)
    : Node(data, classic, std::vector<node_ptr>({lc, rc}))
    {
        assert(static_cast<bool>(lc));
        assert(static_cast<bool>(rc));
    }

    Node(data_t const & data, std::vector<node_ptr> childs)
    : Node(data, subtree, childs)
    {
        assert(childs.size() > 1);
    }

    void add_children(node_ptr node) {
        assert(node_type != whole_duplication);
        childrens.push_back(node);
        (*childrens.rbegin())->set_parent(this);
        node_type = (childrens.size() > 2) ? subtree : classic;
    }

    void replace_children(node_ptr node, size_t ind) {
        assert(node_type != whole_duplication);
        assert(ind < childrens.size());
        childrens[ind] = node;
        childrens[ind]->set_parent(this);
    }

    node_ptr_ref get_most_left_child() {
        assert(!childrens.empty());
        return *childrens.begin();
    }

    node_ptr_ref get_most_right_child() {
        assert(!childrens.empty());
        return *childrens.rbegin();
    }

    data_t const &get_data() const {
        return m_data;
    }

    tag_t get_type() const {
        return node_type;
    }

    Node const *get_parent() const {
        return parent;
    }

    node_const_ptr get_most_left_child() const {
        assert(!childrens.empty());
        return *childrens.cbegin();
    }

    node_const_ptr get_most_right_child() const {
        assert(!childrens.empty());
        return *childrens.crbegin();
    }

    bool has_childrens() const {
        return !childrens.empty();
    }

    bool is_leaf() const {
        return childrens.empty();
    }

    DECLARE_ITERATOR(iter, childrens, begin, begin)

    DECLARE_ITERATOR(iter, childrens, end, end)

    DECLARE_CONST_ITERATOR(citer, childrens, begin, cbegin)

    DECLARE_CONST_ITERATOR(citer, childrens, end, cend)

    DECLARE_CONST_ITERATOR(citer, childrens, cbegin, cbegin)

    DECLARE_CONST_ITERATOR(citer, childrens, cend, cend)

private:
    Node(data_t const & data, tag_t type, std::vector<node_ptr> const & childs)
    : m_data(data)
    , node_type(type)
    , parent(nullptr)
    , childrens(childs)
    {
        for (auto children : childrens) {
            children->set_parent(this);
        }

        if (childrens.size() == 1) {
            node_type = whole_duplication;
        } else if (childrens.size() == 2 || childrens.empty()) {
            node_type = classic;
        }
    }

    void set_parent(Node *new_parent) {
        parent = new_parent;
    }
private:
    data_t m_data;
    tag_t node_type;

    /**
    * If the tree downwards has all the leaf nodes with only one color
    */
    Node *parent;
    std::vector<node_ptr> childrens;
};

}

}

#endif