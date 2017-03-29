#ifndef BRANCH_HPP
#define BRANCH_HPP

namespace structure {

namespace phyl_tree {

template<class mcolor_t>
struct Branch {
    mcolor_t left;
    mcolor_t right;

    Branch(mcolor_t const &l, mcolor_t const &r)
    : left(l)
    , right(r)
    {
        assert(!left.empty());
        assert(!right.empty());
        canonize();
    }

    /**
     * Makes a new node with children colored as the given branch
     */
    structure::phyl_tree::Node<mcolor_t> convert_to_node() const {
        using node_t = typename structure::phyl_tree::Node<mcolor_t>;
        mcolor_t color_center(left, right, mcolor_t::Union);
        node_t node(color_center, std::shared_ptr<node_t>(new node_t(left)),
                    std::shared_ptr<node_t>(new node_t(right)));
        return node;
    }

    /**
     * Checks if two branch multicolors intersect
     */
    bool do_intersect(Branch const &branch) const {
        return !(left.includes(branch.left) || right.includes(branch.left) ||
                 branch.left.includes(left) || branch.right.includes(left));
    }

    /**
     * Make node's children colored as the provided branch
     */
    void fill_node(std::shared_ptr<structure::phyl_tree::Node<mcolor_t>> const &node) const {
        using node_t = typename structure::phyl_tree::Node<mcolor_t>;
        node->add_children(std::make_shared<node_t>(mcolor_t(node->get_data(), left, mcolor_t::Intersection)));
        node->add_children(std::make_shared<node_t>(mcolor_t(node->get_data(), right, mcolor_t::Intersection)));
    }

    /**
    * Suppose we have two branches b1 = Q1|Q2 & b2 = R1|R2
    * Let's say that "b1 left includes b2" if Q1 is a subset of R1 and R2 is a subset of Q2
    * Using the same logic, "b1 right includes b2" if Q2 is a subset of R2 and R1 is a subset of Q1
    * If we have a binary tree node, we can look at it as a branch and use the same logic
    */
    bool node_left_includes(std::shared_ptr<structure::phyl_tree::Node<mcolor_t>> const &node) const {
        return node->get_most_left_child()->get_data().includes(left);
    }

    bool node_right_includes(std::shared_ptr<structure::phyl_tree::Node<mcolor_t>> const &node) const {
        return node->get_most_right_child()->get_data().includes(left);
    }

    bool operator<(Branch const &branch) const {
        return (left < branch.left) || (left == branch.left && right < branch.right);
    }

    bool operator==(Branch const &branch) const {
        return (left == branch.left) && (right == branch.right);
    }

private:
    void canonize() {
        if (left > right) {
            std::swap(left, right);
        }
    }
};

}

}

#endif