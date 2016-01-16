//
// Created by pavel on 11/9/15.
//

#ifndef MGRA_TREE_BUILDER_HPP
#define MGRA_TREE_BUILDER_HPP


namespace structure {

namespace phyl_tree {

template<class tree_t>
struct TreeBuilder {
    using node_t = typename tree_t::node_t;
    using node_ptr = typename tree_t::node_ptr;

    //using genome_number_t = std::unordered_map<std::string, size_t>;
    //using mcolor_name_t = std::map<mcolor_t, std::string>;
    //using branch_t = typename structure::Branch<mcolor_t>;

    static tree_t build_tree(std::string str) {
        assert(!str.empty()); // No string - no tree
        tree_t tree(recursive_build_internal_node(str));
        return tree;
    }

    /*static tree_t build_tree(std::vector<branch_t> &branches, size_t value_of_all_genomes) {
        // No branches - no tree
        assert(!branches.empty());

        // Sort is done, so singleton nodes won't break nodes, for example
        // Imagine we have a = ABC|D, b = AB|CD, c = B|ACD: if we apply a -> c -> b we won't get node AB
        std::sort(std::begin(branches), std::end(branches), [](branch_t const &left, branch_t const &right) {
            auto left_diff = std::abs(
                    static_cast<long>(left.left.size()) - static_cast<long>(left.right.size()));
            auto right_diff = std::abs(
                    static_cast<long>(right.left.size()) - static_cast<long>(right.right.size()));
            return left_diff < right_diff;
        });

        auto branch_iter = branches.begin();
        auto root_node = std::make_shared<node_t>(branch_iter->convert_to_node());

        for (++branch_iter; branch_iter != branches.end(); ++branch_iter) {
            merge_branch_into_node(root_node, *branch_iter);
        }

        tree_t tree(root_node);
        return tree;
    }*/

private:
    /**
     *
     * @return
     */
    static node_ptr recursive_build_internal_node(std::string str) {
        boost::trim(str);
        assert(!str.empty()); // No string - no node

        size_t ind_start_edge = truncate_edge_name(str);
        std::string edge_name = str.substr(ind_start_edge);
        boost::trim(edge_name);
        std::string new_tree = str.substr(0, ind_start_edge - (ind_start_edge != str.length()));
        boost::trim(new_tree);

        size_t ind_start_node = truncate_node_name(new_tree);
        std::string node_name = new_tree.substr(ind_start_node);
        boost::trim(node_name);
        new_tree = str.substr(0, ind_start_node);
        boost::trim(new_tree);

        assert(ind_start_edge >= ind_start_node);

        if (new_tree[0] == '(' && new_tree[new_tree.length() - 1] == ')') {
            std::vector<size_t> indeces = get_indices_of_childrens(new_tree);
            size_t start_position = *indeces.begin();
            std::vector<node_ptr> childs;
            for (auto ind = (++indeces.begin()); ind != indeces.end(); ++ind) {
                childs.push_back(recursive_build_internal_node(new_tree.substr(start_position, *ind - start_position)));
                start_position = *ind + 1;
            }
            return node_ptr(new node_t(node_name, childs));;
        } else if (new_tree[0] != '(' && new_tree[new_tree.length() - 1] != ')') {
            return node_ptr(new node_t(str));
        } else {
            std::cerr << "Bad format input (sub)tree. Check count \'(\' and \')\'" << std::endl;
            exit(1);
        }
    }

    /**
     *
     * @return
     */
    static std::vector<size_t> get_indices_of_childrens(std::string const & str) {
        std::vector<size_t> indices({1});

        int p = 0;
        for (size_t j = 1; j < str.length() - 1; ++j) {
            if (str[j] == '(') {
                ++p;
            } else if (str[j] == ')') {
                --p;
            } else if (str[j] == ',') {
                if (p == 0) {
                    indices.push_back(j);
                }
            }
        }

        if (p != 0 || indices.size() < 1) {
            std::cerr << "Bad format input tree. Check count \'(\' and \')\'" << std::endl;
            exit(1);
        } else {
            indices.push_back(str.length() - 1);
        }

        return indices;
    }
    /**
     *
     * @return
     */
    static size_t truncate_node_name(std::string const &str) {
        int i = 0;

        for (i = (int32_t) (str.length() - 1); str[i] != ')' && str[i] != '}' && i >= 0; --i);

        if (i < 0) {
            return str.length();
        } else {
            return (size_t) (i + 1);
        }
    }

    /**
     *
     * @return
     */
    static size_t truncate_edge_name(std::string const &str) {
        int i = 0;

        for (i = (int32_t) (str.length() - 1); str[i] != ':' && str[i] != ')' && str[i] != '}' && i >= 0; --i);

        if (i < 0 || str[i] == ')' || str[i] == '}') {
            return str.length();
        } else {
            return (size_t) (i + 1);
        }
    }

    /**
     * Merges a given branch into tree
     * @return true if merged correctly, false if found an intersection
     */
    /*
    static bool merge_branch_into_node(node_ptr node, branch_t const &branch) {
        node_ptr current_node = node;

        while (true) {
            if (current_node->is_leaf()) {
                if (current_node->get_data() == branch.left) {
                    // Analogous to completeness check
                    return false;
                }
                branch.fill_node(current_node);
                return true;
            }

            if (branch.node_left_includes(current_node)) {
                current_node = current_node->get_most_left_child();
                continue;
            }

            if (branch.node_right_includes(current_node)) {
                current_node = current_node->get_most_right_child();
                continue;
            }

            // Neither includes left or right, found intersection
            return false;
        }
    }*/

   /**
    * Removes the unnecessary children from the node
    * @return true if the node needs to be pruned itself
    */
    /*static bool does_need_to_prune_node(node_ptr const &node, std::unordered_set<mcolor_t> const &appearing_colors) {
        // Node needs to be pruned in two cases:
        // 1. It doesn't appear in statistics set
        // 2. All its children have to be pruned

        if (node->is_basic_leaf()) {
            return (appearing_colors.count(node->get_data()) == 0);
        }

        if (node->is_subtree_leaf()) {
            return true;
        }

        bool both_children_appear_in_statistics =
                appearing_colors.count(node->get_most_left_child()->get_data()) != 0 &&
                appearing_colors.count(node->get_most_right_child()->get_data()) != 0;
        bool needs_to_be_pruned =
                !both_children_appear_in_statistics &&
                does_need_to_prune_node(node->get_most_left_child(), appearing_colors) &&
                does_need_to_prune_node(node->get_most_right_child(), appearing_colors);

        return needs_to_be_pruned;
    }*/
};

}

}

#endif //MGRA_TREE_BUILDER_HPP
