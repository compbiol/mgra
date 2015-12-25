//
// Created by pavel on 11/9/15.
//

#ifndef MGRA_TREE_BUILDER_HPP
#define MGRA_TREE_BUILDER_HPP

namespace structure {

namespace phyl_tree {

template<class tree_t>
struct TreeBuilder {
    using node_ptr = typename tree_t::node_ptr;
    using node_t = typename tree_t::colored_node_t;
    using mcolor_t = typename node_t::multicolor_t;
    using branch_t = typename structure::Branch<mcolor_t>;
    using genome_number_t = std::unordered_map<std::string, size_t>;
    using mcolor_name_t = std::map<mcolor_t, std::string>;

    static tree_t build_tree(std::string const &str, genome_number_t const &genome_number,
                             mcolor_name_t &mcolor_to_name, size_t value_of_all_genomes) {
        // No string - no tree
        assert(!str.empty());
        tree_t tree(recursive_build_internal_node(str, genome_number, mcolor_to_name));
        if (mcolor_t(tree.get_root()->get_left_child()->get_data(),
                     tree.get_root()->get_right_child()->get_data(),
                     mcolor_t::Union).size() == value_of_all_genomes) {
            tree.set_phylogenetic_root(true);
        }
        return tree;
    }

    static tree_t build_tree(std::vector<branch_t> &branches, size_t value_of_all_genomes) {
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

        if (mcolor_t(tree.get_root()->get_left_child()->get_data(),
                     tree.get_root()->get_right_child()->get_data(),
                     mcolor_t::Union).size() == value_of_all_genomes) {
            tree.set_phylogenetic_root(true);
        }

        return tree;
    }

private:
    /**
     *
     * @return
     */
    static node_ptr recursive_build_internal_node(std::string const &str, genome_number_t const &genome_number,
                                                  mcolor_name_t &mcolor_to_name) {
        size_t ind = truncate_name(str);
        std::string new_tree = str.substr(0, ind);
        node_ptr node;

        if (new_tree[0] == '(') {
            //non-trivial tree
            if (new_tree[new_tree.length() - 1] != ')') {
                std::cerr << "Bad format input (sub)tree. Check count \'(\' and \')\'" << std::endl;
                exit(1);
            }

            int p = 0;
            node_ptr left_child;
            node_ptr right_child;

            for (size_t j = 1; j < new_tree.size() - 1; ++j) {
                if (new_tree[j] == '(' || new_tree[j] == '{') {
                    ++p;
                } else if (new_tree[j] == ')' || new_tree[j] == '}') {
                    --p;
                } else if (new_tree[j] == ',') {
                    if (p == 0) {
                        left_child = recursive_build_internal_node(new_tree.substr(1, j - 1), genome_number,
                                                                   mcolor_to_name);
                        right_child = recursive_build_internal_node(new_tree.substr(j + 1, new_tree.length() - j - 2),
                                                                    genome_number, mcolor_to_name);
                    }
                }

                if (p < 0) {
                    std::cerr << "Bad format input (sub)tree. Check count \'(\' and \')\'" << std::endl;
                    exit(1);
                }
            }

            if (p != 0) {
                std::cerr << "Bad format input (sub)tree. Check count \'(\' and \')\'" << std::endl;
                exit(1);
            }


            node = std::shared_ptr<node_t>(new node_t(mcolor_t(left_child->get_data(), right_child->get_data(),
                                                       mcolor_t::Union), left_child, right_child));
            node->get_left_child()->set_parent(node.get());
            node->get_right_child()->set_parent(node.get());
        } else if (new_tree[0] == '{' && new_tree[new_tree.size() - 1] == '}') {
            node = build_subtree_leaf(new_tree, genome_number);
        } else {
            node = build_leaf(new_tree, genome_number);
        }

        if (ind != str.size() && !node) {
            mcolor_to_name.insert(
                    std::make_pair(node->get_data(), str.substr((size_t) (ind + 1), str.length() - ind - 1)));
        }

        return node;
    }

    /**
     *
     * @return
     */
    static node_ptr build_subtree_leaf(std::string const &str, genome_number_t  const &genome_number) {
        size_t start = 1;
        mcolor_t color;

        for (size_t j = 1; j < str.size(); ++j) {
            if (str[j] == ',' || str[j] == '}') {
                std::string const &name = str.substr(start, j - start);
                auto it = genome_number.find(name);
                if (it == genome_number.end()) {
                    std::cerr << "Unknown genome in (sub)tree: " << name << std::endl;
                    exit(1);
                }
                color.insert(it->second);
                start = j + 1;
            }
        }

        node_ptr node(new node_t(color));
        return node;
    }

    /**
     *
     * @return
     */
    static node_ptr build_leaf(std::string const &str, genome_number_t  const &genome_number) {
        auto it = genome_number.find(str);
        if (it == genome_number.end()) {
            std::cerr << "Unknown genome in (sub)tree: " << str << std::endl;
            exit(1);
        }
        node_ptr node(new node_t(mcolor_t(it->second)));
        return node;
    }

    /**
     *
     * @return
     */
    static size_t truncate_name(std::string const &str) {
        int i = 0;

        for (i = (int) (str.length() - 1); str[i] != ':' && str[i] != ')' && str[i] != '}' && i >= 0; --i);

        if (i < 0) {
            return str.length();
        }

        if (str[i] == ':') {
            return (size_t) i;
        } else {
            return (size_t) (i + 1);
        }
    }

    /**
     * Merges a given branch into tree
     * @return true if merged correctly, false if found an intersection
     */
    static bool merge_branch_into_node(node_ptr const &node, branch_t const &branch) {
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
                current_node = current_node->get_left_child();
                continue;
            }

            if (branch.node_right_includes(current_node)) {
                current_node = current_node->get_right_child();
                continue;
            }

            // Neither includes left or right, found intersection
            return false;
        }
    }
};

}

}

#endif //MGRA_TREE_BUILDER_HPP
