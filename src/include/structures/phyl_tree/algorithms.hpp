//
// Created by pavel on 12/16/15.
//

#ifndef MGRA_TREE_ALGORITHMS_HPP
#define MGRA_TREE_ALGORITHMS_HPP

namespace structure {

namespace phyl_tree {

template<class tree_t, class map_wgd_t>
inline void refine_tree_by_wgds(tree_t & tree, map_wgd_t const & wgds) {
    using node_t = typename tree_t::node_t;
    using node_ptr = typename tree_t::node_ptr;

    std::function<void(node_ptr)> add_wgd_nodes = [&](node_ptr node) -> void {
        if (node->has_childrens()) {
            for (auto children : *node) {
                add_wgd_nodes(children);
            }

            size_t ind = 0;
            for (auto children = node->begin(); children != node->end(); ++children, ++ind) {
                auto wgd = wgds.find(std::make_pair(node->get_data(), (*children)->get_data()));

                if (wgd != wgds.end()) {
                    auto child = *children;
                    node_ptr wgd_node(nullptr);

                    for (auto name = wgd->second.crbegin(); name != wgd->second.crend(); ++name) {
                        wgd_node = node_ptr(new node_t(*name, child));
                        child = wgd_node;
                    }

                    assert(static_cast<bool>(wgd_node));

                    node->replace_children(wgd_node, ind);
                }
            }
        }
    };

    add_wgd_nodes(tree.get_root());
}

template<class tree_t, class genome_number_t, class mcolor_t>
inline void get_T_multicolors(tree_t const &tree, genome_number_t const &genome_number,
                              std::map<mcolor_t, std::string> &color_to_name, mcolor_t &complete_color) {
    using node_t = typename tree_t::node_t;
    using node_ptr = typename tree_t::node_ptr;

    std::function<mcolor_t(node_ptr)> get_multicolors = [&](node_ptr node) -> mcolor_t {
        if (node->is_leaf()) {
            auto it = genome_number.find(node->get_data());

            if (it == genome_number.cend()) {
                ERROR("ERROR: undefined name in phylogenetic tree")
            }

            color_to_name.insert(std::make_pair(mcolor_t(it->second), node->get_data()));
            return mcolor_t(it->second);
        } else if (node->get_type() == node_t::whole_duplication) {
            assert(node->get_parent() != nullptr);
            assert(node->get_most_left_child() == node->get_most_right_child());
            mcolor_t result = get_multicolors(node->get_most_left_child());
            result.doubling();
            color_to_name.insert(std::make_pair(result, node->get_data()));
            return result;
        } else {
            assert(node->has_childrens());

            mcolor_t result;
            for (auto children : *node) {
                mcolor_t color = get_multicolors(children);
                result = mcolor_t(result, color, mcolor_t::Union);
            }

            color_to_name.insert(std::make_pair(result, node->get_data()));
            return result;
        }
    };

    complete_color = get_multicolors(tree.get_root());
}

template<class tree_t, class mcolor_t, class median_mcolor_t>
inline void get_median_colors(tree_t const &tree, mcolor_t const &complete_color,
                              std::unordered_map<std::string, mcolor_t> const & name_mcolor,
                              median_mcolor_t &median_colors) {
    using node_t = typename tree_t::node_t;
    using node_ptr = typename tree_t::node_ptr;

    std::function<void(node_ptr)> get_medians = [&](node_ptr node) -> void {
        for (auto children : *node) {
            get_medians(children);
        }

        if ((node->get_type() == node_t::classic || node->get_type() == node_t::whole_duplication) && !node->is_leaf()) {
            assert(node->has_childrens());
            assert(name_mcolor.count(node->get_most_left_child()->get_data()) != 0);
            assert(name_mcolor.count(node->get_most_right_child()->get_data()) != 0);

            mcolor_t left = name_mcolor.find(node->get_most_left_child()->get_data())->second;
            mcolor_t right = name_mcolor.find(node->get_most_right_child()->get_data())->second;
            mcolor_t current(left, right, mcolor_t::Union);
            mcolor_t result(complete_color, current, mcolor_t::Difference);

            if (node->get_parent() != nullptr) {
                median_colors.push_back(std::make_tuple(left, right, result));
            }
        }
    };

    get_medians(tree.get_root());
}
    /*
    template<class mcolor_t>
inline std::vector<Branch<mcolor_t>> break_tree_into_branches(mcolor_t const & complete_color, BinaryTree<mcolor_t> const & tree) {
    using node_ptr = typename BinaryTree<mcolor_t>::node_ptr;

    std::vector<Branch<mcolor_t>> result;
    std::queue<node_ptr> nodes_to_process;
    nodes_to_process.push(tree.get_root());

    while (!nodes_to_process.empty()) {
        auto node = nodes_to_process.front();
        auto node_color = node->get_data();
        nodes_to_process.pop();

        if (node_color != complete_color) {
            // Otherwise packed compliment is empty
            Branch<mcolor_t> branch;
            branch.init_by_complement(node_color, complete_color);
            result.push_back(branch);
        }

        if (node->has_left_child()) {
            nodes_to_process.push(node->get_left_child());
        }

        if (node->has_right_child()) {
            nodes_to_process.push(node->get_right_child());
        }
    }

    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());

    return result;
}*/

}

}

#endif //MGRA_TREE_ALGORITHMS_HPP
