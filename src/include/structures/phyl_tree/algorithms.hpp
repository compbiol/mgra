//
// Created by pavel on 12/16/15.
//

#ifndef MGRA_TREE_ALGORITHMS_HPP
#define MGRA_TREE_ALGORITHMS_HPP

#include "event/WGD.hpp"

namespace structure {

namespace phyl_tree {

template<class mcolor_t>
std::vector<Branch<mcolor_t>> break_tree_into_branches(mcolor_t const & complete_color, BinaryTree<mcolor_t> const & tree) {
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
}

template<class mcolor_t>
void refine_tree_by_wgds(BinaryTree<mcolor_t> & tree, std::map<Branch<mcolor_t>, event::WGD<mcolor_t>> const & wgds) {
    using node_t = typename BinaryTree<mcolor_t>::colored_node_t;
    using node_ptr = typename BinaryTree<mcolor_t>::node_ptr;

    std::queue<node_ptr> nodes_to_process;

    nodes_to_process.push(tree.get_root());
    while (!nodes_to_process.empty()) {
        auto node = nodes_to_process.front();
        nodes_to_process.pop();

        if (node->has_left_child()) {
            auto wgd = wgds.find(std::make_pair(node->get_data(), node->get_left_child().get_data()));
            if (wgd != wgds.end()) {
                node_ptr wgd_node(new node_t(wgd->get_corresponding_color(),
                                             node->get_left_child()));
                node->set_left_child(wgd_node);
                wgd_node->set_parent(node.get());
                node->get_left_child()->set_parent(wgd_node.get());
            }
            nodes_to_process.push(node->get_left_child());
        }

        if (node->has_right_child()) {
            auto wgd = wgds.find(std::make_pair(node->get_data(), node->get_right_child().get_data()));
            if (wgd != wgds.end()) {
                node_ptr wgd_node(new node_t(wgd->get_corresponding_color(),
                             node->get_right_child()));
                node->set_right_child(wgd_node);
                wgd_node->set_parent(node.get());
                node->get_right_child()->set_parent(wgd_node.get());
            }
            nodes_to_process.push(node->get_right_child());
        }
    }
}

}

}

#endif //MGRA_TREE_ALGORITHMS_HPP
