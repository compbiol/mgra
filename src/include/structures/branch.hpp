#ifndef BRANCH_HPP_
#define BRANCH_HPP_

#include <utility>
#include <iostream>

#include "tree.hpp"

namespace structure {

  template <class mcolor_t>
  struct Branch {
    using branch_t = std::pair<mcolor_t, mcolor_t>;
    using branch_vector = std::vector<branch_t>;
    using tree_t = BinaryTree<mcolor_t>;
    using tree_vector = std::vector<tree_t>;
    using node_t = typename tree_t::colored_node_t;
    using node_ptr = typename tree_t::node_ptr;
    using node_queue = std::queue<node_ptr>;

    /**
    * Merges a given branch into tree
    * @return true if merged correctly, false if found an intersection
    */
    static bool merge_branch_into_node(node_ptr const& node, branch_t const& branch) {
      node_ptr current_node = node;
      while (true) {
        if (current_node->is_leaf()) {
          if (current_node->get_data() == branch.first) {
            // Analogous to completeness check
            return false;
          }
          fill_node_from_branch(current_node, branch);
          return true;
        }

        if (node_left_includes(current_node, branch)) {
          current_node = current_node->get_left_child();
          continue;
        }

        if (node_right_includes(current_node, branch)) {
          current_node = current_node->get_right_child();
          continue;
        }

        // Neither includes left or right, found intersection
        return false;
      }
    }

    /**
    * Checks if two branch multicolors intersect
    */
    static bool do_intersect(branch_t const& left, branch_t const& right) {
      return !(left.first.includes(right.first) || left.second.includes(right.first) ||
          right.first.includes(left.first) || right.second.includes(left.first));
    }

    /**
    * Makes a new node with children colored as the given branch
    */
    static node_ptr node_from_branch(branch_t const& branch) {
      node_ptr new_node = std::make_shared<node_t>(mcolor_t(branch.first, branch.second, mcolor_t::Union));
      fill_node_from_branch(new_node, branch);
      return new_node;
    }

    static node_ptr fold_into_root_node(branch_vector& branches, bool balance_sort = false) {
      // No branches - no tree
      assert(!branches.empty());

      if (balance_sort) {
        // Doubtful fix, sort the branches, so that the root node is most balanced
        std::sort(std::begin(branches), std::end(branches), [](branch_t const& left, branch_t const& right) {
          auto left_diff = std::abs(static_cast<long>(left.first.size()) - static_cast<long>(left.second.size()));
          auto right_diff = std::abs(static_cast<long>(right.first.size()) - static_cast<long>(right.second.size()));
          return left_diff < right_diff;
        });
      }

      auto branch_iter = std::begin(branches);
      auto root_node = node_from_branch(*branch_iter);
      for (; branch_iter != std::end(branches); ++branch_iter) {
        merge_branch_into_node(root_node, *branch_iter);
      }
      root_node->set_name(
          cfg::get().mcolor_to_name(
              mcolor_t(root_node->get_left_child()->get_data(),
                  root_node->get_right_child()->get_data(),
                  mcolor_t::Union)));
      return root_node;
    }

    /**
    * Make node's children colored as the provided branch
    */
    static void fill_node_from_branch(node_ptr const& node, branch_t const& branch) {
      node->set_left_child(std::make_shared<node_t>(mcolor_t(node->get_data(), branch.first, mcolor_t::Intersection)));
      node->get_left_child()->set_parent(node.get());
      node->get_left_child()->set_name(cfg::get().mcolor_to_name(node->get_left_child()->get_data()));
      node->set_right_child(std::make_shared<node_t>(mcolor_t(node->get_data(), branch.second, mcolor_t::Intersection)));
      node->get_right_child()->set_parent(node.get());
      node->get_right_child()->set_name(cfg::get().mcolor_to_name(node->get_right_child()->get_data()));
    }

    /**
    * Suppose we have two branches b1 = Q1|Q2 & b2 = R1|R2
    * Let's say that "b1 left includes b2" if Q1 is a subset of R1 and R2 is a subset of Q2
    * Using the same logic, "b1 right includes b2" if Q2 is a subset of R2 and R1 is a subset of Q1
    * If we have a binary tree node, we can look at it as a branch and use the same logic
    */
    static bool node_left_includes(node_ptr const& node, branch_t const& branch) {
      return node->get_left_child()->get_data().includes(branch.first);
    }

    static bool node_right_includes(node_ptr const& node, branch_t const& branch) {
      return node->get_right_child()->get_data().includes(branch.first);
    }

    static branch_t flip(branch_t const& branch) {
      return branch_t(branch.second, branch.first);
    }

    static branch_vector break_trees_into_branches(
        tree_vector const& trees,
        mcolor_t complete_color = cfg::get().complete_color()) {
      branch_vector result;

      for (auto& tree: trees) {
        node_queue nodes_to_process;
        nodes_to_process.push(tree.get_root());
        while (!nodes_to_process.empty()) {
          auto node = nodes_to_process.front();
          auto node_color = node->get_data();
          nodes_to_process.pop();
          result.push_back(node_color.packed_compliment(complete_color));
          if (!node->is_leaf()) {
            nodes_to_process.push(node->get_left_child());
            nodes_to_process.push(node->get_right_child());
          }
        }
      }

      std::sort(result.begin(), result.end());
      result.erase(std::unique(result.begin(), result.end()), result.end());
      return result;
    }
  };
}

#endif