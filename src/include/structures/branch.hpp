#ifndef BRANCH_HPP__
#define BRANCH_HPP__

#include <utility>
#include <iostream>

#include "tree.hpp"

namespace structure {

  template <class mcolor_t>
  struct Branch {
    using branch_t = std::pair<mcolor_t, mcolor_t>;
    using tree_t = BinaryTree<mcolor_t>;
    using node_t = typename tree_t::colored_node_t;
    using node_ptr = typename tree_t::node_ptr;

    /**
    * Merges a given branch into tree
    * @return true if merged correctly, false if found an intersection
    */
    static bool merge_branch_into_node(node_ptr const& node, branch_t const& branch) {
      node_ptr current_node = node;
      while (!current_node->is_complete()) {
        if (current_node->is_leaf()) {
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
      // Found a complete node, no need to go deeper
      // have not merged the branch
      return false;
    }


    /**
    * Checks if two branch multicolors intersect
    */
    static bool do_intersect(branch_t const& left, branch_t const& right) {
      return !left.first.includes(right.first) && !left.second.includes(right.first);
    }

    /**
    * Makes a new node with children colored as the given branch
    */
    static node_ptr node_from_branch(branch_t const& branch) {
      node_ptr new_node = std::make_shared<node_t>(mcolor_t(branch.first, branch.second, mcolor_t::Union));
      fill_node_from_branch(new_node, branch);
      return new_node;
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


  private:

  };
}

#endif