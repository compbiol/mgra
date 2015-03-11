#ifndef GREEDY_TREE_BUILDER_HPP__
#define GREEDY_TREE_BUILDER_HPP__

#include <memory>
#include <cassert>
#include "structures/tree.hpp"

namespace algo {

  template <class tree_t>
  struct GreedyTreeBuilder {
    using node_t = typename tree_t::colored_node_t;
    using mcolor_t = typename node_t::multicolor_t;
    using node_ptr = typename tree_t::node_ptr;
    using branch_t = std::pair<mcolor_t, mcolor_t>;
    using statistic_t = std::pair<branch_t, size_t>;
    using statistic_vector = std::vector<statistic_t>;

    // Vector should be sorted in descending format
    GreedyTreeBuilder(statistic_vector statistics): m_branch_statistics(statistics) {
      validate_statistics();

      if (m_branch_statistics.empty()) {
        ERROR("No branch statistics supplied")
        exit(4);
      }

      // Take the best statistic as a root node
      auto iter = std::begin(m_branch_statistics);
      construct_root_node(iter->first);
      iter++;

      for (; iter != std::end(m_branch_statistics); ++iter) {
        node_ptr current_node = m_root_node;
        bool need_to_recurse_further = !current_node->is_complete();
        ColorRelationship relationship = Leaf;
        while (need_to_recurse_further && relationship != Intersects) {
          relationship = evaluateRelationship(current_node, iter->first);
          switch (relationship) {
            case SubsetOfLeft:
              current_node = current_node->get_left_child();
              break;
            case SubsetOfRight:
              current_node = current_node->get_right_child();
              break;
            case Leaf:
              fill_node_from_branch(current_node, iter->first);
              need_to_recurse_further = false;
              break;
            default: break;
          }
          need_to_recurse_further = need_to_recurse_further && !current_node->is_complete();
        }
      }
    }

    node_ptr get_result() {
      return m_root_node;
    }

  private:
    // Describes how branch's colors relate to the current node's children
    enum ColorRelationship {
      Leaf,
      Intersects,
      SubsetOfLeft,
      SubsetOfRight
    };

    void validate_statistics() {
      // Maximum size_t value
      size_t previous_value = static_cast<size_t>(-1);
      for (auto& statistic: m_branch_statistics) {
        assert(previous_value >= statistic.second &&
            "Works only on sorted statistics");
        previous_value = statistic.second;
      }
    }

    void construct_root_node(branch_t const& best_branch) {
      m_root_node = std::make_shared<node_t>(mcolor_t(best_branch.first,
          best_branch.second, mcolor_t::Union));
      fill_node_from_branch(m_root_node, best_branch);
    }

    /**
    * Make node's children colored as the provided branch
    */
    void fill_node_from_branch(node_ptr const& node, branch_t const& branch) {
      node->set_left_child(std::make_shared<node_t>(mcolor_t(node->get_data(), branch.first, mcolor_t::Intersection)));
      node->get_left_child()->set_parent(node.get());
      node->set_right_child(std::make_shared<node_t>(mcolor_t(node->get_data(), branch.second, mcolor_t::Intersection)));
      node->get_right_child()->set_parent(node.get());
    }

    /**
    * Evaluates the branch's left color's relationship with the current node children's ones
    */
    ColorRelationship evaluateRelationship(node_ptr const& node, branch_t const& branch) {
      if (node->is_leaf()) {
        return Leaf;
      }
      if (node->get_left_child()->get_data().includes(branch.first)) {
        return SubsetOfLeft;
      }
      if (node->get_right_child()->get_data().includes(branch.first)) {
        return SubsetOfRight;
      }
      return Intersects;
    }

    statistic_vector& m_branch_statistics;
    node_ptr m_root_node;
  };
}

#endif