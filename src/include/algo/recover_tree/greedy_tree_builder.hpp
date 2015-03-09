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

    // Vector should sorted in the descending format
    GreedyTreeBuilder(statistic_vector statistics): m_branch_statistics(statistics) {
      validate_statistics();

      if (m_branch_statistics.empty()) {
        ERROR("No branch statistics supplied")
        exit(4);
      }

      // Take the best statistic as a root node
      auto iter = m_branch_statistics.begin();
      construct_root_node(*iter);
      iter++;

      node_ptr current_node = m_root_node;

      for (;iter != m_branch_statistics.end(); ++iter) {
        bool need_to_recurse_further = true;
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

    void construct_root_node(statistic_t const& best_statistic) {
      m_root_node = std::make_shared<node_t>(mcolor_t(best_statistic.first.first,
          best_statistic.first.second, mcolor_t::Union));
      fill_node_from_branch(m_root_node, best_statistic.first);
    }

    void fill_node_from_branch(node_ptr const& node, branch_t const& branch) {
      node->set_left_child(std::make_shared<node_t>(mcolor_t(node->get_data(), branch.first, mcolor_t::Intersection)));
      node->get_left_child()->set_parent(node.get());
      node->set_right_child(std::make_shared<node_t>(mcolor_t(node->get_data(), branch.second, mcolor_t::Intersection)));
      node->get_right_child()->set_parent(node.get());
    }

    /**
    * Evaluates the branch's left color's realtionship with the current node children's ones
    */
    ColorRelationship evaluateRelationship(node_ptr const& node, branch_t const& branch) {
      // Check only left child, tree is binary
      if (!node->has_left_child()) {
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