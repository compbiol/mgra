#ifndef GREEDY_RECOVER_TREE_ALGORITHM_HPP__
#define GREEDY_RECOVER_TREE_ALGORITHM_HPP__

#include "recover_tree_algorithm.hpp"
#include "graph/graph_pack.hpp"

namespace algo {

  template <class graph_pack_t>
  struct GreedyRecoverTreeAlgorithm : RecoverTreeAlgorithm<graph_pack_t> {
    using mcolor_t = typename RecoverTreeAlgorithm<graph_pack_t>::mcolor_t;
    using tree_t = typename RecoverTreeAlgorithm<graph_pack_t>::tree_t;
    using node_t = typename tree_t::colored_node_t;
    using node_ptr = std::shared_ptr<node_t>;
    using node_vector = std::vector<node_ptr>;
    using tree_ptr = typename RecoverTreeAlgorithm<graph_pack_t>::tree_ptr;
    using tree_vector = typename RecoverTreeAlgorithm<graph_pack_t>::tree_vector;
    using branch_t = std::pair<mcolor_t, mcolor_t>;
    using statistic_t = std::pair<branch_t, size_t>;
    using statistic_vector = std::vector<statistic_t>;

    GreedyRecoverTreeAlgorithm(graph_pack_t& graph_pack) : m_graph_pack(graph_pack) {
    }

    tree_vector recover_trees() {
      std::map<branch_t, size_t> edges_statistics = m_graph_pack.stats.multiedges_count;
      statistic_vector color_edges_pairs(std::begin(edges_statistics), std::end(edges_statistics));

      std::sort(std::begin(color_edges_pairs), std::end(color_edges_pairs),
          // Descending by number of edges sort
          [](statistic_t left, statistic_t right) {
            return left.second > right.second;
          });

      return build_trees(color_edges_pairs);
    }

    std::vector<tree_ptr> build_trees(statistic_vector& branch_statistics) {
      // Vector should be sorted in descending format
      validate_statistics(branch_statistics);

      if (branch_statistics.empty()) {
        ERROR("No branch statistics supplied")
        exit(4);
      }

      // Take the best statistic as a root node
      auto iter = std::begin(branch_statistics);
      auto root_node = construct_root_node(iter->first);
      iter++;

      for (; iter != std::end(branch_statistics); ++iter) {
        node_ptr current_node = root_node;
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
            default:
              break;
          }
          need_to_recurse_further = need_to_recurse_further && !current_node->is_complete();
        }
      }
      return {std::make_shared<tree_t>(root_node)};
    }

  private:
    /**
    * Describes how branch's colors relate to the current node's children
    */
    enum ColorRelationship {
      Leaf,
      Intersects,
      SubsetOfLeft,
      SubsetOfRight
    };

    /**
    * Checks if statistics are in descending order
    */
    void validate_statistics(statistic_vector const& branch_statistics) {
      // Maximum size_t value
      size_t previous_value = static_cast<size_t>(-1);
      for (auto const& statistic: branch_statistics) {
        assert(previous_value >= statistic.second &&
            "Works only on sorted statistics");
        previous_value = statistic.second;
      }
    }

    node_ptr construct_root_node(branch_t const& best_branch) {
      node_ptr root_node = std::make_shared<node_t>(mcolor_t(best_branch.first,
          best_branch.second, mcolor_t::Union));
      fill_node_from_branch(root_node, best_branch);
      return root_node;
    }

    /**
    * Make node's children colored as the provided branch
    */
    void fill_node_from_branch(node_ptr const& node, branch_t const& branch) {
      node->set_left_child(std::make_shared<node_t>(mcolor_t(node->get_data(), branch.first, mcolor_t::Intersection)));
      node->get_left_child()->set_parent(node.get());
      node->get_left_child()->set_name(cfg::get().mcolor_to_name(node->get_left_child()->get_data()));
      node->set_right_child(std::make_shared<node_t>(mcolor_t(node->get_data(), branch.second, mcolor_t::Intersection)));
      node->get_right_child()->set_parent(node.get());
      node->get_right_child()->set_name(cfg::get().mcolor_to_name(node->get_right_child()->get_data()));
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

    graph_pack_t& m_graph_pack;
  };
}

#endif