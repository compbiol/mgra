#ifndef BRUTEFORCE_RECOVER_TREE_ALGORITHM_HPP__
#define BRUTEFORCE_RECOVER_TREE_ALGORITHM_HPP__

#include "recover_tree_algorithm.hpp"
#include "structures/branch.hpp"
#include "graph/graph_pack.hpp"

namespace algo {

  template <class graph_pack_t>
  struct BruteforceRecoverTreeAlgorithm : RecoverTreeAlgorithm<graph_pack_t> {
    using mcolor_t = typename RecoverTreeAlgorithm<graph_pack_t>::mcolor_t;
    using tree_t = typename RecoverTreeAlgorithm<graph_pack_t>::tree_t;
    using node_t = typename tree_t::colored_node_t;
    using node_ptr = std::shared_ptr<node_t>;
    using node_vector = std::vector<node_ptr>;
    using tree_ptr = typename RecoverTreeAlgorithm<graph_pack_t>::tree_ptr;
    using tree_vector = typename RecoverTreeAlgorithm<graph_pack_t>::tree_vector;
    using branch_t = std::pair<mcolor_t, mcolor_t>;
    using branch_vector = std::vector<branch_t>;
    using class_t = std::pair<branch_vector, size_t>;
    using class_vector = std::vector<class_t>;
    using statistic_t = std::pair<branch_t, size_t>;
    using statistic_vector = std::vector<statistic_t>;
    using BranchHelper = structure::Branch<mcolor_t>;

    BruteforceRecoverTreeAlgorithm(graph_pack_t& graph_pack) : m_graph_pack(graph_pack) {
    }

    tree_vector recover_trees() {
      std::map<branch_t, size_t> edges_statistics = m_graph_pack.stats.multiedges_count;
      // No statistics - no trees
      assert(!edges_statistics.empty());
      statistic_vector color_edges_pairs(std::begin(edges_statistics), std::end(edges_statistics));

      std::sort(std::begin(color_edges_pairs), std::end(color_edges_pairs),
          // Descending by number of edges sort
          [](statistic_t left, statistic_t right) {
            return left.second > right.second;
          });

      return build_trees(color_edges_pairs);
    }

    tree_vector build_trees(statistic_vector& branch_statistics) {
      class_vector tree_classes;

      for (auto& statistic: branch_statistics) {
        auto& new_branch = statistic.first;
        bool has_been_added = false;
        for (auto& cls: tree_classes) {
          bool new_branch_conflicts = false;
          for (auto& branch: cls.first) {
            if (BranchHelper::do_intersect(branch, new_branch)) {
              new_branch_conflicts = true;
              break;
            }
          }
          if (!new_branch_conflicts) {
            cls.first.push_back(new_branch);
            cls.second += statistic.second;
            has_been_added = true;
          }
        }
        if (!has_been_added) {
          tree_classes.push_back(class_t({statistic.first}, statistic.second));
        }
      }

      std::sort(std::begin(tree_classes), std::end(tree_classes),
          // Descending by class score
          [](class_t left, class_t right) {
            return left.second > right.second;
          });

      tree_vector results;
      results.reserve(tree_classes.size());

      std::transform(std::begin(tree_classes), std::end(tree_classes), std::back_inserter(results),
          [](class_t& cls_to_fold) {
            return fold_into_tree(cls_to_fold.first);
          });

      return results;
    }

  private:
    static tree_ptr fold_into_tree(branch_vector& branches) {
      // No branches - no tree
      assert(!branches.empty());
      auto branch_iter = std::begin(branches);
      auto root_node = BranchHelper::node_from_branch(*branch_iter);
      for (; branch_iter != std::end(branches); ++branch_iter) {
        BranchHelper::merge_branch_into_node(root_node, *branch_iter);
      }
      return std::make_shared<tree_t>(root_node);
    }

    graph_pack_t& m_graph_pack;
  };

}

#endif