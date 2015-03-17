//
// Created by Nikita Kartashov on 17/03/2015.
//

#ifndef _MGRA_DYNAMIC_RECOVER_TREE_ALGORITHM_HPP_
#define _MGRA_DYNAMIC_RECOVER_TREE_ALGORITHM_HPP_

#include "recover_tree_algorithm.hpp"
#include "structures/branch.hpp"
#include "graph/graph_pack.hpp"
#include "structures/mcolor_hash.hpp"
#include "structures/color_column.hpp"

namespace algo {

  template <class graph_pack_t>
  struct DynamicRecoverTreeAlgorithm : RecoverTreeAlgorithm<graph_pack_t> {
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
    using color_set = std::unordered_set<mcolor_t>;
    using pyramid_t = structure::ColorColumn<mcolor_t>;

    DynamicRecoverTreeAlgorithm(graph_pack_t& graph_pack) : m_graph_pack(graph_pack) {
    }

    tree_vector recover_trees() {
      auto statistics = StatisticsProducer<graph_pack_t>::make_statistics(m_graph_pack);

      pyramid_t color_pyramid;
      for (statistic_t& statistic: statistics) {
        color_pyramid.insert(std::make_pair(statistic.first.first, statistic.second));
      }

      for (size_t level_index = 1; level_index < color_pyramid.size(); ++level_index) {
        
      }

      return {};
    }


  private:
    graph_pack_t& m_graph_pack;
  };
}

#endif //_MGRA_DYNAMIC_RECOVER_TREE_ALGORITHM_HPP_
