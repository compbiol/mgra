//
// Created by Nikita Kartashov on 17/03/2015.
//

#ifndef MGRA_DYNAMIC_RECOVER_TREE_ALGORITHM_HPP_
#define MGRA_DYNAMIC_RECOVER_TREE_ALGORITHM_HPP_


#include "recover_tree_algorithm.hpp"
#include "structures/branch.hpp"
#include "graph/graph_pack.hpp"
#include "structures/mcolor_hash.hpp"
#include "structures/color_column.hpp"
#include "structures/mcolor_info.hpp"
#include "scoreboard.hpp"
#include "statistics_based_recover_tree_algorithm.hpp"

namespace algo {

  template <class graph_pack_t>
  struct DynamicRecoverTreeAlgorithm : StatisticsBasedRecoverTreeAlgorithm<graph_pack_t> {
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
    using pyramid_t = structure::ColorColumn<mcolor_t, structure::McolorInfo<mcolor_t>>;
    using position_t = typename pyramid_t::position_t;
    using position_vector = typename pyramid_t::position_vector;
    using mcolor_info_t = typename pyramid_t::color_info_t;
    using statistic_producer_ptr = std::shared_ptr<StatisticsProducer<graph_pack_t> > const&;

    DynamicRecoverTreeAlgorithm(graph_pack_t& graph_pack,
                                statistic_producer_ptr statistic_producer,
                                size_t returned_trees = 1) :
        StatisticsBasedRecoverTreeAlgorithm<graph_pack_t>(graph_pack, statistic_producer),
        m_returned_trees(returned_trees) {
    }

    tree_vector recover_trees() {
      auto statistics = this->m_statistic_producer->make_statistics();

      pyramid_t color_pyramid;
      for (statistic_t& statistic: statistics) {
        color_pyramid.insert(std::make_pair(statistic.first.first, statistic.second));
        color_pyramid.insert(std::make_pair(statistic.first.second, statistic.second));
      }
      color_pyramid.insert(std::make_pair(cfg::get().complete_color(), 0));

      Scoreboard<position_t> scoreboard(m_returned_trees);

      // Bypass the level 1, because it's made of singleton colors
      for (size_t current_level_index = 2; current_level_index < color_pyramid.size(); ++current_level_index) {
        auto& current_level = color_pyramid[current_level_index];

        // Iterate on all colors in the current level, building subtrees on them
        for (size_t current_color_index = 0; current_color_index != current_level.size(); ++current_color_index) {
          auto& current_weighted_color = current_level[current_color_index];
          auto& current_color = current_weighted_color.first;

          // Try to break off a color and check if the remainder color is present
          for (size_t subcolor_level_index = current_level_index - 1;
               subcolor_level_index != 0; --subcolor_level_index) {
            auto& subcolor_level = color_pyramid[subcolor_level_index];

            for (size_t subcolor_index = 0; subcolor_index != subcolor_level.size(); ++subcolor_index) {
              auto& weighted_subcolor = subcolor_level[subcolor_index];
              auto& subcolor = weighted_subcolor.first;

              // If subcolor is not in the current_color, try again
              if (!current_color.includes(subcolor)) {
                continue;
              }

              auto remainder_color = mcolor_t(current_color, subcolor, mcolor_t::Difference);
              auto remainder_color_index = color_pyramid.find(remainder_color);

              // If remainder is not present, try again
              if (remainder_color_index == color_pyramid.NOT_FOUND) {
                continue;
              }

              // Found subcolor + remainder_color pair, compute scores, compare to existing
              auto subcolor_position = std::make_pair(subcolor_level_index, subcolor_index);
              auto remainder_color_position = std::make_pair(remainder_color.size(), remainder_color_index);
              auto& current_info = current_weighted_color.second;
              auto& subcolor_info = weighted_subcolor.second;
              auto& remainder_info = color_pyramid.get_info(remainder_color_position);
              auto new_info = current_info.make_from_subtree_info(
                  subcolor_info,
                  remainder_info,
                  subcolor_position,
                  remainder_color_position);
              current_info.update_if_better(new_info);
              scoreboard.update(std::make_pair(new_info.get_whole_tree_score(),
                                               position_t(current_level_index, current_color_index)));
            }
          }
        }
      }

      std::vector<branch_vector> branches_to_fold;
      branches_to_fold.reserve(m_returned_trees);
      for (size_t i = 0; i < m_returned_trees; ++i) {
        branches_to_fold.push_back(recover_branches(scoreboard[i].second, color_pyramid));
      }

      tree_vector results;

      std::transform(std::begin(branches_to_fold), std::end(branches_to_fold), std::back_inserter(results),
                     [](branch_vector& branches) {
                       auto root_node = BranchHelper::fold_into_root_node(branches);
                       return std::make_shared<tree_t>(root_node);
                     });

      return results;
    }

  private:
    branch_vector recover_branches(position_t const& root_position, pyramid_t const& pyramid) {
      auto& root_node_info = pyramid[root_position].second;
      auto complete_color = cfg::get().complete_color();
      std::queue<position_t> position_queue;
      position_queue.push(root_node_info.get_subtree_positions().first);
      position_queue.push(root_node_info.get_subtree_positions().second);
      branch_vector result;
      while (!position_queue.empty()) {
        auto position = position_queue.front();
        position_queue.pop();
        auto& color = pyramid[position].first;
        auto& next_info = pyramid[position].second;
        if (next_info.has_subtrees()) {
          position_queue.push(next_info.get_subtree_positions().first);
          position_queue.push(next_info.get_subtree_positions().second);
        }
        result.push_back(branch_t(color, mcolor_t(complete_color, color, mcolor_t::Difference)));
      }

      return result;
    }

    const size_t m_returned_trees;
  };
}

#endif // MGRA_DYNAMIC_RECOVER_TREE_ALGORITHM_HPP_
