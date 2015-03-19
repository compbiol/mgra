//
// Created by Nikita Kartashov on 17/03/2015.
//

#ifndef _MGRA_DYNAMIC_RECOVER_TREE_ALGORITHM_HPP_
#define _MGRA_DYNAMIC_RECOVER_TREE_ALGORITHM_HPP_

#include <queue>

#include "recover_tree_algorithm.hpp"
#include "structures/branch.hpp"
#include "graph/graph_pack.hpp"
#include "structures/mcolor_hash.hpp"
#include "structures/color_column.hpp"
#include "structures/mcolor_info.hpp"

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
    using pyramid_t = structure::ColorColumn<mcolor_t, structure::McolorInfo<mcolor_t>>;
    using position_t = typename pyramid_t::position_t;
    using position_vector = typename pyramid_t::position_vector;
    using mcolor_info_t = typename pyramid_t::color_info_t;

    DynamicRecoverTreeAlgorithm(graph_pack_t& graph_pack) : m_graph_pack(graph_pack) {
    }

    tree_vector recover_trees() {
      auto statistics = StatisticsProducer<graph_pack_t>::make_statistics(m_graph_pack);

      pyramid_t color_pyramid;
      for (statistic_t& statistic: statistics) {
        color_pyramid.insert(std::make_pair(statistic.first.first, statistic.second));
      }

      Scoreboard scoreboard(3);

      // Bypass the level 1, because it's made of singleton colors
      for (size_t current_level_index = 2; current_level_index < color_pyramid.size(); ++current_level_index) {
        auto& current_level = color_pyramid[current_level_index];

        // Iterate on all colors in the current level, building subtrees on them
        for (size_t current_color_index = 0; current_color_index != current_level.size(); ++current_color_index) {
          auto& current_weighted_color = current_level[current_color_index];
          auto& current_color = current_weighted_color.first;

          // Try to break off a color and check if the remainder color is present
          for (size_t subcolor_level_index = current_level_index - 1; subcolor_level_index != 0; --subcolor_level_index) {
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

      auto branches = recover_branches(scoreboard.top().second, color_pyramid);
      return {std::make_shared<tree_t>(BranchHelper::fold_into_root_node(branches))};
    }

    struct Scoreboard {
      using scored_position_t = std::pair<size_t, position_t>;

      Scoreboard(size_t max_scores) {
        assert(max_scores != 0);
        m_scores.resize(max_scores);
      }

      scored_position_t const& top() const {
        return m_scores[0];
      }

      scored_position_t const& operator[](size_t index) const {
        return m_scores[index];
      }

      size_t min_score() const {
        return m_scores.back().first;
      }

      size_t max_score() const {
        return m_scores.front().first;
      }

      bool update(scored_position_t scored_position) {
        if (scored_position.first <= min_score()) {
          return false;
        }
        if (scored_position.first > max_score()) {
          std::rotate(std::begin(m_scores), std::begin(m_scores) + 1, std::end(m_scores));
          m_scores[0] = scored_position;
          return true;
        }

        // Lower bound cannot point to the last element, it would have been rejected by 1st if clause
        auto lower_bound = std::lower_bound(m_scores.rbegin(), m_scores.rend(), scored_position).base();
        std::rotate(lower_bound, lower_bound + 1, std::end(m_scores));

        return true;
      }

    private:
      std::vector<scored_position_t> m_scores;
    };

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

    graph_pack_t& m_graph_pack;
  };
}

#endif //_MGRA_DYNAMIC_RECOVER_TREE_ALGORITHM_HPP_
