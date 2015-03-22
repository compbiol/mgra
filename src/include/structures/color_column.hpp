//
// Created by Nikita Kartashov on 17/03/2015.
//

#ifndef MGRA_COLOR_PYRAMID_HPP_
#define MGRA_COLOR_PYRAMID_HPP_

#include <vector>
#include <unordered_map>

namespace structure {

  /**
  * Holds levels of color, with info specified for them
  * Level 0 should always be empty
  */
  template <class mcolor_t, class mcolor_info_t>
  struct ColorColumn {
    using scored_color_t = std::pair<mcolor_t, size_t>;
    using weighted_color_t = std::pair<mcolor_t, mcolor_info_t>;
    using level_t = std::vector<weighted_color_t>;
    using level_vector = std::vector<level_t>;
    using branch_t = std::pair<mcolor_t, mcolor_t>;
    using statistic_t = std::pair<branch_t, mcolor_t>;
    using position_t = typename mcolor_info_t::position_t;
    using position_vector = std::vector<position_t>;
    using color_info_t = mcolor_info_t;

    static const size_t NOT_FOUND = static_cast<size_t>(-1);

    ColorColumn() {
      m_levels.push_back({});
    }

    level_t& operator[](size_t index) {
      assert(index != 0);
      return m_levels[index];
    }

    weighted_color_t const& operator[](position_t const& position) const {
      assert(position.first != 0);
      return m_levels[position.first][position.second];
    }

    void insert(scored_color_t scored_color) {
      auto& color = scored_color.first;
      auto color_size = color.size();
      assert(color_size != 0);
      extend_to_size(color_size);
      m_levels[color_size].push_back(std::make_pair(color, mcolor_info_t(scored_color.second)));
    }

    size_t size() const {
      return m_levels.size();
    }

    bool exists(mcolor_t& color) const {
      return find(color) != NOT_FOUND;
    }

    size_t find(mcolor_t& color) const {
      auto color_size = color.size();
      if (m_levels.size() < (color_size + 1)) {
        return NOT_FOUND;
      }
      for (size_t i = 0; i != m_levels[color_size].size(); ++i) {
        if (m_levels[color_size][i].first == color) {
          return i;
        }
      }
      return NOT_FOUND;
    }

    mcolor_info_t& get_info(position_t const& position) {
      assert(position.first > 0);
      return m_levels[position.first][position.second].second;
    }

  private:
    void extend_to_size(size_t n) {
      while (m_levels.size() <= n) {
        m_levels.push_back({});
      }
    }

    level_vector m_levels;
  };
}

#endif //_MGRA_COLOR_PYRAMID_HPP_
