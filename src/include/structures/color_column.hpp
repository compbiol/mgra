//
// Created by Nikita Kartashov on 17/03/2015.
//

#ifndef _MGRA_COLOR_PYRAMID_HPP_
#define _MGRA_COLOR_PYRAMID_HPP_

#include <vector>
#include <unordered_map>

namespace structure {

  template <class mcolor_t>
  struct ColorColumn {
    using level_t = std::unordered_map<mcolor_t, size_t>;
    using level_vector = std::vector<level_t>;
    using branch_t = std::pair<mcolor_t, mcolor_t>;
    using statistic_t = std::pair<branch_t, mcolor_t>;
    using weighted_color_t = std::pair<mcolor_t, size_t>;

    ColorColumn() {
      m_levels.push_back(level_t());
    }

    level_t& operator[](size_t index) {
      return m_levels[index];
    }

    void insert(weighted_color_t color) {
      auto color_size = color.first.size();
      extend_to_size(color_size);
      m_levels[color_size].insert(color);
    }

    void erase(mcolor_t const& color) {
      auto color_size = color.first.size();
      m_levels[color_size].erase(color);
    }

    bool exists(mcolor_t const& color) const {
      auto color_size = color.first.size();
      return m_levels[color_size].count(color);
    }

    size_t size() const {
      return m_levels.size();
    }

  private:
    void extend_to_size(size_t n) {
      while (m_levels.size() <= n) {
        m_levels.push_back(level_t());
      }
    }

    level_vector m_levels;
  };
}

#endif //_MGRA_COLOR_PYRAMID_HPP_
