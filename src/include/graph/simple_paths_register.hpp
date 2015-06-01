//
// Created by Nikita Kartashov on 02/06/2015.
//

#ifndef MGRA_SIMPLE_PATHS_REGISTER_HPP
#define MGRA_SIMPLE_PATHS_REGISTER_HPP


template <class mcolor_t, class vertex_t>
struct SimplePathsRegister {
  using color_pair_t = std::pair<mcolor_t, mcolor_t>;
  using scored_color_pairs_t = std::vector<std::pair<color_pair_t, size_t>>;

  struct SimplePath {
    struct ExtensibleVertex {
      ExtensibleVertex(vertex_t vertex, mcolor_t outcoming_edge_color) :
          m_vertex(vertex), m_outcoming_edge_color(outcoming_edge_color) { }

      void extend(vertex_t vertex, mcolor_t outcoming_edge_color) {
        m_vertex = vertex;
        m_outcoming_edge_color = outcoming_edge_color;
      }

      mcolor_t& outcoming_edge_color() {
        return m_outcoming_edge_color;
      }

      bool equals(vertex_t vertex) {
        return m_vertex == vertex;
      }

    private:
      vertex_t m_vertex;
      mcolor_t m_outcoming_edge_color;
    };

    SimplePath(vertex_t left, vertex_t right, color_pair_t& colors) :
        m_colors(colors),
        m_length(1),
        m_left(left, m_colors.second),
        m_right(right, m_colors.second) {
    }

    void extend(vertex_t old_vertex,
                vertex_t new_vertex,
                mcolor_t& color) {
      auto vertex_to_extend = m_left.equals(old_vertex) ? m_left : m_right;
      if (vertex_to_extend.outcoming_edge_color() != color) {
        return;
      }
      vertex_to_extend.extend(new_vertex, alternate(color));
      ++m_length;
    }

    color_pair_t colors() const {
      return m_colors;
    }

    size_t length() const {
      return m_length;
    }

  private:
    mcolor_t& alternate(mcolor_t color) {
      if (color == m_colors.first) {
        return m_colors.second;
      }
      return m_colors.first;
    }

    color_pair_t m_colors;
    size_t m_length;
    ExtensibleVertex m_left;
    ExtensibleVertex m_right;
  };

  void new_vertex_pair(vertex_t left, vertex_t right, color_pair_t& colors) {
    auto left_iter = m_ends.find(left);
    auto right_iter = m_ends.find(right);

    if (left_iter != m_ends.end() && right_iter != m_ends.end()) {
      return;
    }

    if (left_iter == m_ends.end() && right_iter == m_ends.end()) {
      m_ends[left] = m_simple_paths.size();
      m_ends[right] = m_simple_paths.size();
      m_simple_paths.push_back(SimplePath(left, right, colors));
    } else {
      auto found_end = left_iter == m_ends.end() ? right_iter : left_iter;
      auto new_end = found_end->first == left ? right : left;
      const auto path_index = found_end->second;
      m_simple_paths[path_index].extend(found_end->first, new_end, colors.first);
      m_ends.erase(found_end);
      m_ends[new_end] = path_index;
    }
  }

  template <class inserter_t>
  void get_scored_color_pairs(inserter_t& out) {
    std::transform(m_simple_paths.begin(),
    m_simple_paths.end(), out, [](SimplePath const& path) {
          return std::make_pair(path.colors(), path.length());
        });
    m_ends.clear();
    m_simple_paths.clear();
  }
private:
  std::unordered_map<vertex_t, size_t> m_ends;
  std::vector<SimplePath> m_simple_paths;
};

#endif //MGRA_SIMPLE_PATHS_REGISTER_HPP
