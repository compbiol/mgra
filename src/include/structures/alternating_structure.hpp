//
// Created by Nikita Kartashov on 06/06/2015.
//

#ifndef MGRA_ALTERNATING_STRUCTURE_HPP
#define MGRA_ALTERNATING_STRUCTURE_HPP

#include <deque>

namespace structure {

  template <class vertex_t, class graph_pack_t>
  struct AlternatingStructure {
    using mcolor_t = typename graph_pack_t::mcolor_type;

    AlternatingStructure(vertex_t vertex) : m_is_cycle(false) {
      m_vertices.push_back(vertex);
    }

    void add_to_right(vertex_t vertex) {
      m_vertices.push_back(vertex);
    }

    void add_to_left(vertex_t vertex) {
      m_vertices.push_front(vertex);
    }

    void promote_to_cycle() {
      m_is_cycle = true;
    }

    bool is_cycle() const {
      return m_is_cycle;
    }

    size_t length() const {
      if (m_vertices.size() < 2) {
        return 0;
      }
      auto result = m_vertices.size();
      if (!m_is_cycle) {
        // Path has 1 edges less than vertices
        result -= 1;
      }
      return result;
    }

    bool is_color_alternating(graph_pack_t& graph_pack) {
      if (m_vertices.size() < 2) {
        return false;
      }
      m_colors.resize(2);
      m_colors[0] = graph_pack.get_all_multicolor_edge(m_vertices[0], m_vertices[1]);
      m_colors[1] = graph_pack.multicolors.get_complement_color(m_colors[0]);
      for (size_t i = 0; i != (m_vertices.size() - 1); ++i) {
        auto current_vertex = m_vertices[i];
        auto next_vertex = m_vertices[i + 1];
        auto current_color = graph_pack.get_all_multicolor_edge(current_vertex, next_vertex);
        if (current_color != m_colors[i % 2]) {
          return false;
        }
      }

      if (m_is_cycle) {
        return graph_pack.get_all_multicolor_edge(m_vertices[0], m_vertices[m_vertices.size() - 1]) ==
               m_colors[1];
      }

      return true;
    }

    std::pair<std::pair<mcolor_t, mcolor_t>, size_t> structure_digest() const {
      return std::make_pair(mcolor_t::pack(std::make_pair(m_colors[0], m_colors[1])), length());
    }

  private:
    std::deque<vertex_t> m_vertices;
    bool m_is_cycle;
    std::vector<mcolor_t> m_colors;
  };
}

#endif //MGRA_ALTERNATING_STRUCTURE_HPP
