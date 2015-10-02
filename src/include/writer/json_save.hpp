#ifndef JSON_SAVE_HPP
#define JSON_SAVE_HPP

#include "writer/json_graph.hpp"
#include "writer/json_history.hpp"

namespace writer {

template<class graph_pack_t>
struct SaveJson {
  using mcolor_t = typename graph_pack_t::mcolor_type;
  using edge_t = typename graph_pack_t::edge_t;
  using mularcs_t = typename graph_pack_t::mularcs_t;

  void open(std::string const & path) {
    m_path = path;
  }

  void save_jsons(graph_pack_t const & graph, size_t stage, size_t round);

private:
  std::string m_path;
  writer::GraphJson<graph_pack_t> json_graph;
  writer::HistoryJson<graph_pack_t> json_history;
};

template<class graph_pack_t>
void SaveJson<graph_pack_t>::save_jsons(graph_pack_t const & graph, size_t stage, size_t round) {
  std::string jsons_dir = path::append_path(m_path, "stage" + std::to_string(stage) + "_round" + std::to_string(round));
  path::make_dir(jsons_dir);

  std::unordered_map<vertex_t, int> vertex_to_id;
  int id = 1;
  for (vertex_t const & x : graph.graph) {
      vertex_to_id.insert(std::make_pair(x, id));
      ++id;
  }
  vertex_to_id.insert(std::make_pair("oo", id));
  vertex_to_id.insert(std::make_pair("o0o", id));

  json_graph.save_json_graph(jsons_dir, graph, vertex_to_id);
  json_history.save_json_history(jsons_dir, graph, vertex_to_id);
}

}
#endif // JSON_SAVE_HPP
