#ifndef JSON_GRAPH_HPP
#define JSON_GRAPH_HPP

namespace writer {

template<class graph_pack_t>
struct GraphJson {
  using mcolor_t = typename graph_pack_t::mcolor_t;
  using mularcs_t = typename graph_pack_t::mularcs_t;

  GraphJson()
  : infinity(0)
  {
  }

  // Save multiply breakpoint graph in json file
  void save_json_graph(std::string const & path_to_save, graph_pack_t const & graph, std::unordered_map<vertex_t, int> const & vertex_to_id);

private:
  size_t infinity;

  // Save multiedge in json file
  void describe_multiedge(std::ofstream & json, std::unordered_map<vertex_t, int> const & vertex_to_id, vertex_t const & x, mcolor_t const & color, vertex_t const & y);

  //Save vertex in json file
  void describe_vertex(std::ofstream & json, std::unordered_map<vertex_t, int> const & vertex_to_id, vertex_t const & v);

};

template<class graph_pack_t>
void GraphJson<graph_pack_t>::save_json_graph(std::string const & path_to_save, graph_pack_t const & graph_pack, std::unordered_map<vertex_t, int> const & vertex_to_id) {
  std::ofstream json(path::append_path(path_to_save, "graph.json"));

  //describe edges
  json << "{" << std::endl;
  json << "\t\"breakpoint_graph\": {" << std::endl;
  json << "\t\t\"edges\": [" << std::endl;

  std::unordered_set<vertex_t> processed;
  infinity = graph_pack.graph.size() + 1;
  bool need_comma = false;
  for (vertex_t const & x : graph_pack.graph) {
    mularcs_t const & mularcs = graph_pack.get_all_adjacent_multiedges(x);

    for (auto const & edge : mularcs) {
      if (processed.find(edge.first) != processed.end()) { //this edge is already processed
        continue;
      }
      if (need_comma) {
          json << "," << std::endl;
      }
      describe_multiedge(json, vertex_to_id, x, edge.second, edge.first);
      need_comma = true;
    }
    processed.insert(x);
  }

  json << std::endl;
  json << "\t\t]," << std::endl;
  json << "\t\t\"vertices\": [" << std::endl;

  //desribe vertices
  need_comma = false;
  for (vertex_t const & x : graph_pack.graph) {
      if (need_comma) {
          json << "," << std::endl;
      }
      describe_vertex(json, vertex_to_id, x);
      need_comma = true;
  }

  //describe infinity vertex?
  json << "," << std::endl;
  json << "\t\t\t{" << std::endl;
  json << "\t\t\t\t\"name\": oo," << std::endl;
  json << "\t\t\t\t\"v_id\": " << infinity << std::endl;
  json << "\t\t\t}" << std::endl;

  json << "\t\t]" << std::endl;
  json << "\t}" << std::endl;
  json << "}" << std::endl;
  json.close();
}

template<class graph_pack_t>
void GraphJson<graph_pack_t>::describe_multiedge(std::ofstream & json, std::unordered_map<vertex_t, int> const & vertex_to_id, vertex_t const & x, mcolor_t const & color, vertex_t const & y) {
  json << "\t\t\t{" << std::endl;

  json << "\t\t\t\t\"vertex1_id\": " << vertex_to_id.at(x) << "," << std::endl;
  if (y == Infty) {
    json << "\t\t\t\t\"vertex2_id\": " << infinity << "," << std::endl;
  } else {
    json << "\t\t\t\t\"vertex2_id\": " << vertex_to_id.at(y) << "," << std::endl;
  }

  json << "\t\t\t\t\"multicolor\": [" << std::endl;
  bool need_comma = false;
  for (auto const & match : color) {
    for (size_t i = 0; i < match.second; ++i) {
      if (need_comma) {
          json << "," << std::endl;
      }
      json << "\t\t\t\t\t" << match.first;
      need_comma = true;
    }
  }
  json << std::endl;
  json << "\t\t\t\t]" << std::endl;
  json << "\t\t\t}";
}

template<class graph_pack_t>
void GraphJson<graph_pack_t>::describe_vertex(std::ofstream & json, std::unordered_map<vertex_t, int> const & vertex_to_id, vertex_t const & v) {
    json << "\t\t\t{" << std::endl;
    json << "\t\t\t\t\"name\": \"" << v << "\"," << std::endl;
    json << "\t\t\t\t\"v_id\": " << vertex_to_id.at(v) << std::endl;
    json << "\t\t\t}";
}

}

#endif // JSON_GRAPH_HPP
