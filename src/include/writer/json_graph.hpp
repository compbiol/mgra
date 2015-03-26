#ifndef JSON_GRAPH_HPP
#define JSON_GRAPH_HPP

namespace writer {

template<class graph_pack_t>
struct GraphJson {
  using mcolor_t = typename graph_pack_t::mcolor_type;
  using edge_t = typename graph_pack_t::edge_t;
  using mularcs_t = typename graph_pack_t::mularcs_t;

  GraphJson()
  : infinity(0)
  {
  }

  void open(std::string const & path) {
    m_path = path;
  }

  // Save multiply breakpoint graoh in json file
  void save_json_graph(graph_pack_t const & graph, size_t stage);

private:
  std::string m_path;
  size_t infinity;
  std::unordered_map<vertex_t, int> vertex_to_id;

  // Save multiedge in json file
  void describe_multiedge(std::ofstream & json, vertex_t const & x, mcolor_t const & color, vertex_t const & y);

  //Save vertex in json file
  void describe_vertex(std::ofstream & json, vertex_t const & v);

};

template<class graph_pack_t>
void GraphJson<graph_pack_t>::save_json_graph(graph_pack_t const & graph_pack, size_t stage) {
  std::string jsonname = "stage" + std::to_string(stage) + ".json";
  std::ofstream json(path::append_path(m_path, jsonname));

  //fill map for vertices
  int id = 1;
  for (vertex_t const & x : graph_pack.graph) {
      vertex_to_id.insert(std::make_pair(x, id));
      ++id;
  }

  //describe edges
  json << "\"breakpoint_graph\": {" << std::endl;
  json << "\t\"edges\": [" << std::endl;

  std::unordered_set<vertex_t> processed;
  infinity = graph_pack.graph.size();
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
      describe_multiedge(json, x, edge.second, edge.first);
      need_comma = true;
    }
    processed.insert(x);
  }

  json << std::endl;
  json << "\t]," << std::endl;
  json << "\t\"vertices\": [" << std::endl;

  //desribe vertices
  need_comma = false;
  for (vertex_t const & x : graph_pack.graph) {
      if (need_comma) {
          json << "," << std::endl;
      }
      describe_vertex(json, x);
      need_comma = true;
  }

  //describe infinity vertex?
  json << "," << std::endl;
  json << "\t\t{" << std::endl;
  json << "\t\t\t\"name\": \"Infinity\"," << std::endl;
  json << "\t\t\t\"v_id\": " << infinity << std::endl;
  json << "\t\t}" << std::endl;

  json << "\t]" << std::endl;
  json << "}" << std::endl;
  json.close();
}

template<class graph_pack_t>
void GraphJson<graph_pack_t>::describe_multiedge(std::ofstream & json, vertex_t const & x, mcolor_t const & color, vertex_t const & y) {
  json << "\t\t{" << std::endl;

  json << "\t\t\t\"vertex1_id\": " << vertex_to_id[x] << "," << std::endl;
  if (y == Infty) {
    json << "\t\t\t\"vertex2_id\": " << infinity << "," << std::endl;
  } else {
    json << "\t\t\t\"vertex2_id\": " << vertex_to_id[y] << "," << std::endl;
  }

  json << "\t\t\t\"multicolor\": [" << std::endl;
  bool need_comma = false;
  for (auto const & match : color) {
    for (size_t i = 0; i < match.second; ++i) {
      if (need_comma) {
          json << "," << std::endl;
      }
      json << "\t\t\t\t" << match.first;
      need_comma = true;
    }
  }
  json << std::endl;
  json << "\t\t\t]" << std::endl;
  json << "\t\t}";
}

template<class graph_pack_t>
void GraphJson<graph_pack_t>::describe_vertex(std::ofstream & json, vertex_t const & v) {
    json << "\t\t{" << std::endl;
    json << "\t\t\t\"name\": \"" << v << "\"," << std::endl;
    json << "\t\t\t\"v_id\": " << vertex_to_id[v] << std::endl;
    json << "\t\t}";
}

}

#endif // JSON_GRAPH_HPP
