#ifndef JSON_READER_H
#define JSON_READER_H

#include "defined.h"
#include "graph/multigraph.hpp"
#include "graph/history.hpp"
#include "event/TwoBreak.hpp"
#include "event/Clone.hpp"
#include "structures/mcolor.hpp"
#include "structures/mularcs.hpp"
#include "io/path_helper.hpp"
#include "structures/config_struct.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>


namespace reader {

template<class mcolor_t>
struct ReadJson {
  using edge_t = std::pair<vertex_t, vertex_t>;
  using mularc_t = structure::Mularcs<mcolor_t>;

  void read_jsons(std::string const & path_to_graph, std::string const & path_to_history, MultiGraph & mulgraph, History<mcolor_t> & history);

private:
  std::unordered_map<int, vertex_t> id_to_vertex;

  void read_json_graph(std::string const & path_to_file, MultiGraph & mulgraph);
  void read_json_history(std::string const & path_to_file, History<mcolor_t> & history);
};

template<class mcolor_t>
void ReadJson<mcolor_t>::read_jsons(std::string const & path_to_graph, std::string const & path_to_history, MultiGraph & mulgraph, History<mcolor_t> & history) {
    read_json_graph(path_to_graph, mulgraph);
    read_json_history(path_to_history, history);
}

template<class mcolor_t>
void ReadJson<mcolor_t>::read_json_graph(std::string const & path_to_file, MultiGraph & mulgraph) {
      std::ifstream input(path_to_file);
      boost::property_tree::ptree pt;
      boost::property_tree::read_json(input, pt);
      for (auto vertex : pt.get_child("breakpoint_graph.vertices")) {
          vertex_t v = vertex.second.get<vertex_t>("name");
          int id = vertex.second.get<int>("v_id");
          id_to_vertex.insert(std::make_pair(id, v));
          mulgraph.add_vertex(v);
      }
      for (auto edge : pt.get_child("breakpoint_graph.edges")) {
          vertex_t u = id_to_vertex.at(edge.second.get<int>("vertex1_id"));
          vertex_t v = id_to_vertex.at(edge.second.get<int>("vertex2_id"));
          for (auto color : edge.second.get_child("multicolor")) {
              size_t c = color.second.get_value<size_t>();
              mulgraph.add_edge(c, u, v);
          }
      }
  }
template<class mcolor_t>
void ReadJson<mcolor_t>::read_json_history(std::string const & path_to_file, History<mcolor_t> & history) {
      std::ifstream input(path_to_file);
      boost::property_tree::ptree pt;
      boost::property_tree::read_json(input, pt);
      for (auto twobreak : pt.get_child("history.two_breaks")) {
          vertex_t u1 = id_to_vertex.at(twobreak.second.get<int>("vertex1_id"));
          vertex_t v1 = id_to_vertex.at(twobreak.second.get<int>("vertex2_id"));
          vertex_t u2 = id_to_vertex.at(twobreak.second.get<int>("vertex3_id"));
          vertex_t v2 = id_to_vertex.at(twobreak.second.get<int>("vertex4_id"));
          edge_t e1 = std::make_pair(u1, v1);
          edge_t e2 = std::make_pair(u2, v2);
          std::map<size_t, size_t> mcolor;
          for (auto color : twobreak.second.get_child("multicolor")) {
              size_t c = color.second.get_value<size_t>();
              mcolor.insert(std::make_pair(c, c));
          }
          history.save_twobreak(event::TwoBreak<mcolor_t>(e1, e2, structure::Mcolor(mcolor)));
      }
      for (auto clone : pt.get_child("history.clone")) {
          vertex_t mother_vertex = id_to_vertex.at(clone.second.get<int>("mother_vertex_id"));
          bool is_pseudo_mother_vertex = clone.second.get<bool>("is_pseudo_mother_vertex");
          vertex_t central_u = id_to_vertex.at(clone.second.get<int>("central_vertex1_id"));
          vertex_t central_v = id_to_vertex.at(clone.second.get<int>("central_vertex2_id"));
          edge_t central_e = std::make_pair(central_u, central_v);
          mularc_t fathers;
          for (auto father_edge : clone.second.get_child("father_edges")) {
              vertex_t v = id_to_vertex.at(father_edge.second.get<int>("vertex_id"));
              std::map<size_t, size_t> mcolor;
              for (auto color : father_edge.second.get_child("multicolor")) {
                  size_t c = color.second.get_value<size_t>();
                  mcolor.insert(std::make_pair(c, c));
              }
              fathers.insert(v, structure::Mcolor(mcolor));
          }
          std::map<size_t, size_t> mcolor;
          for (auto color : clone.second.get_child("multicolor")) {
              size_t c = color.second.get_value<size_t>();
              mcolor.insert(std::make_pair(c, c));
          }

          //I'm not sure that it's really color of mother arc
          auto mother_arc = std::make_pair(mother_vertex, structure::Mcolor(mcolor));
          history.save_clone(event::Clone<mcolor_t>(central_e, fathers, mother_arc, is_pseudo_mother_vertex));
      }

  }
}

#endif // JSON_READER_H
