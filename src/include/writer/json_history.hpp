#ifndef JSON_HISTORY_HPP
#define JSON_HISTORY_HPP

namespace writer {

template<class graph_pack_t>
struct HistoryJson {
  using mcolor_t = typename graph_pack_t::mcolor_type;
  using twobreak_t = typename graph_pack_t::twobreak_t;
  using clone_t =  typename graph_pack_t::clone_t;
  using type_action = typename History<mcolor_t>::type_action;

  // Save history in json file
  void save_json_history(std::string const & path_to_save, graph_pack_t const & graph, std::unordered_map<vertex_t, int> const & vertex_to_id);

private:
  // Save twobreak in json file
  void describe_two_break(std::ofstream & json, std::unordered_map<vertex_t, int> const & vertex_to_id, twobreak_t const & action, size_t action_id);

  //Save clone in json file
  void describe_clone(std::ofstream & json, std::unordered_map<vertex_t, int> const & vertex_to_id, clone_t const & action, size_t action_id);

};

template<class graph_pack_t>
void HistoryJson<graph_pack_t>::save_json_history(std::string const & path_to_save, graph_pack_t const & graph_pack, std::unordered_map<vertex_t, int> const & vertex_to_id) {
  std::ofstream json(path::append_path(path_to_save, "history.json"));

  std::list<std::pair<type_action, size_t> > complete_history = graph_pack.history.get_complete_history();
  std::vector<twobreak_t> twobreak_history = graph_pack.history.get_twobreak_history();
  std::vector<clone_t> clone_history = graph_pack.history.get_clone_history();

  //fill operation ids
  std::vector<size_t> twobreak_ids;
  std::vector<size_t> clone_ids;
  size_t counter = 1;

  for (auto action = complete_history.rbegin(); action != complete_history.rend(); ++action) {
      if (action->first == type_action::twobreak_action) {
          twobreak_ids.push_back(counter);
      } else if (action->first == type_action::clone_action) {
          clone_ids.push_back(counter);
      }
      ++counter;
  }

  //describe two breaks
  json << "\"history\": {" << std::endl;
  json << "\t\"two_breaks\": [" << std::endl;
  bool need_comma = false;

  for (size_t i = 0; i < twobreak_history.size(); ++i){
      if (need_comma) {
          json << "," << std::endl;
      }
      describe_two_break(json, vertex_to_id, twobreak_history[i], twobreak_ids[i]);
      need_comma = true;
  }

  json << std::endl;
  json << "\t]," << std::endl;

  json << "\t\"clone\": [" << std::endl;

  //desribe clone
  need_comma = false;
  for (size_t i = 0; i < clone_history.size(); ++i) {
      if (need_comma) {
          json << "," << std::endl;
      }
      describe_clone(json, vertex_to_id, clone_history[i], clone_ids[i]);
      need_comma = true;
  }

  json << std::endl;
  json << "\t]" << std::endl;
  json << "}" << std::endl;
  json.close();
}

template<class graph_pack_t>
void HistoryJson<graph_pack_t>::describe_two_break(std::ofstream & json, std::unordered_map<vertex_t, int> const & vertex_to_id, twobreak_t const & action, size_t action_id) {
  json << "\t\t{" << std::endl;

  json << "\t\t\t\"operation_id\": " << action_id << "," << std::endl;

  json << "\t\t\t\"vertex1_id\": " << vertex_to_id.at(action.get_vertex(0)) << "," << std::endl;
  json << "\t\t\t\"vertex2_id\": " << vertex_to_id.at(action.get_vertex(1)) << "," << std::endl;
  json << "\t\t\t\"vertex3_id\": " << vertex_to_id.at(action.get_vertex(2)) << "," << std::endl;
  json << "\t\t\t\"vertex4_id\": " << vertex_to_id.at(action.get_vertex(3)) << "," << std::endl;

  json << "\t\t\t\"multicolor\": [" << std::endl;
  bool need_comma = false;
  for (auto const & match : action.get_multicolor()) {
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
void HistoryJson<graph_pack_t>::describe_clone(std::ofstream & json, std::unordered_map<vertex_t, int> const & vertex_to_id, clone_t const & action, size_t action_id) {
    json << "\t\t{" << std::endl;
    json << "\t\t\t\"operation_id\": " << action_id << "," << std::endl;
    json << "\t\t\t\"moter_vertex_id\": " << vertex_to_id.at(action.get_mother_arc().first) << "," << std::endl;
    json << "\t\t\t\"is_pseudo_mother_vertex\": " << action.is_have_pseudo_vertex() << "," << std::endl;
    json << "\t\t\t\"central_vertex1_id\": " << vertex_to_id.at(action.get_central_edge().first) << "," << std::endl;
    json << "\t\t\t\"central_vertex2_id\": " << vertex_to_id.at(action.get_central_edge().second) << "," << std::endl;
    json << "\t\t\t\"father_edges\": [" << std::endl;
    bool need_comma1 = false;
    for (auto const & edge : action.get_fathers()) {
        if (need_comma1) {
            json << "," << std::endl;
        }
        json << "\t\t\t\t{" << std::endl;
        json << "\t\t\t\t\t\"vertex_id\": " << vertex_to_id.at(edge.first) << "," << std::endl;
        json << "\t\t\t\t\t\"multicolor\": [" << std::endl;
        need_comma1 = true;
        bool need_comma2 = false;
        for (auto const & match : edge.second) {
            for (size_t i = 0; i < match.second; ++i) {
                if (need_comma2) {
                    json << "," << std::endl;
                }
                json << "\t\t\t\t\t\t" << match.first;
                need_comma2 = true;
            }
        }
        json << std::endl;
        json << "\t\t\t\t\t]" << std::endl;
        json << "\t\t\t\t}";
    }

    json << std::endl;
    json << "\t\t\t]," << std::endl;
    json << "\t\t\t\"multicolor\": [" << std::endl;
    need_comma1 = false;
    for (auto const & match : action.get_mother_arc().second) {
        for (size_t i = 0; i < match.second; ++i) {
            if (need_comma1) {
                json << "," << std::endl;
            }
            json << "\t\t\t\t" << match.first;
            need_comma1 = true;
        }
    }
    json << std::endl;
    json << "\t\t\t]" << std::endl;
    json << "\t\t}";
}

}

#endif // JSON_HISTORY_HPP
