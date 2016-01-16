#ifndef BRUTE_FORCE_WITH_BLOSSOM_HPP
#define BRUTE_FORCE_WITH_BLOSSOM_HPP

#include "blossom5/PerfectMatching.h"

namespace algo {

template<class graph_pack_t>
struct ProcessWithBlossomV : public algo::AbsStage<graph_pack_t> {
  using int_edge_t = std::pair<int, int>;
  using weight_edge_t = std::tuple<int, int, int>;
  
  using mcolor_t = typename graph_pack_t::mcolor_t;
  
  using edge_t = typename graph_pack_t::edge_t;  
  using arc_t = typename graph_pack_t::arc_t; 
  using mularcs_t = typename graph_pack_t::mularcs_t; 
  using twobreak_t = typename graph_pack_t::twobreak_t;
  using transform_t = typename graph_pack_t::transform_t;

  ProcessWithBlossomV()
  : AbsStage<graph_pack_t>("Total resolve component with score function and blossomV", "wrap_blossom", post_stage_t, 3) 
  , max_number(0)
  , super_infinity("oooo")
  { 
  }

  bool run(graph_pack_t & graph_pack) override;

private:
  size_t process_component(graph_pack_t & graph_pack, std::set<vertex_t> const & vertex_set);
    
  void init_process_tc_edges(graph_pack_t const & graph_pack, std::set<vertex_t> const & vertex_set);
  std::vector<weight_edge_t> init(graph_pack_t const & graph_pack, std::set<vertex_t> const & vertex_set);

  int vertex_to_number(vertex_t const & a) {
    if (this->vertex2num.count(a) != 0) { 
      return this->vertex2num.find(a)->second;
    } else { 
      this->vertex2num.insert(std::make_pair(a, max_number));
      this->num2vertex.insert(std::make_pair(max_number, a));
      ++max_number;
      return this->vertex2num.find(a)->second;
    } 
  }

private:
  int max_number;
  vertex_t super_infinity;
  
  std::multimap<int, edge_t> queue_tc_edges; 
  std::map<edge_t, int> tc_edges;
  std::map<edge_t, int> score_edges;

  std::unordered_map<int, vertex_t> num2vertex;
  std::unordered_map<vertex_t, int> vertex2num;  

  std::unordered_set<int> pseudo_infinity_verteces;  
  std::unordered_set<int> double_pseudo_infinity_verteces;  

private:
  DECL_LOGGER("WrapperBlossomV");
};

template<class graph_pack_t>
bool ProcessWithBlossomV<graph_pack_t>::run(graph_pack_t & graph_pack) {
  size_t number_rear = 0; 

  utility::equivalence<vertex_t> components = graph_pack.split_on_components();
  std::map<vertex_t, std::set<vertex_t> > const & classes = components.get_eclasses<std::set<vertex_t> >();

  for (auto const & vertex_set : classes) { 
    std::set<vertex_t> component = vertex_set.second; 
    number_rear += process_component(graph_pack, component);
  }

  return (number_rear != 0);
} 

template<class graph_pack_t>
size_t ProcessWithBlossomV<graph_pack_t>::process_component(graph_pack_t & graph_pack, std::set<vertex_t> const & vertex_set) {
  size_t number_rear = 0;

  init_process_tc_edges(graph_pack, vertex_set);
  
  if (!tc_edges.empty()) {
    edge_t edge = queue_tc_edges.begin()->second;
    transform_t twobreaks = graph_pack.take_edge_on_color(edge.first, graph_pack.multicolors.get_complete_color(), edge.second);      
    for (twobreak_t const & twobreak : twobreaks) { 
      graph_pack.apply(twobreak);
      ++number_rear;
    }      
 } else {
    auto graph_edges = init(graph_pack, vertex_set);

    PerfectMatching pm(num2vertex.size(), graph_edges.size());
    
    for (auto const & edge : graph_edges) { 
      int x = 0; 
      int y = 0; 
      int score = 0;
      std::tie(x, y, score) = edge;
      pm.AddEdge(x, y, score);
    }

    INFO("Start Blossom V algorithm.")
    pm.Solve();
    INFO("Finish Blossom V algorithm.")

    auto get_min_edge_in_matching_lambda = [&]() -> edge_t { 
      edge_t min_edge;
      int score_min_edge = std::numeric_limits<int>::max();

      for (vertex_t const & v : vertex_set) {
        assert(vertex2num.find(v) != vertex2num.end());
        int num_x = vertex2num.find(v)->second;
        int num_y = pm.GetMatch(num_x);
        assert(num2vertex.find(num_y) != num2vertex.end());
        vertex_t u = num2vertex.find(num_y)->second;
        if (pseudo_infinity_verteces.count(num_y) != 0) {
          u = Infty;
        }

        auto iter = score_edges.find(edge_t(v, u)); 
        if (iter == score_edges.end()) {
          iter = score_edges.find(edge_t(u, v));
        } 

        int score = std::numeric_limits<int>::max();

        if (iter != score_edges.end()) {
          score = iter->second;
        }
      
        if (score < score_min_edge) { 
          score_min_edge = score; 
          min_edge = edge_t(v, u);
        }          
      } 

      return min_edge;
    };

    edge_t edge = get_min_edge_in_matching_lambda();
    transform_t twobreaks = graph_pack.take_edge_on_color(edge.first, graph_pack.multicolors.get_complete_color(), edge.second);      
    for (twobreak_t const & twobreak : twobreaks) { 
      graph_pack.apply(twobreak);
      ++number_rear;
    }      
  } 

  return number_rear;
}

template<class graph_pack_t>
void ProcessWithBlossomV<graph_pack_t>::init_process_tc_edges(graph_pack_t const & graph_pack, 
      std::set<vertex_t> const & vertex_set) {
  queue_tc_edges.clear();
  tc_edges.clear(); 
  std::unordered_set<vertex_t> processed;

  for (vertex_t const & v : vertex_set) {
    mularcs_t mularcs = graph_pack.get_all_adjacent_multiedges(v);

    for (auto const & arc : mularcs) { 
      if ((processed.count(arc.first) == 0) && (graph_pack.is_contain_T_consistent_color(v, arc.first))) {
        int score = graph_pack.calculate_cost(v, arc.first);
        queue_tc_edges.insert(std::make_pair(score, edge_t(v, arc.first)));
        tc_edges.insert(std::make_pair(edge_t(v, arc.first), score));
      } 
    }
    processed.insert(v);
  }

}

template<class graph_pack_t>
std::vector<std::tuple<int, int, int> > ProcessWithBlossomV<graph_pack_t>::init(graph_pack_t const & graph_pack,
  std::set<vertex_t> const & vertex_set) {
  max_number = 0; 

  score_edges.clear();
  tc_edges.clear();
  queue_tc_edges.clear();
  vertex2num.clear();
  num2vertex.clear(); 
  pseudo_infinity_verteces.clear();
  double_pseudo_infinity_verteces.clear();

  std::vector<weight_edge_t> graph_edges;
   //Create verteces and edges for graph
  std::unordered_set<vertex_t> processed;
  for (vertex_t const & v : vertex_set) {
    mularcs_t mularcs = graph_pack.get_all_adjacent_multiedges(v);
    for (auto const & arc : mularcs) { 
      if (processed.count(arc.first) == 0) { 
        int num_x = vertex_to_number(v);  
        int num_y = -1; 
        if (arc.first == Infty) { 
          vertex_t pseudo_vertex = "oo" + std::to_string(pseudo_infinity_verteces.size()) + "oo";
          num_y = vertex_to_number(pseudo_vertex);
          this->pseudo_infinity_verteces.insert(num_y);
          vertex_t double_pseudo_vertex = "ooo" + std::to_string(double_pseudo_infinity_verteces.size()) + "ooo";
          int num_z = vertex_to_number(double_pseudo_vertex);
          this->double_pseudo_infinity_verteces.insert(num_z);

          graph_edges.push_back(std::make_tuple(num_y, num_z, 200));   
        } else { 
          num_y = vertex_to_number(arc.first);
        }        
        assert(num_y != -1);

        int score = graph_pack.calculate_cost(v, arc.first);     
        score_edges.insert(std::make_pair(edge_t(v, arc.first), score));
        graph_edges.push_back(std::make_tuple(num_x, num_y, score)); 
      } 
    }
    processed.insert(v);
  }

  if (vertex2num.size() % 2 != 0 && !pseudo_infinity_verteces.empty()) {   
    int num_y = vertex_to_number(super_infinity);
    for (int num_x : double_pseudo_infinity_verteces) { 
      graph_edges.push_back(std::make_tuple(num_x, num_y, 200)); 
    } 
  }

  //std::cerr << vertex_set.size() << " " << pseudo_infinity_verteces.size() << " " << double_pseudo_infinity_verteces.size() << " " << vertex2num.size() << std::endl;

  assert((vertex2num.size() % 2 == 0) && (max_number % 2 == 0) && (num2vertex.size() % 2 == 0));

  return graph_edges;
}

} 
#endif