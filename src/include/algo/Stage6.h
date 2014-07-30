#ifndef STAGE6_H_ 
#define STAGE6_H_ 

#include "blossom5/PerfectMatching.h"

template<class graph_t>
struct Algorithm<graph_t>::BruteForce : public Algorithm<graph_t>::Stage {
  typedef Stage base;
  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  
  explicit BruteForce(std::shared_ptr<graph_t> const & graph)
  : Stage(graph) 
  {
  }
  
  BruteForce(std::shared_ptr<graph_t> const & graph, bool flag = true)
  : Stage(graph) 
  , use_library(flag)
  {
  }
  
  bool do_action() override;

  std::string get_name() override {
    return "Resolve component with score function."; 
  }

private:
  size_t take_edge_on_color(vertex_t const & x, mcolor_t const & color, vertex_t const & y);
  
  std::multimap<size_t, arc_t> create_minimal_matching(std::set<vertex_t> const & vertex_set); 
  std::multimap<size_t, arc_t> wrapper_create_minimal_matching(std::set<vertex_t> const & vertex_set);

private: 
  bool use_library;
};

template<class graph_t>
bool Algorithm<graph_t>::BruteForce::do_action() { 
    size_t number_break = 0;

  std::map<vertex_t, std::set<vertex_t> > const & classes = graph->split_on_components().get_eclasses();
  for (auto const & vertex_set : classes) { 
 
    //std::cerr << "Component have size " << vertex_set.second.size() << std::endl;

    if (vertex_set.second.size() <= max_size_component) { 

      /*std::cerr << "Start process component" << std::endl;
      for (auto const v : vertex_set.second) {
        std::cerr << v << " "; 
      } 
      std::cerr << std::endl;*/

      std::multimap<size_t, arc_t> matching = create_minimal_matching(vertex_set.second);
      /*for (auto const & edge : matching) { 
        std::cerr << "(" << edge.second.first << ", " << edge.second.second << "):" << edge.first << std::endl;
      }*/

      //std::cerr << "Start tested matching" << std::endl;
      //std::multimap<size_t, arc_t> matching = wrapper_create_minimal_matching(vertex_set.second);
      /*for (auto const & edge : matching) { 
        std::cerr << "(" << edge.second.first << ", " << edge.second.second << "):" << edge.first << std::endl;
      }*/
      
      //std::cerr << "Start process minimal matching" << std::endl;
      auto edge = matching.cbegin()->second;
      for (auto weight_edge: matching) { 
        if (weight_edge.second.first != Infty && weight_edge.second.second != Infty) { 
          edge = weight_edge.second;
          break;
        }
      }
      number_break += take_edge_on_color(edge.first, graph->get_complete_color(), edge.second);
    } 
  }  
  return (number_break != 0);
} 

template<class graph_t>
std::multimap<size_t, arc_t> Algorithm<graph_t>::wrapper_create_minimal_matching(std::set<vertex_t> const & vertex_set) {
  std::map<arc_t, std::pair<size_t, mcolor_t> > weight_edges; 
  
  std::unordered_map<std::string, int> vertex2num; 
  std::unordered_map<int, std::string> num2vertex;
  int max_number = 0; 
  
  for(vertex_t const & v : vertex_set) { 
    mularcs_t mularcs = graph->get_adjacent_multiedges(v);
    num2vertex.insert(std::make_pair(max_number, v)) ;
    vertex2num.insert(std::make_pair(v, max_number));
    ++max_number;

    for (auto const & arc: mularcs) {
      if (weight_edges.count(std::make_pair(v, arc.first)) == 0 && weight_edges.count(std::make_pair(arc.first, v)) == 0) { 
        size_t sz = graph->calculate_cost(v, arc.first);
        weight_edges.insert(std::make_pair(std::make_pair(v, arc.first), std::make_pair(sz, arc.second)));
      }
    }  
  }

  std::multimap<size_t, arc_t> minimal_matching; 
  std::unordered_set<std::string> processed;

  for (auto const &edge : weight_edges) { 
    if (graph->is_postponed_deletions(edge.first.first, edge.first.second) || (graph->is_T_consistent_color(edge.second.second) && !graph->is_vec_T_consistent_color(edge.second.second)) || edge.second.second.includes(graph->get_root_color())) { 
      processed.insert({edge.first.first, edge.first.second});
      processed.erase(Infty);
      minimal_matching.insert(std::make_pair(edge.second.first, edge.first));
    }
  }

  PerfectMatching pm(vertex_set.size(), weight_edges.size());
  for (auto const & edge : weight_edges) {
    if (processed.count(edge.first.first) == 0 && processed.count(edge.first.second) == 0) {
      if (edge.first.first == Infty) {
        num2vertex.insert(std::make_pair(max_number, edge.first.first));
        pm.AddEdge(max_number, vertex2num.find(edge.first.second)->second, edge.second.first);
        ++max_number;
      } else if (edge.first.second == Infty) { 
        num2vertex.insert(std::make_pair(max_number, edge.first.second));
        pm.AddEdge(vertex2num.find(edge.first.second)->second, max_number, edge.second.first);
        ++max_number;
      } else { 
        pm.AddEdge(vertex2num.find(edge.first.first)->second, vertex2num.find(edge.first.second)->second, edge.second.first);
      } 
    } 
  }

  pm.Solve();

  for (vertex_t const & v : vertex_set) { 
    if (processed.count(v) == 0) {
      int result = pm.GetMatch(vertex2num.find(v)->second);
      std::string u = num2vertex.find(result)->second;
      processed.insert({u, v});

      auto iter = weight_edges.find(std::make_pair(v, u));
      if (iter != weight_edges.end()) { 
        minimal_matching.insert(std::make_pair(iter->second.first, std::make_pair(v, u)));
      } else {
        iter = weight_edges.find(std::make_pair(u, v));
        if (iter != weight_edges.end()) { 
          minimal_matching.insert(std::make_pair(iter->second.first, std::make_pair(u, v)));
        } else {
          std::cerr << "ERROR: problem with minimal matching" << std::endl; 
          assert(false);
        }
      }
    }
  }
  
  return minimal_matching;  
} 

template<class graph_t>
std::multimap<size_t, arc_t> Algorithm<graph_t>::create_minimal_matching(std::set<vertex_t> const & vertex_set) {
  std::map<arc_t, std::pair<size_t, mcolor_t> > weight_edges; 

  for(vertex_t const & v : vertex_set) { 
    mularcs_t mularcs = graph->get_adjacent_multiedges(v);
    for (auto const & arc: mularcs) {
      if (weight_edges.count(std::make_pair(v, arc.first)) == 0 && weight_edges.count(std::make_pair(arc.first, v)) == 0) { 
        size_t sz = graph->calculate_cost(v, arc.first);
      	//std::cerr << "Calculate cost " << v << " " << arc.first << " have " << sz << std::endl;
      	weight_edges.insert(std::make_pair(std::make_pair(v, arc.first), std::make_pair(sz, arc.second)));
      }
    } 
  } 

  std::unordered_set<vertex_t> processed;
  std::set<arc_t> current; 
  //std::cerr << "Size " << vertex_set.size() << std::endl;
  for (auto const &edge : weight_edges) { 
    bool is_have_tc = false; 
    auto edge_colors = graph->split_color(edge.second.second);
    
    for (auto const & color : edge_colors) { 
      if (graph->is_T_consistent_color(color) && !graph->is_vec_T_consistent_color(color)) { 
        is_have_tc = true; 
        break;
      }
    }

    if (graph->is_postponed_deletion(edge.first.first, edge.first.second) || is_have_tc) { //(graph->is_T_consistent_color(edge.second.second) && !graph->is_vec_T_consistent_color(edge.second.second)) || edge.second.second.includes(graph->get_root_color())) { 
      current.insert(edge.first); 
      processed.insert({edge.first.first, edge.first.second});
      processed.erase(Infty);
    } 
  } 

  std::set<std::set<arc_t> > answer;
  if (processed.size() == vertex_set.size()) { 
    answer.insert(current);
  } else { 
    std::function<void()> findLambda = [&] () -> void { 
      for (auto const &e : weight_edges) { 
        if (processed.count(e.first.first) == 0 && processed.count(e.first.second) == 0) { 
          processed.insert({e.first.first, e.first.second});
          processed.erase(Infty);
      	  current.insert(e.first);
          findLambda();
      	  current.erase(e.first);
      	  processed.erase(e.first.first);
      	  processed.erase(e.first.second);
        }
      }
      if (processed.size() == vertex_set.size()) { 
        answer.insert(current);
      } 
    }; 
 
    for (auto const &edge : weight_edges) { 
      if (processed.count(edge.first.first) == 0 && processed.count(edge.first.second) == 0) { 	
        current.insert(edge.first); 
        processed.insert({edge.first.first, edge.first.second});
        processed.erase(Infty);
        findLambda();  
        current.erase(edge.first);
        processed.erase(edge.first.first);
        processed.erase(edge.first.second);
      }
    }
  } 

  std::multimap<size_t, arc_t> minimal_matching; 
  size_t min = std::numeric_limits<size_t>::max() / 4;
  for (auto const & match: answer) { 
    size_t weight = 0; 
    //auto max_edge = *match.cbegin(); 
    for (auto const &edge: match) { 
      weight += weight_edges.find(edge)->second.first;
      //if (weight_edges.find(max_edge)->second.first < weight_edges.find(edge)->second.first) { 
        //max_edge = edge;
      //} 
    } 
    if ((weight/* - weight_edges.find(max_edge)->second.first*/) < min) { 
      min = weight; 
      minimal_matching.clear();
      for (auto const &edge: match) { 
        minimal_matching.insert(std::make_pair(weight_edges.find(edge)->second.first, edge));
      } 
      //minimal_matching.erase(max_edge);
    } 
  } 
 
  /*std::cerr << "Find minimal matching " << minimal_matching.size() << std::endl;
  for (auto const & edge : minimal_matching) {
    std::cerr << "(" << edge.second.first  << ", " << edge.second.second << ") " << edge.first << std::endl; 
  }*/

  return minimal_matching;
} 


template<class graph_t>
size_t Algorithm<graph_t>::take_edge_on_color(vertex_t const & x, mcolor_t const & color, vertex_t const & y) {
  size_t num_rear = 0;

  //std::cerr << "Start to work with (" << x << ", " << y << ") " << genome_match::mcolor_to_name(color) << std::endl;  

  if (y == Infty) {
    mcolor_t need_color(color, graph->get_edge_multicolor(x, y), mcolor_t::Difference);    
    mularcs_t mularcs_x = graph->get_adjacent_multiedges_with_info(x, false);
    mularcs_x.erase(y);
  
    for (auto const & arc : mularcs_x) {
      if (need_color.includes(arc.second) && graph->is_vec_T_consistent_color(arc.second)) {
      	graph->apply(twobreak_t(x, arc.first, Infty, Infty, arc.second));
      	++num_rear; 
      } 
    } 
  } else { 
    mcolor_t need_color(color, graph->get_edge_multicolor(x, y), mcolor_t::Difference);
    //std::cerr << "Target color " << genome_match::mcolor_to_name(need_color) << " " << graph->is_vec_T_consistent_color(need_color) << std::endl; 
    //std::cerr << " edge: " << genome_match::mcolor_to_name(graph->get_edge_multicolor(x, y)) << std::endl;
  
    mularcs_t mularcs_x = graph->get_adjacent_multiedges_with_info(x, false);
    mularcs_x.erase(y);
  
    mularcs_t mularcs_y = graph->get_adjacent_multiedges_with_info(y, false);
    mularcs_y.erase(x);

    /*std::cerr << "mulacrs_x vertex have " << std::endl; 
    for (auto const & l : mularcs_x) { 
      std::cerr << genome_match::mcolor_to_name(l.second) << " " << graph->is_vec_T_consistent_color(l.second) << " " << l.first << std::endl;
    }

    std::cerr << "mularcs_y vertex have " << std::endl; 
    for (auto const & r : mularcs_y) {
      std::cerr << genome_match::mcolor_to_name(r.second) << " " << graph->is_vec_T_consistent_color(r.second) << " " << r.first << std::endl;
    }*/
  

    typedef std::pair<std::pair<vertex_t, structure::Mcolor>, size_t> colacr_t;
    utility::equivalence<colacr_t> equiv; 

    for (auto const & arc_x : mularcs_x) { 
      for (auto const & arc_y : mularcs_y) {
        if (color.includes(arc_x.second) && color.includes(arc_y.second)) {  
          mcolor_t color(arc_x.second, arc_y.second, mcolor_t::Intersection);
          if (color.size() > 0) {   
            equiv.addrel(std::make_pair(arc_x, 0), std::make_pair(arc_y, 1));
          }
        } 
      } 
    }  
    equiv.update();

    std::map<colacr_t, std::set<colacr_t> > const & classes = equiv.get_eclasses<std::set<colacr_t> >();   
    for (auto const & color_set : classes) { 
      std::multimap<mcolor_t, vertex_t> left; 
      std::multimap<mcolor_t, vertex_t> right; 
      for (auto const & color : color_set.second) { 
        if (color.second == 0) {
          left.insert(std::make_pair(color.first.second, color.first.first));
        } else {
          right.insert(std::make_pair(color.first.second, color.first.first));
        } 
      }
	
      /*
      std::cerr << "left vertex have " << std::endl; 
      for (auto const & l : left) { 
        std::cerr << genome_match::mcolor_to_name(l.first) << " " << graph->is_vec_T_consistent_color(l.first) << " " << l.second << std::endl;
      }

      std::cerr << "right vertex have " << std::endl; 
      for (auto const & r : right) {
        std::cerr << genome_match::mcolor_to_name(r.first) << " " << graph->is_vec_T_consistent_color(r.first) << " " << r.second << std::endl;
      }
      */

      if (left.size() == 1 && right.size() == 1 && left.begin()->first == right.begin()->first) {
        assert(graph->is_vec_T_consistent_color(left.begin()->first)); 
        assert(graph->is_vec_T_consistent_color(right.begin()->first)); 
        if (need_color.includes(left.begin()->first)) {
          //std::cerr << "Good 2-break: " << x << " " << left.begin()->second << " " << y << " " << right.begin()->second <<  " " << genome_match::mcolor_to_name(left.begin()->first) << std::endl;
      	  graph->apply(twobreak_t(x, left.begin()->second, y, right.begin()->second, left.begin()->first));
      	  ++num_rear; 
        } 
      } else if (left.size() == 1 || right.size() == 1) {  
        if (left.size() == 1 && need_color.includes(left.begin()->first)) {
          //std::cerr << "Left and go recursevly" << std::endl;
          num_rear += take_edge_on_color(y, left.begin()->first, right.begin()->second);
          assert(graph->is_vec_T_consistent_color(left.begin()->first)); 
          //std::cerr << "Left 2-break: " << x << " " << left.begin()->second << " " << y << " " << right.begin()->second <<  " " << genome_match::mcolor_to_name(left.begin()->first) << std::endl;
          graph->apply(twobreak_t(x, left.begin()->second, y, right.begin()->second, left.begin()->first));
          ++num_rear; 
        } else if (right.size() == 1 && need_color.includes(right.begin()->first)) {
          //std::cerr << "Right and go recursevly" << std::endl;
          num_rear += take_edge_on_color(x, right.begin()->first, left.begin()->second);
      	  assert(graph->is_vec_T_consistent_color(right.begin()->first));  
          //std::cerr << "Left 2-break: " << x << " " << left.begin()->second << " " << y << " " << right.begin()->second <<  " " << genome_match::mcolor_to_name(right.begin()->first) << std::endl;
          graph->apply(twobreak_t(x, left.begin()->second, y, right.begin()->second, right.begin()->first));
          ++num_rear; 
        }   
      } else { 
        assert(left.size() == 1 || right.size() == 1); 
      } 
    } 
  } 

  return num_rear;
} 


template<class graph_t>
bool Algorithm<graph_t>::stage6() { 
  size_t number_break = 0;

  utility::equivalence<vertex_t> components = graph->split_on_components();
  std::map<vertex_t, std::set<vertex_t> > const & classes = components.get_eclasses<std::set<vertex_t> >();
  for (auto const & vertex_set : classes) { 
 
    //std::cerr << "Component have size " << vertex_set.second.size() << std::endl;

    if (vertex_set.second.size() <= max_size_component) { 

      std::cerr << "Start process component" << std::endl;
      for (auto const v : vertex_set.second) {
      	std::cerr << v << " "; 
      } 
      std::cerr << std::endl;

      std::multimap<size_t, arc_t> matching = create_minimal_matching(vertex_set.second);
      /*for (auto const & edge : matching) { 
        std::cerr << "(" << edge.second.first << ", " << edge.second.second << "):" << edge.first << std::endl;
      }*/

      //std::cerr << "Start tested matching" << std::endl;
      //std::multimap<size_t, arc_t> matching = wrapper_create_minimal_matching(vertex_set.second);
      /*for (auto const & edge : matching) { 
        std::cerr << "(" << edge.second.first << ", " << edge.second.second << "):" << edge.first << std::endl;
      }*/
      
      //std::cerr << "Start process minimal matching" << std::endl;
      assert(!matching.empty());
      auto edge = matching.cbegin()->second;
      for (auto weight_edge: matching) { 
        if (weight_edge.second.first != Infty && weight_edge.second.second != Infty) { 
          edge = weight_edge.second;
          break;
        }
      } 

      number_break += take_edge_on_color(edge.first, graph->get_complete_color(), edge.second);
    } 
  }  

  return (number_break != 0);
} 

#endif
