#ifndef TOTAL_BRUTE_FORCE_HPP
#define TOTAL_BRUTE_FORCE_HPP

#include "blossom5/PerfectMatching.h"

template<class graph_t>
struct Algorithm<graph_t>::BruteForceTwo : public Algorithm<graph_t>::Stage {
  typedef Stage base;

  typedef std::pair<int, int> int_edge_t;
  typedef std::tuple<int, int, int> weight_edge_t;
  typedef std::pair<vertex_t, vertex_t> edge_t;
  
  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::twobreak_t twobreak_t;

  BruteForceTwo(std::shared_ptr<graph_t> const & graph)
  : Stage(graph) 
  , max_number(0)
  , super_infinity("ooo")
  { 
  }

  bool do_action() override;

  std::string get_name() override {
    return "Total resolve component with score function."; 
  }

  void init_process_tc_edges(std::set<vertex_t> const & vertex_set);

private:
  size_t process_component(std::set<vertex_t> & vertex_set);
  
  std::vector<weight_edge_t> init(std::set<vertex_t> & vertex_set);
  
  std::list<twobreak_t> take_edge_on_color(vertex_t const & x, mcolor_t const & color, vertex_t const & y);

  void update_process_tc_edges(edge_t const & central, std::list<twobreak_t> const & twobreaks);

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

  std::multimap<int, edge_t> queue_score_edges;
  std::multimap<int, edge_t> queue_tc_edges; 
  std::map<edge_t, int> tc_edges;
  std::map<edge_t, int> score_edges;

  vertex_t super_infinity;
  std::unordered_map<int, vertex_t> num2vertex;
  std::unordered_map<vertex_t, int> vertex2num;  

  std::unordered_set<int> pseudo_infinity_verteces;  
  std::unordered_set<int> double_pseudo_infinity_verteces;  
};

template<class graph_t>
bool Algorithm<graph_t>::BruteForceTwo::do_action() {
  size_t number_rear = 0; 

  utility::equivalence<vertex_t> components = this->graph->split_on_components();
  std::map<vertex_t, std::set<vertex_t> > const & classes = components.get_eclasses<std::set<vertex_t> >();

  for (auto const & vertex_set : classes) { 
    std::set<vertex_t> component = vertex_set.second; 
    //std::cerr << "Start work with new component " << std::endl;
    number_rear += process_component(component);
  }

  return (number_rear != 0);
} 

template<class graph_t>
size_t Algorithm<graph_t>::BruteForceTwo::process_component(std::set<vertex_t> & vertex_set) {
  tc_edges.clear(); 
  queue_tc_edges.clear();
  size_t number_rear = 0;

  init_process_tc_edges(vertex_set);

  edge_t edge; 
    
  auto main_stage_lambda = [&]() -> void {
    std::list<twobreak_t> twobreaks = take_edge_on_color(edge.first, this->graph->get_complete_color(), edge.second);      
    number_rear += twobreaks.size(); 

    vertex_set.erase(edge.first);
    vertex_set.erase(edge.second);

    //update_process_tc_edges(edge, twobreaks);        
  };

    //std::cerr << "Start process tc edges" << std::endl; 
  if (!tc_edges.empty()) {
    auto iter = queue_tc_edges.begin();
    edge = iter->second;

    tc_edges.erase(edge);
    queue_tc_edges.erase(iter);

    main_stage_lambda();    
  } else {
    auto graph_edges = init(vertex_set);

    PerfectMatching pm(num2vertex.size(), graph_edges.size());
    for (auto const & edge : graph_edges) { 
      int x = 0; 
      int y = 0; 
      int score = 0;
      std::tie(x, y, score) = edge;

      //std::cerr << x << " " << y << std::endl;
      //std::cerr << num2vertex.find(x)->second << " " << num2vertex.find(y)->second << std::endl; 

      pm.AddEdge(x, y, score);
    }
    pm.Solve();

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

    edge = get_min_edge_in_matching_lambda();

    main_stage_lambda();
  } 

  //std::cerr << "End all " << std::endl;

  return number_rear;
}

template<class graph_t>
void Algorithm<graph_t>::BruteForceTwo::init_process_tc_edges(std::set<vertex_t> const & vertex_set) {
  std::unordered_set<vertex_t> processed;

  for (vertex_t const & v : vertex_set) {
    mularcs_t mularcs = this->graph->get_all_adjacent_multiedges(v);

    for (auto const & arc : mularcs) { 
      if ((processed.count(arc.first) == 0) 
          && (this->graph->is_contain_T_consistent_color(v, arc.first))) { //|| this->graph->is_postponed_deletion(v, arc.first))) { 
        int score = this->graph->calculate_cost(v, arc.first);
        queue_tc_edges.insert(std::make_pair(score, edge_t(v, arc.first)));
        tc_edges.insert(std::make_pair(edge_t(v, arc.first), score));
      } 
    }
    processed.insert(v);
  }
}

template<class graph_t>
std::vector<std::tuple<int, int, int> > Algorithm<graph_t>::BruteForceTwo::init(std::set<vertex_t> & vertex_set) {
  max_number = 0; 
  tc_edges.clear();
  queue_tc_edges.clear();
  vertex2num.clear();
  num2vertex.clear(); 
  pseudo_infinity_verteces.clear();
  double_pseudo_infinity_verteces.clear();

  std::vector<weight_edge_t> graph_edges;

/*
  for (auto iter = vertex_set.begin(); iter != vertex_set.end();) {
    mularcs_t mularcs = this->graph->get_all_adjacent_multiedges(*iter);
    if (mularcs.size() == 1 && this->graph->get_complete_color() == mularcs.begin()->second) { 
      vertex_set.erase(iter++);
    } else { 
      ++iter;
    }
  } 
*/

  std::unordered_set<vertex_t> processed;
   //Create verteces and edges for graph
  for (vertex_t const & v : vertex_set) {
    mularcs_t mularcs = this->graph->get_all_adjacent_multiedges(v);

    for (auto const & arc : mularcs) { 
      if (processed.count(arc.first) == 0) { 
        //Prepare verteces
        int num_x = vertex_to_number(v);  
        int num_y = -1; 
        if (arc.first == Infty) { 
          vertex_t pseudo_vertex = "o" + std::to_string(pseudo_infinity_verteces.size()) + "o";
          num_y = vertex_to_number(pseudo_vertex);
          this->pseudo_infinity_verteces.insert(num_y);
          vertex_t double_pseudo_vertex = "oo" + std::to_string(double_pseudo_infinity_verteces.size()) + "oo";
          int num_z = vertex_to_number(double_pseudo_vertex);
          this->double_pseudo_infinity_verteces.insert(num_z);

          graph_edges.push_back(std::make_tuple(num_y, num_z, 0));   
        } else { 
          num_y = vertex_to_number(arc.first); 
        }        
        assert(num_y != -1);
        int score = this->graph->calculate_cost(v, arc.first);
        score_edges.insert(std::make_pair(edge_t(v, arc.first), score));
        graph_edges.push_back(std::make_tuple(num_x, num_y, score)); 
      } 
    }
    processed.insert(v);
  }

  if (vertex2num.size() % 2 != 0 && !pseudo_infinity_verteces.empty()) { 
    int num_y = vertex_to_number(super_infinity);

    for (int num_x : double_pseudo_infinity_verteces) { 
      graph_edges.push_back(std::make_tuple(num_x, num_y, 0)); 
    } 
  }

  //std::cerr << "Component size " << vertex_set.size() << " " << vertex2num.size() << " " << pseudo_infinity_verteces.size() << std::endl;  
  assert((vertex2num.size() % 2 == 0) && (max_number % 2 == 0) && (num2vertex.size() % 2 == 0));

  return graph_edges;
}

template<class graph_t>
void Algorithm<graph_t>::BruteForceTwo::update_process_tc_edges(edge_t const & central, std::list<twobreak_t> const & twobreaks) {
  std::unordered_set<vertex_t> viewed_vertex;
  
  for (twobreak_t const & twobreak : twobreaks) { 
    viewed_vertex.insert({twobreak.get_vertex(0), twobreak.get_vertex(1), twobreak.get_vertex(2), twobreak.get_vertex(3)});
  }

  viewed_vertex.erase(Infty);
  viewed_vertex.erase(central.first);
  viewed_vertex.erase(central.second);

  std::unordered_set<vertex_t> processed;

  for (vertex_t const & viewed : viewed_vertex) { 
    mularcs_t mularcs = this->graph->get_all_adjacent_multiedges(viewed);

    for (auto const & arc : mularcs) { 
      if (processed.count(arc.first) == 0 
        && (this->graph->is_contain_T_consistent_color(viewed, arc.first))) { // || this->graph->is_postponed_deletion(viewed, arc.first))) { 
        auto iter = tc_edges.find(edge_t(viewed, arc.first)); 
        if (iter == tc_edges.end()) {
          iter = tc_edges.find(edge_t(arc.first, viewed));
        } 

        int new_score = this->graph->calculate_cost(viewed, arc.first);
        if (iter == tc_edges.end()) {
          queue_tc_edges.insert(std::make_pair(new_score, edge_t(viewed, arc.first)));
          tc_edges.insert(std::make_pair(edge_t(viewed, arc.first), new_score));        
        } else { 
          if (new_score != iter->second) { 
            auto range = queue_tc_edges.equal_range(iter->second);
            
            for (auto it = range.first; it != range.second; ++it) { 
              if (it->second == iter->first) { 
                queue_tc_edges.erase(it);
                break;
              }
            }

            queue_tc_edges.insert(std::make_pair(new_score, iter->first));
            iter->second = new_score;
          } 
        }
      }
    }
    processed.insert(viewed);
  }  
} 


template<class graph_t>
std::list<typename graph_t::twobreak_t> Algorithm<graph_t>::BruteForceTwo::take_edge_on_color(vertex_t const & x, mcolor_t const & color, vertex_t const & y) {
  std::list<twobreak_t> result;
  //std::cerr << "Start to work with (" << x << ", " << y << "): " << genome_match::mcolor_to_name(color) << std::endl;   
  //size_t num_rear = 0; 

  if (y == Infty) {
    mcolor_t need_color(color, this->graph->get_all_multicolor_edge(x, y), mcolor_t::Difference);    
    mularcs_t mularcs_x = this->graph->get_all_adjacent_multiedges_with_info(x, false);
    mularcs_x.erase(y);
  
    for (auto const & arc : mularcs_x) {
      if (need_color.includes(arc.second) && this->graph->is_vec_T_consistent_color(arc.second)) {
        twobreak_t twobreak(x, arc.first, Infty, Infty, arc.second);
        this->graph->apply(twobreak);
        //++num_rear;
        result.push_back(twobreak);
      } 
    }
  } else { 
    mcolor_t need_color(color, this->graph->get_all_multicolor_edge(x, y), mcolor_t::Difference);
    
    //std::cerr << "Target color " << genome_match::mcolor_to_name(need_color) << " " << this->graph->is_vec_T_consistent_color(need_color) << std::endl; 
    //std::cerr << " edge: " << genome_match::mcolor_to_name(this->graph->get_all_multicolor_edge(x, y)) << std::endl;
  
    mularcs_t mularcs_x = this->graph->get_all_adjacent_multiedges_with_info(x, false);
    mularcs_x.erase(y);
  
    mularcs_t mularcs_y = this->graph->get_all_adjacent_multiedges_with_info(y, false);
    mularcs_y.erase(x);

/*
    std::cerr << "mulacrs_x vertex have " << std::endl; 
    for (auto const & l : mularcs_x) { 
      std::cerr << genome_match::mcolor_to_name(l.second) << " " << this->graph->is_vec_T_consistent_color(l.second) << " " << l.first << std::endl;
    }

    std::cerr << "mularcs_y vertex have " << std::endl; 
    for (auto const & r : mularcs_y) {
      std::cerr << genome_match::mcolor_to_name(r.second) << " " << this->graph->is_vec_T_consistent_color(r.second) << " " << r.first << std::endl;
    }  
*/
    typedef std::pair<std::pair<vertex_t, structure::Mcolor>, size_t> colacr_t;
    utility::equivalence<colacr_t> equiv; 

    for (auto const & arc_x : mularcs_x) { 
      for (auto const & arc_y : mularcs_y) {
        if (need_color.includes(arc_x.second) && need_color.includes(arc_y.second)) {  
          mcolor_t color(arc_x.second, arc_y.second, mcolor_t::Intersection);
          if (color.size() > 0) {   
            equiv.addrel(std::make_pair(arc_x, 0), std::make_pair(arc_y, 1));
          }
        }
      } 
    }  
    equiv.update();
    
    std::map<colacr_t, std::set<colacr_t> > const & classes = equiv.get_eclasses<std::set<colacr_t> >();   
   
    mcolor_t temp_need_color = need_color;
    for (auto const & color_set : classes) {
      mularcs_t left; 
      mularcs_t right; 

      for (auto const & color : color_set.second) { 
        if (color.second == 0) {
          left.insert(color.first.first, color.first.second);
        } else {
          right.insert(color.first.first, color.first.second);
        }
      }

      if (left.size() == 1) {
        temp_need_color = mcolor_t(temp_need_color, left.begin()->second, mcolor_t::Difference); 
      } else if (right.size() == 1) { 
        temp_need_color = mcolor_t(temp_need_color, right.begin()->second, mcolor_t::Difference); 
      }
    } 

    if (!temp_need_color.empty()) { 
      //std::cerr << "Non-empty need color " << genome_match::mcolor_to_name(temp_need_color) << std::endl;
      return result; 
    }

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
        std::cerr << genome_match::mcolor_to_name(l.first) << " " << this->graph->is_vec_T_consistent_color(l.first) << " " << l.second << std::endl;
      }

      std::cerr << "right vertex have " << std::endl; 
      for (auto const & r : right) {
        std::cerr << genome_match::mcolor_to_name(r.first) << " " << this->graph->is_vec_T_consistent_color(r.first) << " " << r.second << std::endl;
      }
    */  

      if (left.size() == 1 && right.size() == 1 && left.begin()->first == right.begin()->first) {
        assert(this->graph->is_vec_T_consistent_color(left.begin()->first)); 
        assert(this->graph->is_vec_T_consistent_color(right.begin()->first)); 
        if (need_color.includes(left.begin()->first)) {
          //std::cerr << "Good 2-break: " << x << " " << left.begin()->second << " " << y << " " << right.begin()->second <<  " " << genome_match::mcolor_to_name(left.begin()->first) << std::endl;
          twobreak_t twobreak(x, left.begin()->second, y, right.begin()->second, left.begin()->first);
          this->graph->apply(twobreak);
          //++num_rear;
          result.push_back(twobreak);
        } 
      } else if (left.size() == 1 || right.size() == 1) {  
        if (left.size() == 1 && need_color.includes(left.begin()->first)) {
          //std::cerr << "Left and go recursevly right" << std::endl;
          assert(this->graph->is_vec_T_consistent_color(left.begin()->first)); 
          
          std::list<twobreak_t> local_twobreaks; 
          //size_t local_twobreaks = 0; 
          auto iter_right = right.begin();
          while ((local_twobreaks.empty())/*(local_twobreaks == 0)*/ && (iter_right != right.end())) { 
            //std::cerr << "See on " << iter_right->second << std::endl;
            local_twobreaks = take_edge_on_color(y, left.begin()->first, iter_right->second);
            if (local_twobreaks.empty()) { //} == 0) {
              ++iter_right;
            } 
          } 
          
          if (iter_right != right.end()) {
            //std::cerr << "Left 2-break: " << x << " " << left.begin()->second << " " << y << " " << iter_right->second <<  " " << genome_match::mcolor_to_name(left.begin()->first) << std::endl;
            twobreak_t twobreak(x, left.begin()->second, y, iter_right->second, left.begin()->first);
            this->graph->apply(twobreak);
           // num_rear += (local_twobreaks + 1);
            result.insert(result.end(), local_twobreaks.begin(), local_twobreaks.end());
            result.push_back(twobreak);
          }
          
        } else if (right.size() == 1 && need_color.includes(right.begin()->first)) {
          //std::cerr << "Right and go recursevly left" << std::endl;
          assert(this->graph->is_vec_T_consistent_color(right.begin()->first));  
          
          std::list<twobreak_t> local_twobreaks;
          //size_t local_twobreaks = 0; 
          auto iter_left = left.begin();
          while ((local_twobreaks.empty())/*(local_twobreaks == 0)*/ && (iter_left != left.end())) {
            //std::cerr << "See on " << iter_left->second << std::endl; 
            local_twobreaks = take_edge_on_color(x, right.begin()->first, iter_left->second);
            if (local_twobreaks.empty()) { // == 0) {
              ++iter_left;
            } 
          }

          if (iter_left != left.end()) {
            //std::cerr << "Right 2-break: " << x << " " << iter_left->second << " " << y << " " << right.begin()->second <<  " " << genome_match::mcolor_to_name(right.begin()->first) << std::endl;
            twobreak_t twobreak(x, iter_left->second, y, right.begin()->second, right.begin()->first);
            this->graph->apply(twobreak);
            //num_rear += (local_twobreaks + 1);
            result.insert(result.end(), local_twobreaks.begin(), local_twobreaks.end());
            result.push_back(twobreak);
          }             
        } else { 
          ;//std::cerr << "Se are here " << std::endl;
        }   
      } else { 
        assert(left.size() == 1 || right.size() == 1); 
      } 
    } 
  } 

  return result; //num_rear;
}

#endif