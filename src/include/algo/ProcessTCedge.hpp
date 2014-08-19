#ifndef PROCESS_TC_EDGES_HPP
#define PROCESS_TC_EDGES_HPP

template<class graph_t>
struct Algorithm<graph_t>::ProcessTCedges : public Algorithm<graph_t>::Stage {
  typedef Stage base;
  
  typedef std::pair<vertex_t, vertex_t> edge_t;
  
  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::twobreak_t twobreak_t;

  ProcessTCedges(std::shared_ptr<graph_t> const & graph)
  : Stage(graph) 
  { 
  }

  bool do_action() override;

  std::string get_name() override {
    return "Resolve all T-consistent edges."; 
  }

private:
  size_t process_tc_edges(std::set<vertex_t> const & vertex_set);
  
  void init_process_tc_edges(std::set<vertex_t> const & vertex_set);
  
  void update_process_tc_edges(edge_t const & central, std::list<twobreak_t> const & twobreaks);

private: 
  //go to graph
  std::list<twobreak_t> take_edge_on_color(vertex_t const & x, mcolor_t const & color, vertex_t const & y);

private:
  std::multimap<int, edge_t> queue_tc_edges; 
  std::map<edge_t, int> tc_edges;
};

template<class graph_t>
bool Algorithm<graph_t>::ProcessTCedges::do_action() {
  size_t number_rear = 0; 

  utility::equivalence<vertex_t> components = this->graph->split_on_components();
  std::map<vertex_t, std::set<vertex_t> > const & classes = components.get_eclasses<std::set<vertex_t> >();

  for (auto const & vertex_set : classes) { 
    number_rear += process_tc_edges(vertex_set.second);
  }

  //FIXME remove
  std::unordered_set<vertex_t> processed;
  for (vertex_t const & v : *this->graph) { 
    mularcs_t mularcs = this->graph->get_all_adjacent_multiedges(v);

    for (auto const & arc : mularcs) { 
      if (processed.count(arc.first) == 0) { 
        if (this->graph->is_contain_T_consistent_color(v, arc.first)) {  
          std::cerr << v << " " << arc.first << std::endl;
        } 
      } 
    }
    processed.insert(v);
  }

  return (number_rear != 0);
} 

template<class graph_t>
size_t Algorithm<graph_t>::ProcessTCedges::process_tc_edges(std::set<vertex_t> const & vertex_set) {
  size_t number_rear = 0;
  
  init_process_tc_edges(vertex_set); 

  while (!tc_edges.empty()) { 
    auto iter = queue_tc_edges.begin();
    auto edge = iter->second;
    
    std::list<twobreak_t> twobreaks = take_edge_on_color(edge.first, this->graph->get_complete_color(), edge.second);      
    number_rear += twobreaks.size(); 

    tc_edges.erase(edge);
    queue_tc_edges.erase(iter);

    update_process_tc_edges(edge, twobreaks);
  }

  return number_rear;
}

template<class graph_t>
void Algorithm<graph_t>::ProcessTCedges::init_process_tc_edges(std::set<vertex_t> const & vertex_set) {
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
void Algorithm<graph_t>::ProcessTCedges::update_process_tc_edges(edge_t const & central, std::list<twobreak_t> const & twobreaks) {
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
std::list<typename graph_t::twobreak_t> Algorithm<graph_t>::ProcessTCedges::take_edge_on_color(vertex_t const & x, mcolor_t const & color, vertex_t const & y) {
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
  
    if (need_color.empty()) { 
      return result;
    }

    mularcs_t mularcs_x = this->graph->get_all_adjacent_multiedges_with_info(x, false);
    mularcs_x.erase(y);
  
    mularcs_t mularcs_y = this->graph->get_all_adjacent_multiedges_with_info(y, false);
    mularcs_y.erase(x);

    /*std::cerr << "mulacrs_x vertex have " << std::endl; 
    for (auto const & l : mularcs_x) { 
      std::cerr << genome_match::mcolor_to_name(l.second) << " " << this->graph->is_vec_T_consistent_color(l.second) << " " << l.first << std::endl;
    }

    std::cerr << "mularcs_y vertex have " << std::endl; 
    for (auto const & r : mularcs_y) {
      std::cerr << genome_match::mcolor_to_name(r.second) << " " << this->graph->is_vec_T_consistent_color(r.second) << " " << r.first << std::endl;
    } */ 

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
    
      /*std::cerr << "left vertex have " << std::endl; 
      for (auto const & l : left) { 
        std::cerr << genome_match::mcolor_to_name(l.first) << " " << this->graph->is_vec_T_consistent_color(l.first) << " " << l.second << std::endl;
      }

      std::cerr << "right vertex have " << std::endl; 
      for (auto const & r : right) {
        std::cerr << genome_match::mcolor_to_name(r.first) << " " << this->graph->is_vec_T_consistent_color(r.first) << " " << r.second << std::endl;
      }*/
      

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