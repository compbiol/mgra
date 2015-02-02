#ifndef BRUTE_FORCE_HPP 
#define BRUTE_FORCE_HPP

namespace algo { 

template<class graph_pack_t>
struct BruteForce : public algo::AbsStage<graph_pack_t> {
  using mcolor_t = typename graph_pack_t::mcolor_type;
  
  using edge_t = typename graph_pack_t::edge_t;  
  using arc_t = typename graph_pack_t::arc_t; 
  using mularcs_t = typename graph_pack_t::mularcs_t; 
  using twobreak_t = typename graph_pack_t::twobreak_t;

  BruteForce(size_t size_component)
  : AbsStage<graph_pack_t>("Resolve component with score function", "brute_force", 3) 
  , max_size_component(size_component)
  {
  }
    
  bool run(graph_pack_t & graph_pack) override;
  
private:
  std::multimap<size_t, edge_t> create_minimal_matching(graph_pack_t const & graph_pack, 
                                                        std::set<vertex_t> const & vertex_set); 
  
private: 
  size_t max_size_component;

private:
  DECL_LOGGER("BruteForce");
};

template<class graph_pack_t>
bool BruteForce<graph_pack_t>::run(graph_pack_t & graph_pack) { 
  size_t number_rear = 0; // number of rearrangements 

  utility::equivalence<vertex_t> components = graph_pack.split_on_components();
  std::map<vertex_t, std::set<vertex_t> > const & classes = components.get_eclasses<std::set<vertex_t> >();

  for (auto const & vertex_set : classes) { 
     if ((vertex_set.second.size() <= max_size_component)) { 
      std::multimap<size_t, edge_t> matching;
      matching = create_minimal_matching(graph_pack, vertex_set.second);
      auto edge = matching.cbegin()->second;
      auto result = graph_pack.take_edge_on_color(edge.first, graph_pack.multicolors.get_complete_color(), edge.second);

      for (twobreak_t const & twobreak : result) { 
        graph_pack.apply(twobreak);
        ++number_rear;
      }      
    } 
  }

  return (number_rear != 0);
} 

template<class graph_pack_t>
std::multimap<size_t, typename graph_pack_t::edge_t> 
      BruteForce<graph_pack_t>::create_minimal_matching(graph_pack_t const & graph_pack, 
                                                  std::set<vertex_t> const & vertex_set) {
  std::map<edge_t, std::pair<size_t, mcolor_t> > weight_edges; 

  for (vertex_t const & v : vertex_set) { 
    mularcs_t mularcs = graph_pack.get_all_adjacent_multiedges(v);
    for (arc_t const & arc: mularcs) {
      if (weight_edges.count(std::make_pair(v, arc.first)) == 0 && weight_edges.count(std::make_pair(arc.first, v)) == 0) { 
      	weight_edges.insert(std::make_pair(std::make_pair(v, arc.first), std::make_pair(graph_pack.calculate_cost(v, arc.first), arc.second)));
      }
    } 
  } 

  std::unordered_set<vertex_t> processed;
  std::set<edge_t> current; 
  for (auto const & edge : weight_edges) { 
    bool is_have_tc = false; 
    auto edge_colors = graph_pack.split_color(edge.second.second);
    
    for (auto const & color : edge_colors) { 
      if (graph_pack.multicolors.is_T_consistent_color(color) && !graph_pack.multicolors.is_vec_T_consistent_color(color)) { 
        is_have_tc = true; 
        break;
      }
    }

    if (graph_pack.is_postponed_deletion(edge.first.first, edge.first.second) || is_have_tc) { 
      current.insert(edge.first); 
      processed.insert(edge.first.first); processed.insert(edge.first.second);
      processed.erase(Infty);
    } 
  } 

  std::set<std::set<edge_t> > answer;
  if (processed.size() == vertex_set.size()) { 
    answer.insert(current);
  } else { 
    std::function<void()> findLambda = [&] () -> void { 
      for (auto const &e : weight_edges) { 
        if (processed.count(e.first.first) == 0 && processed.count(e.first.second) == 0) { 
          processed.insert(e.first.first); processed.insert(e.first.second);
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
        processed.insert(edge.first.first); processed.insert(edge.first.second);
        processed.erase(Infty);
        findLambda();  
        current.erase(edge.first);
        processed.erase(edge.first.first);
        processed.erase(edge.first.second);
      }
    }
  } 

  std::multimap<size_t, edge_t> minimal_matching; 
  size_t min = std::numeric_limits<size_t>::max() / 4;
  for (auto const & match: answer) { 
    size_t weight = 0; 
    for (auto const & edge: match) { 
      weight += weight_edges.find(edge)->second.first;
    } 

    if (weight < min) { 
      min = weight; 
      minimal_matching.clear();
      for (auto const &edge: match) { 
        minimal_matching.insert(std::make_pair(weight_edges.find(edge)->second.first, edge));
      } 
    }

  } 
 
  /*
  TRACE("Find minimal matching");
  for (auto const & edge : minimal_matching) {
    TRACE("(" << edge.second.first  << ", " << edge.second.second << ") " << edge.first); 
  }
  */

  return minimal_matching;
} 

} 


/*template<class graph_t>
size_t Algorithm<graph_t>::BruteForce::take_edge_on_color(vertex_t const & x, mcolor_t const & color, vertex_t const & y) {
  size_t num_rear = 0;

  //std::cerr << "Start to work with (" << x << ", " << y << ") " << genome_match::mcolor_to_name(color) << std::endl;  
  std::list<twobreak_t> result; 

  if (y == Infty) {
    mcolor_t need_color(color, graph_pack.get_all_multicolor_edge(x, y), mcolor_t::Difference);    
    mularcs_t mularcs_x = this->graph->get_all_adjacent_multiedges_with_info(x, false);
    mularcs_x.erase(y);
  
    for (auto const & arc : mularcs_x) {
      if (need_color.includes(arc.second) && this->graph->multicolors.is_vec_T_consistent_color(arc.second)) {
        twobreak_t twobreak(x, arc.first, Infty, Infty, arc.second);
        this->graph->apply(twobreak);
        result.push_back(twobreak);
        ++num_rear; 
      } 
    } 
  } else { 
    mcolor_t need_color(color, this->graph->get_all_multicolor_edge(x, y), mcolor_t::Difference);
  
    mularcs_t mularcs_x = this->graph->get_all_adjacent_multiedges_with_info(x, false);
    mularcs_x.erase(y);
  
    mularcs_t mularcs_y = this->graph->get_all_adjacent_multiedges_with_info(y, false);
    mularcs_y.erase(x);

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
  
      
      if (left.size() == 1 && right.size() == 1 && left.begin()->first == right.begin()->first) {
        assert(this->graph->multicolors.is_vec_T_consistent_color(left.begin()->first)); 
        assert(this->graph->multicolors.is_vec_T_consistent_color(right.begin()->first)); 
        if (need_color.includes(left.begin()->first)) {
          //std::cerr << "Good 2-break: " << x << " " << left.begin()->second << " " << y << " " << right.begin()->second <<  " " << genome_match::mcolor_to_name(left.begin()->first) << std::endl;
          this->graph->apply(twobreak_t(x, left.begin()->second, y, right.begin()->second, left.begin()->first));
          ++num_rear; 
        } 
      } else if (left.size() == 1 || right.size() == 1) {  
        if (left.size() == 1 && need_color.includes(left.begin()->first)) {
          //std::cerr << "Left and go recursevly right" << std::endl;
          num_rear += take_edge_on_color(y, left.begin()->first, right.begin()->second);
          assert(this->graph->multicolors.is_vec_T_consistent_color(left.begin()->first)); 
          //std::cerr << "Left 2-break: " << x << " " << left.begin()->second << " " << y << " " << right.begin()->second <<  " " << genome_match::mcolor_to_name(left.begin()->first) << std::endl;
          this->graph->apply(twobreak_t(x, left.begin()->second, y, right.begin()->second, left.begin()->first));
          ++num_rear; 
        } else if (right.size() == 1 && need_color.includes(right.begin()->first)) {
          //std::cerr << "Right and go recursevly left" << std::endl;
          num_rear += take_edge_on_color(x, right.begin()->first, left.begin()->second);
          assert(this->graph->multicolors.is_vec_T_consistent_color(right.begin()->first));  
          //std::cerr << "Right 2-break: " << x << " " << left.begin()->second << " " << y << " " << right.begin()->second <<  " " << genome_match::mcolor_to_name(right.begin()->first) << std::endl;
          this->graph->apply(twobreak_t(x, left.begin()->second, y, right.begin()->second, right.begin()->first));
          ++num_rear; 
        } else { 
          std::cerr << "Se are here " << std::endl;
        }   
      } else { 
        assert(left.size() == 1 || right.size() == 1); 
      } 
    } 
  } 

  return num_rear;
}
*/

#if 0 
template<class graph_t>
std::list<typename graph_t::twobreak_t> Algorithm<graph_t>::ProcessWithBlossomV::take_edge_on_color(vertex_t const & x, mcolor_t const & color, vertex_t const & y) {
  std::list<twobreak_t> result;
  //std::cerr << "Start to work with (" << x << ", " << y << "): " << genome_match::mcolor_to_name(color) << std::endl;   

  if (y == Infty) {
    mcolor_t need_color(color, graph_pack.get_all_multicolor_edge(x, y), mcolor_t::Difference);    
    mularcs_t mularcs_x = graph_pack.get_all_adjacent_multiedges_with_info(x, false);
    mularcs_x.erase(y);
  
    for (arc_t const & arc : mularcs_x) {
      if (need_color.includes(arc.second) && graph_pack.multicolors.is_vec_T_consistent_color(arc.second)) {
        twobreak_t twobreak(x, arc.first, Infty, Infty, arc.second);
        this->graph->apply(twobreak);
        result.push_back(twobreak);
      } 
    }
  } else { 
    mcolor_t need_color(color, graph_pack.get_all_multicolor_edge(x, y), mcolor_t::Difference);
    
    //std::cerr << "Target color " << genome_match::mcolor_to_name(need_color) << " " << this->graph->multicolors.is_vec_T_consistent_color(need_color) << std::endl; 
    //std::cerr << " edge: " << genome_match::mcolor_to_name(this->graph->get_all_multicolor_edge(x, y)) << std::endl;
  
    mularcs_t mularcs_x = graph_pack.get_all_adjacent_multiedges_with_info(x, false);
    mularcs_x.erase(y);
  
    mularcs_t mularcs_y = graph_pack.get_all_adjacent_multiedges_with_info(y, false);
    mularcs_y.erase(x);

/*
    std::cerr << "mulacrs_x vertex have " << std::endl; 
    for (auto const & l : mularcs_x) { 
      std::cerr << genome_match::mcolor_to_name(l.second) << " " << this->graph->multicolors.is_vec_T_consistent_color(l.second) << " " << l.first << std::endl;
    }

    std::cerr << "mularcs_y vertex have " << std::endl; 
    for (auto const & r : mularcs_y) {
      std::cerr << genome_match::mcolor_to_name(r.second) << " " << this->graph->multicolors.is_vec_T_consistent_color(r.second) << " " << r.first << std::endl;
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
        std::cerr << genome_match::mcolor_to_name(l.first) << " " << this->graph->multicolors.is_vec_T_consistent_color(l.first) << " " << l.second << std::endl;
      }

      std::cerr << "right vertex have " << std::endl; 
      for (auto const & r : right) {
        std::cerr << genome_match::mcolor_to_name(r.first) << " " << this->graph->multicolors.is_vec_T_consistent_color(r.first) << " " << r.second << std::endl;
      }
    */  

      if (left.size() == 1 && right.size() == 1 && left.begin()->first == right.begin()->first) {
        assert(this->graph->multicolors.is_vec_T_consistent_color(left.begin()->first)); 
        assert(this->graph->multicolors.is_vec_T_consistent_color(right.begin()->first)); 
        if (need_color.includes(left.begin()->first)) {
          //std::cerr << "Good 2-break: " << x << " " << left.begin()->second << " " << y << " " << right.begin()->second <<  " " << genome_match::mcolor_to_name(left.begin()->first) << std::endl;
          twobreak_t twobreak(x, left.begin()->second, y, right.begin()->second, left.begin()->first);
          this->graph->apply(twobreak);
          result.push_back(twobreak);
        } 
      } else if (left.size() == 1 || right.size() == 1) {  
        if (left.size() == 1 && need_color.includes(left.begin()->first)) {
          //std::cerr << "Left and go recursevly right" << std::endl;
          assert(this->graph->multicolors.is_vec_T_consistent_color(left.begin()->first)); 
          
          std::list<twobreak_t> local_twobreaks; 

          auto iter_right = right.begin();
          while ((local_twobreaks.empty()) && (iter_right != right.end())) { 
            //std::cerr << "See on " << iter_right->second << std::endl;
            local_twobreaks = take_edge_on_color(y, left.begin()->first, iter_right->second);
            if (local_twobreaks.empty()) { 
              ++iter_right;
            } 
          } 
          
          if (iter_right != right.end()) {
            //std::cerr << "Left 2-break: " << x << " " << left.begin()->second << " " << y << " " << iter_right->second <<  " " << genome_match::mcolor_to_name(left.begin()->first) << std::endl;
            twobreak_t twobreak(x, left.begin()->second, y, iter_right->second, left.begin()->first);
            this->graph->apply(twobreak);
            result.insert(result.end(), local_twobreaks.begin(), local_twobreaks.end());
            result.push_back(twobreak);
          }          
        } else if (right.size() == 1 && need_color.includes(right.begin()->first)) {
          //std::cerr << "Right and go recursevly left" << std::endl;
          assert(this->graph->multicolors.is_vec_T_consistent_color(right.begin()->first));  
          
          std::list<twobreak_t> local_twobreaks;
          auto iter_left = left.begin();
          while ((local_twobreaks.empty()) && (iter_left != left.end())) {
            //std::cerr << "See on " << iter_left->second << std::endl; 
            local_twobreaks = take_edge_on_color(x, right.begin()->first, iter_left->second);
            if (local_twobreaks.empty()) { 
              ++iter_left;
            } 
          }

          if (iter_left != left.end()) {
            //std::cerr << "Right 2-break: " << x << " " << iter_left->second << " " << y << " " << right.begin()->second <<  " " << genome_match::mcolor_to_name(right.begin()->first) << std::endl;
            twobreak_t twobreak(x, iter_left->second, y, right.begin()->second, right.begin()->first);
            this->graph->apply(twobreak);
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

  return result; 
}
#endif

#endif
