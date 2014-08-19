#ifndef SIMPLE_PATH_STAGE_HPP
#define SIMPLE_PATH_STAGE_HPP

template<class graph_t>
struct Algorithm<graph_t>::ProcessSimplePath : public Algorithm<graph_t>::Stage { 
  typedef std::list<vertex_t> path_t;
  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::twobreak_t twobreak_t;
  

  explicit ProcessSimplePath(std::shared_ptr<graph_t> const & graph)
  : Stage(graph) 
  {
  }

  bool do_action() override;
  
  std::string get_name() override { 
    return "Process good and simple paths.";
  }

private:
  size_t process_simple_path(path_t& path); 
};

template<class graph_t>
bool Algorithm<graph_t>::ProcessSimplePath::do_action() { 
  bool isChanged = false; 
  size_t number_rear = 0; // number of rearrangements 

  do {
    number_rear = 0; 
    for (vertex_t const & v: *(this->graph)) {  
      if (this->graph->is_simple_vertex(v)) { 
        path_t path({v});
        std::unordered_set<vertex_t> processed({v, Infty}); // we count oo as already processed

        auto const find_simple_path_lambda = [&] (vertex_t const & prev, vertex_t const & cur, bool is_next) -> vertex_t { 
          vertex_t previous = prev;
          vertex_t current = cur;
          bool is_continue = true; 

          while (is_continue) {
            is_continue = false; 

            if (is_next) { 
              path.push_front(current);
            } else { 
              path.push_back(current);
            } 

            if (processed.find(current) == processed.end() && !this->graph->is_duplication_vertex(current)) {     
              processed.insert(current);
              mcolor_t previous_color = this->graph->get_all_multicolor_edge(current, previous); 
              mularcs_t&& new_edges = this->graph->get_all_adjacent_multiedges(current);
              new_edges.erase(previous);
              mcolor_t next_color = new_edges.union_multicolors();
      
              if ((this->graph->degree_vertex(current) == 2) && this->graph->get_complement_color(previous_color) == next_color
                  && !this->graph->is_postponed_deletion(current, new_edges.cbegin()->first)) {

                auto const check_lambda = [&] (vertex_t const & v) -> bool {
                  bool flag = true; 
                  std::set<mcolor_t> const & colors = this->graph->get_all_multicolor_edge_with_info(current, v, false); 
                  for (auto color = colors.cbegin(); color != colors.cend() && flag; ++color) {
                    flag = this->graph->is_vec_T_consistent_color(*color);
                  }     
                  return flag;
                };

                bool nedge = check_lambda(new_edges.cbegin()->first); 
                bool pedge = check_lambda(previous);
                if (nedge || pedge) { 
                  previous = current;
                  current = new_edges.cbegin()->first; 
                  is_continue = true;
                }
              } 
            }
          } 
          return current;
        };

        //std::cerr << std::endl << "start " << v << std::endl; 
        //std::cerr << "go to " << current.cbegin()->first << genome_match::mcolor_to_name(current.cbegin()->second) << std::endl; 
        mularcs_t const & current = this->graph->get_all_adjacent_multiedges(v);
        vertex_t const & last = find_simple_path_lambda(v, current.cbegin()->first, true);
        if (last != v) { 
          //std::cerr << "go to " << current.crbegin()->first << genome_match::mcolor_to_name(current.crbegin()->second) << std::endl; 
          find_simple_path_lambda(v, current.crbegin()->first, false);
        }         

        number_rear += process_simple_path(path);
      }
    } 

    if (number_rear != 0) { 
      isChanged = true;
    } 
  } while (number_rear > 0); 
   
  return isChanged;
} 

template<class graph_t>
size_t Algorithm<graph_t>::ProcessSimplePath::process_simple_path(path_t& path) {
  size_t number_rear = 0;

  if (path.size() >= 4 || (path.size() == 3 && *path.begin() == *path.rbegin())) {
    /*std::cerr << std::endl << "Processing a path of length " << path.size() - 1 << std::endl;
    std::cerr << "path:\t" << *path.begin();
    for(auto ip = ++path.begin(); ip != path.end(); ++ip) {
      std::cerr << " -- " << *ip;
    }
    std::cerr << std::endl;*/
    auto const count_lambda = [&] (vertex_t const & v) -> std::pair<size_t, bool> {
      size_t vtc = 0;
      bool tc = false;
      std::set<mcolor_t> const & colors = this->graph->get_all_multicolor_edge_with_info(*(++path.begin()), v, false); 
      for (mcolor_t const & color : colors) {
        if (this->graph->is_vec_T_consistent_color(color)) { 
          ++vtc;
        } else { 
          tc = true;
        }
      }     
      return std::make_pair(vtc, tc);
    };

    mcolor_t process_color; 
    auto first_edge = count_lambda(*path.begin());
    auto second_edge = count_lambda(*(++++path.begin())); 
    std::set<mcolor_t> process_colors;

    if ((!first_edge.second && second_edge.second) 
      || (!first_edge.second && first_edge.first == std::min(first_edge.first, second_edge.first))) { 
      process_color = this->graph->get_all_multicolor_edge(*(++path.begin()), *path.begin());
      process_colors = this->graph->get_all_multicolor_edge_with_info(*(++path.begin()), *path.begin(), false);
    } else if ((first_edge.second && !second_edge.second)
      || (!second_edge.second && second_edge.first== std::min(first_edge.first, second_edge.first))) { 
      process_color = this->graph->get_all_multicolor_edge(*(++path.begin()), *(++++path.begin()));
      process_colors = this->graph->get_all_multicolor_edge_with_info(*(++path.begin()), *(++++path.begin()), false);
    } 
    
    //std::cerr << "Process color " << genome_match::mcolor_to_name(process_color) << std::endl;

    if ((path.size() % 2 != 0) && (*path.begin() != *path.rbegin())) {
      if (process_color != this->graph->get_all_multicolor_edge(*(++path.begin()), *path.begin())) { 
        //std::cerr << "Erase first" << std::endl;
        path.erase(path.begin());
      } else {
        //std::cerr << "Erase second" << std::endl;
        path.erase(--path.end());
      }
    }

    if (*path.begin() == *path.rbegin()) {
      if (path.size() % 2 == 0) {
        if (process_color == this->graph->get_all_multicolor_edge(*(++path.begin()), *path.begin())) {
          //std::cerr << "... semi-cycle, fusion applied" << std::endl;
          vertex_t const & self_v = *(path.begin());
          vertex_t const & x0 = *(++path.begin());
          vertex_t const & y0 = *(++path.rbegin());

          std::set<mcolor_t> const & colors = this->graph->get_all_multicolor_edge_with_info(x0, self_v, false);
          for (mcolor_t const & color : colors) {
            this->graph->apply(twobreak_t(self_v, x0, self_v, y0, color));
            ++number_rear;
          }

          path.erase(--path.end());
          *path.begin() = y0;
        } else {
          //std::cerr << "... semi-cycle, fission applied" << std::endl;
          vertex_t const & self_v = *(path.begin());
          vertex_t const & y0 = *(++path.rbegin());
          vertex_t const & y1 = *(++++path.rbegin());

          std::set<mcolor_t> const & colors = this->graph->get_all_multicolor_edge_with_info(y0, y1, false);
          for (mcolor_t const & color : colors) {
            this->graph->apply(twobreak_t(y0, y1, self_v, self_v, color));
            ++number_rear;
          }

          path.erase(--path.end());
          *path.rbegin() = self_v;
        }

        if (path.size() < 4) { 
          return number_rear;
        }
      } else { 
        ;//std::cerr << "... cycle" << std::endl;
      }
    }
  
    mcolor_t current_color = this->graph->get_all_multicolor_edge(*(++path.begin()), *path.begin());
    while (process_color != current_color) {
      // multicolor of (z1,z2). N.B.: x2 is NOT oo
      //std::cerr << "... multicolors of first and second multiedges: ";
      if (*path.begin() == *path.rbegin()) {
        //std::cerr << "... rotating" << std::endl;
        path.push_back(*path.begin());  
        path.erase(path.begin());
      } else {
        if (*path.begin() == Infty && *path.rbegin() != Infty) {
          //std::cerr << "... flipping" << std::endl;
          for(auto ip = ++path.begin(); ip != path.end();) {
      	    path.push_front(*ip);
      	    path.erase(ip++);
      	  }
        }

        if (*path.rbegin() == Infty) {
          //std::cerr << "... extending beyond oo" << std::endl;
          path.push_back(Infty);
          path.erase(path.begin());
        } else {
          //std::cerr << "... truncating ??" << std::endl;
          path.erase(path.begin());
          path.erase(--path.end());
          if (path.size() < 4) { 
            return number_rear;
          }
        }
      }
      current_color = this->graph->get_all_multicolor_edge(*(++path.begin()), *path.begin());
    }

    // x1 -- x2 -- x3 -- ... -- x2k
    // results in:
    // x2 == x3   x4 == x5 ... x(2k-2) == x(2k-1)   x1 -- x2k

    auto z3 = path.begin();
    auto z0 = z3++;
    auto z1 = z3++;
    auto z2 = z3++;
    
    while (z3 != path.end()) {
      for (mcolor_t const & color : process_colors) {
      	this->graph->apply(twobreak_t(*z0, *z1, *z3, *z2, color));
      	++number_rear;     
      } 
      z1 = z3++;
      if (z3 != path.end()) { 
        z2 = z3++;
      }
    }

    //std::cerr << "... resolved with " << number_rear << " 2-breaks" << std::endl;
  }
  return number_rear;
}

#endif
