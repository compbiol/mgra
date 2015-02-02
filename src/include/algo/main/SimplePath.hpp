#ifndef SIMPLE_PATH_STAGE_HPP
#define SIMPLE_PATH_STAGE_HPP

namespace algo { 

template<class graph_pack_t>
struct ProcessSimplePath : public algo::AbsStage<graph_pack_t> { 
  using path_t = std::list<vertex_t>;
  using mcolor_t = typename graph_pack_t::mcolor_type;
  using mularcs_t = typename graph_pack_t::mularcs_t; 
  using twobreak_t = typename graph_pack_t::twobreak_t;
  
  explicit ProcessSimplePath(size_t max_round)
  : AbsStage<graph_pack_t>("Start process good and simple paths", "simple_path", max_round) 
  {
  }

  bool run(graph_pack_t& graph_pack) override;
  
private:
  size_t process_simple_path(graph_pack_t& graph_pack, path_t& path); 

private:
  DECL_LOGGER("SimplePathStage");
};

template<class graph_pack_t>
bool ProcessSimplePath<graph_pack_t>::run(graph_pack_t& graph_pack) { 
  bool isChanged = false; 
  size_t number_rear = 0; // number of rearrangements 

  do {
    number_rear = 0; 
    for (vertex_t const & v: graph_pack.graph) {  
      if (graph_pack.is_simple_vertex(v)) { 
        path_t path; path.push_back(v);
        std::unordered_set<vertex_t> processed;
        processed.insert(v); processed.insert(Infty); // we count oo as already processed

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

            if (processed.find(current) == processed.end() && !graph_pack.is_duplication_vertex(current)) {     
              processed.insert(current);
              mcolor_t previous_color = graph_pack.get_all_multicolor_edge(current, previous); 
              mularcs_t&& new_edges = graph_pack.get_all_adjacent_multiedges(current);
              new_edges.erase(previous);
              mcolor_t next_color = new_edges.union_multicolors();
      
              if ((graph_pack.graph.degree_vertex(current) == 2) && graph_pack.multicolors.get_complement_color(previous_color) == next_color
                  && !graph_pack.is_postponed_deletion(current, new_edges.cbegin()->first)) {

                auto const check_lambda = [&] (vertex_t const & v) -> bool {
                  bool flag = true; 
                  std::set<mcolor_t> const & colors = graph_pack.get_all_multicolor_edge_with_info(current, v, false); 
                  for (auto color = colors.cbegin(); color != colors.cend() && flag; ++color) {
                    flag = graph_pack.multicolors.is_vec_T_consistent_color(*color);
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

        mularcs_t const & current = graph_pack.get_all_adjacent_multiedges(v);
        vertex_t const & last = find_simple_path_lambda(v, current.cbegin()->first, true);
        if (last != v) { 
          find_simple_path_lambda(v, current.crbegin()->first, false);
        }         

        number_rear += process_simple_path(graph_pack, path);
      }
    } 

    if (number_rear != 0) { 
      isChanged = true;
    } 
  } while (number_rear > 0); 
   
  return isChanged;
} 

template<class graph_pack_t>
size_t ProcessSimplePath<graph_pack_t>::process_simple_path(graph_pack_t& graph_pack, 
      path_t& path) {
  size_t number_rear = 0;

  if (path.size() >= 4 || (path.size() == 3 && *path.begin() == *path.rbegin())) {
    {
      std::stringstream ss; 
      ss << "\nProcessing a path of length " << path.size() - 1 << "\npath:\t" << *path.begin();
      for(auto ip = ++path.cbegin(); ip != path.cend(); ++ip) ss << " -- " << *ip;
      std::string temp = ss.str();
      TRACE(temp) 
    } 

    auto const count_lambda = [&] (vertex_t const & v) -> std::pair<size_t, bool> {
      size_t vtc = 0;
      bool tc = false;
      std::set<mcolor_t> const & colors = graph_pack.get_all_multicolor_edge_with_info(*(++path.begin()), v, false); 
      for (mcolor_t const & color : colors) {
        if (graph_pack.multicolors.is_vec_T_consistent_color(color)) { 
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
      process_color = graph_pack.get_all_multicolor_edge(*(++path.begin()), *path.begin());
      process_colors = graph_pack.get_all_multicolor_edge_with_info(*(++path.begin()), *path.begin(), false);
    } else if ((first_edge.second && !second_edge.second)
      || (!second_edge.second && second_edge.first== std::min(first_edge.first, second_edge.first))) { 
      process_color = graph_pack.get_all_multicolor_edge(*(++path.begin()), *(++++path.begin()));
      process_colors = graph_pack.get_all_multicolor_edge_with_info(*(++path.begin()), *(++++path.begin()), false);
    } 
    
    //std::cerr << "Process color " << genome_match::mcolor_to_name(process_color) << std::endl;

    if ((path.size() % 2 != 0) && (*path.begin() != *path.rbegin())) {
      if (process_color != graph_pack.get_all_multicolor_edge(*(++path.begin()), *path.begin())) { 
        TRACE("Erase first")
        path.erase(path.begin());
      } else {
        TRACE("Erase second")
        path.erase(--path.end());
      }
    }

    if (*path.begin() == *path.rbegin()) {
      if (path.size() % 2 == 0) {
        if (process_color == graph_pack.get_all_multicolor_edge(*(++path.begin()), *path.begin())) {
          TRACE("... semi-cycle, fusion applied");
          vertex_t const & self_v = *(path.begin());
          vertex_t const & x0 = *(++path.begin());
          vertex_t const & y0 = *(++path.rbegin());

          std::set<mcolor_t> const & colors = graph_pack.get_all_multicolor_edge_with_info(x0, self_v, false);
          for (mcolor_t const & color : colors) {
            graph_pack.apply(twobreak_t(self_v, x0, self_v, y0, color));
            ++number_rear;
          }

          path.erase(--path.end());
          *path.begin() = y0;
        } else {
          TRACE("... semi-cycle, fission applied");
          vertex_t const & self_v = *(path.begin());
          vertex_t const & y0 = *(++path.rbegin());
          vertex_t const & y1 = *(++++path.rbegin());

          std::set<mcolor_t> const & colors = graph_pack.get_all_multicolor_edge_with_info(y0, y1, false);
          for (mcolor_t const & color : colors) {
            graph_pack.apply(twobreak_t(y0, y1, self_v, self_v, color));
            ++number_rear;
          }

          path.erase(--path.end());
          *path.rbegin() = self_v;
        }

        if (path.size() < 4) { 
          return number_rear;
        }
      } else { 
        TRACE("... cycle");
      }
    }
  
    mcolor_t current_color = graph_pack.get_all_multicolor_edge(*(++path.begin()), *path.begin());
    while (process_color != current_color) {
      // multicolor of (z1,z2). N.B.: x2 is NOT oo
      TRACE("... multicolors of first and second multiedges: ");
      if (*path.begin() == *path.rbegin()) {
        TRACE("... rotating");
        path.push_back(*path.begin());  
        path.erase(path.begin());
      } else {
        if (*path.begin() == Infty && *path.rbegin() != Infty) {
          TRACE("... flipping");
          for(auto ip = ++path.begin(); ip != path.end();) {
      	    path.push_front(*ip);
      	    path.erase(ip++);
      	  }
        }

        if (*path.rbegin() == Infty) {
          TRACE("... extending beyond oo");
          path.push_back(Infty);
          path.erase(path.begin());
        } else {
          TRACE("... truncating ??");
          path.erase(path.begin());
          path.erase(--path.end());
          if (path.size() < 4) { 
            return number_rear;
          }
        }
      }
      current_color = graph_pack.get_all_multicolor_edge(*(++path.begin()), *path.begin());
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
      	graph_pack.apply(twobreak_t(*z0, *z1, *z3, *z2, color));
      	++number_rear;     
      } 
      z1 = z3++;
      if (z3 != path.end()) { 
        z2 = z3++;
      }
    }

    TRACE("... resolved with " << number_rear << " 2-breaks");
  }
  return number_rear;
}

} 
#endif
