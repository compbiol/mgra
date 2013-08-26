#ifndef STAGE7_H_
#define STAGE7_H_

template<class graph_t>
bool Algorithm<graph_t>::stage7() {
  bool isChanged = false; 
  size_t number_rear = 0; // number of rearrangements 

  do {
    number_rear = 0; 
    for (const auto& v: *graph) {  
      if (graph->is_simple_vertex(v)) { 
	path_t path({v});
	std::unordered_set<vertex_t> processed({v, Infty}); // we count oo as already processed
	Mularcs<Mcolor> current = graph->get_adjacent_multiedges(v);

	Mcolor vec_color; 
        if (graph->is_vec_T_consistent_color(current.cbegin()->second)) {
          vec_color = current.cbegin()->second;
        } else if (graph->is_vec_T_consistent_color(current.cbegin()->second)) {
          vec_color = (++current.cbegin())->second;
        } else {
          continue;
        } 
	for(auto im = current.cbegin(); im != current.cend(); ++im) {	
	  bool is_next = (im == current.cbegin()); 
	  vertex_t current = find_less_simple_path(path, processed, v, im->first, vec_color, is_next);
	  if (current == v) { 
	    break; // got a cycle from x to x, cannot extend it 
	  }  		    
	}

        number_rear += convert_less_simple_path(path);
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
vertex_t Algorithm<graph_t>::how_many_paths(const Mularcs<Mcolor>& mularcs, const Mcolor& target_color) {
  //size_t num = 0;
  vertex_t v;  
  for (const auto& arc : mularcs) {
    if (!graph->get_adjacent_multiedges(arc.first).get_vertex(target_color).empty()) { 
      v = arc.first; 
    } 
  } 
  /*if (num != 1) {
    return vertex_t();  
  } */
  return v; 
} 

template<class graph_t>
vertex_t Algorithm<graph_t>::find_less_simple_path(path_t& path, std::unordered_set<vertex_t>& processed, const vertex_t& prev, const vertex_t& cur, Mcolor vec_color, bool is_next) { 
  vertex_t previous  = prev;
  vertex_t current = cur;
  bool stop = true; 
  auto testColorLambda = [&] (const Mularcs<Mcolor>& mularcs) -> bool {
    bool flag = true;
    std::for_each(mularcs.cbegin(), mularcs.cend(), [&] (const std::pair<vertex_t, Mcolor>& arc) -> bool {
      if (!graph->is_vec_T_consistent_color(arc.second)) {
        flag = false;
      } 
    }); 
    return flag;  
  }; 

  while (stop) {
    stop = false; 

    if (is_next) { 
      path.push_front(current);
    } else { 
      path.push_back(current);
    } 

    if (processed.find(current) == processed.end() && !graph->is_duplication_vertex(current)) {     
      processed.insert(current);
      Mularcs<Mcolor> new_edges = graph->get_adjacent_multiedges(current);
      Mcolor previous_color = new_edges.get_multicolor(previous); 
      new_edges.erase(previous);
      
      if (new_edges.size() == 2) { // && testColorLambda(new_edges)) { 
         if (new_edges.union_multicolors() == vec_color) {
           const auto& arc_f = *new_edges.cbegin();
           const auto& arc_s = *(++new_edges.cbegin()); 

           Mularcs<Mcolor> mul_f = graph->get_adjacent_multiedges(arc_f.first); 
           Mularcs<Mcolor> mul_s = graph->get_adjacent_multiedges(arc_s.first);   
          
           const vertex_t& x = mul_f.get_vertex(previous_color);
           const vertex_t& y = mul_s.get_vertex(previous_color);
        
           if (!x.empty() && mul_f.size() == 3 && processed.count(x) == 0) { 
	     if (is_next) { 
               path.push_front(arc_f.first);
             } else { 
               path.push_back(arc_f.first);
             } 
             processed.insert(arc_f.first);
	     previous = arc_f.first;
             current = x; 
             stop = true; 
           } else if (x.empty() && !y.empty() && processed.count(y) == 0 && mul_s.size() == 3) { 
             if (is_next) { 
               path.push_front(arc_s.first);
             } else { 
               path.push_back(arc_s.first);
             } 
	     processed.insert(arc_s.first);
	     previous = arc_s.first;
             current = y;
             stop = true; 
           }  	
           /*vertex_t v;
           do { 
             new_edges.erase(v);
             v = how_many_paths(new_edges, previous_color);
           } while (new_edges.size() != 0 && processed.count(v) != 0 
			&& graph->get_adjacent_multiedges(v).get_multicolors() != graph->get_adjacent_multiedges(current).get_multicolors());

           if (!v.empty()) {
             Mularcs<Mcolor> mularcs_s = graph->get_adjacent_multiedges(v);
             const vertex_t& u = mularcs_s.get_vertex(previous_color);
             mularcs_s.erase(u);

             if (testColorLambda(mularcs_s) && mularcs_s.union_multicolors() == vec_color) { 
               if (is_next) { 
                 path.push_front(v);
               } else { 
                 path.push_back(v);
               } 
            
               processed.insert(v);
	       previous = v;
               current = u;
               stop = true;
             } 
           }*/
         } 
      } else if (new_edges.size() == 1 && graph->get_complement_color(previous_color) == new_edges.cbegin()->second) { 
        if (graph->is_T_consistent_color(new_edges.cbegin()->second) && graph->is_T_consistent_color(previous_color)) { 
 	  previous = current;
	  current = new_edges.cbegin()->first; 
	  stop = true; 
	}  
      }
    }
  }  
             
           

  return current;
} 

template<class graph_t>
size_t Algorithm<graph_t>::convert_less_simple_path(path_t& path) {
  size_t num_rear = 0;
  if (path.size() >= 4 || (path.size() == 3 && *path.begin() == *path.rbegin())) {
        //std::cerr << "path size " << path.size() << std::endl;
        std::cerr << std::endl << "Processing a path of length " << path.size() - 1 << std::endl;
        std::cerr << "path:\t" << *path.begin();
        for(auto ip = ++path.begin(); ip != path.end(); ++ip) {
          std::cerr << " -- " << *ip;
        }
        std::cerr << std::endl;
        vertex_t previous = *path.begin(); 
        vertex_t current = *(++path.begin()); 
        vertex_t next = *(++++path.begin()); 
        for(auto ip = (++++++path.begin()); ip != path.end(); ++ip) {
	  Mularcs<Mcolor> mul_f = graph->get_adjacent_multiedges(current);
          if (mul_f.size() == 3) {
            mul_f.erase(previous);
            mul_f.erase(next);  
            Mularcs<Mcolor> mul_s = graph->get_adjacent_multiedges(next);
            mul_s.erase(current); 
            mul_s.erase(*ip);
            if (mul_f.begin()->second != mul_s.begin()->second) {
              std::cerr << mul_f.begin()->first << " " << mul_s.begin()->first << " " 
              << genome_match::mcolor_to_name(mul_f.begin()->second) << " " << genome_match::mcolor_to_name(mul_s.begin()->second) << std::endl;
            } 
            assert(mul_f.begin()->second == mul_s.begin()->second);
            twobreak_t br2(current, mul_f.begin()->first, next, mul_s.begin()->first, mul_f.begin()->second); 
            graph->apply_two_break(br2);
            ++num_rear;
          }    
          previous = current; 
          current = next; 
          next = *ip;
        }
        std::cerr << num_rear << std::endl;
  } 
  return num_rear; 
} 

#endif
