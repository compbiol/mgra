#ifndef STAGE_2_H_ 
#define STAGE_2_H_ 

/*
can incident multiedges of x form multicolor Q (current don't check T-consistent formation)
if return false, then Q cannot be formed
if true - who knows...
*/
template<class graph_t>
bool Algorithm<graph_t>::canformQ(const vertex_t& x, const Mcolor& Q) const {
  if (x == Infty) {
    return canformQoo;
  }

  // color Q can be formed if some it adjacent multicolors form a partition of Q
  // OR 
  // if every intersection Q \cap QQ = \emptyset or QQ.

  Mularcs<Mcolor> mularcs = graph->get_adjacent_multiedges(x, split_bad_colors);

  for(const auto &arc : mularcs) { 
    if ((viewed_edges2.count(std::make_pair(x, arc.first)) != 0) || (viewed_edges2.count(std::make_pair(arc.first, x)) != 0)) { 
      return true;
    }
    Mcolor color(Q, arc.second, Mcolor::Intersection); 
    if (color.size() > 0 && color.size() < arc.second.size()) { 
      return false;
    } 
  }
  return true;
}

// test "mobility" of central edge
// can it be ever find neighboring edge of the same multicolor
template<class graph_t>
bool Algorithm<graph_t>::is_mobil_edge(const vertex_t& x, const vertex_t& y, const Mularcs<Mcolor>& mularcs_x, const Mularcs<Mcolor>& mularcs_y) const {
  bool mobilQ = false;

  auto arcs = mularcs_x.equal_range(y); 

  for (auto jc = arcs.first; (jc != arcs.second) && !mobilQ; ++jc) { 
    if (graph->is_vec_T_consistent_color(jc->second)) { //cental sub-edge
      const Mcolor& QQ = jc->second; // color of central sub-edge (QQ is sub-multicolor of Q)

      for(auto ix = mularcs_x.cbegin(); (ix != mularcs_x.cend()) && !mobilQ; ++ix) { 
	if (ix->first != y) { 
	  mobilQ = canformQ(ix->first, QQ);
	  //std::cerr << "MOBIL: " << ix->first << " canForm: " << genome_match::mcolor_to_name(QQ) << " " << QQ.begin()->first<< std::endl;
	} 
      }
	// FIXME : NEAD LAMBDA
      if (!mobilQ) { 
	for(auto iy = mularcs_y.cbegin(); (iy != mularcs_y.cend()) && !mobilQ; ++iy) { 
	  mobilQ = canformQ(iy->first, QQ);
	  //std::cerr << "MOBIL: " << y << "-" << iy->first << " canForm: " << genome_match::mcolor_to_name(QQ) << std::endl;
	}
      }
    } 
  }
  return mobilQ; 
} 

template<class graph_t> 
bool Algorithm<graph_t>::stage2() { 
  bool isChanged = false; 
  size_t number_rear = 0; // number of rearrangements 

  do {
    number_rear = 0; 
	
    for(const auto &x : *graph) {  
      if (graph->is_duplication_vertex(x) || graph->is_indel_vertex(x)) { 
	continue; 
      } 

      Mularcs<Mcolor> mularcs = graph->get_adjacent_multiedges(x);
      Mularcs<Mcolor> mularcs_x = graph->get_adjacent_multiedges(x, split_bad_colors);	
      
      bool found = false;
      for(auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
	const vertex_t& y = im->first; // Q == im->second - color of central edge

if ((x == "134t" && y == "134h") || (x == "134h" && y == "134t")) { 
	std::cerr << x << " " << y << " " << (viewed_edges1.count(std::make_tuple(x, y, im->second)) != 0) << " " << genome_match::mcolor_to_name(im->second) << std::endl;
} 

        if (y != Infty && !graph->is_duplication_vertex(y) && !graph->is_indel_vertex(y)) { 
	  Mularcs<Mcolor> mularcs_y = graph->get_adjacent_multiedges(y, split_bad_colors);
	  mularcs_y.erase(x);
	  if (viewed_edges1.count(std::make_tuple(x, y, im->second)) != 0 || !is_mobil_edge(x, y, mularcs_x, mularcs_y)) {
	    //std::cerr << "NOT MOBIL" << std::endl;
	    for (const auto &arc : mularcs_x) { 
              vertex_t v = mularcs_y.get_vertex(arc.second);
	      //std::cerr << " Sub-multiedge " << v << " " << genome_match::mcolor_to_name(arc.second) << std::endl;
	      if (arc.first != y && graph->is_vec_T_consistent_color(arc.second) && !v.empty()) { 
	        graph->apply_two_break(TwoBreak<Mcolor>(x, arc.first, y, v, arc.second));
	        found = true;
	        ++number_rear;
	      }
	    } 
	  }
	}  
      }
    } 

    if (number_rear != 0) { 
      isChanged = true;
    } 
  } while (number_rear > 0); 
   
  return isChanged;
} 
 
#endif
