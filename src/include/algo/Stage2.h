#ifndef STAGE_2_H_ 
#define STAGE_2_H_ 

/*
canformQ(x,Q) говорит, можно ли из мультицветов мультиребер инцидентных вершине x образовать мультицвет Q.
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

  Mularcs<Mcolor> M = graph.get_adjacent_multiedges(x, split_bad_colors);
    
  for(auto im = M.cbegin(); im != M.cend(); ++im) { 
    Mcolor C(Q, im->second, Mcolor::Intersection); 
    if (C.size() > 0 && C.size() < im->second.size()) { 
      return false;
    } 
  }
  return true;
}

// test "mobility" of central edge
// can it be ever find neighboring edge of the same multicolor
template<class graph_t>
bool Algorithm<graph_t>::is_mobil_edge(const vertex_t& y, const Mularcs<Mcolor>& Cx, const Mularcs<Mcolor>& Cy) const {
  bool mobilQ = false;

  for(auto jc = Cx.cbegin(); (jc != Cx.cend()) && !mobilQ; ++jc) {
    if (jc->first == y && graph.is_vec_T_color(jc->second)) { //cental sub-edge
      const Mcolor& QQ = jc->second; // color of central sub-edge (QQ is sub-multicolor of Q)
      for(auto ix = Cx.cbegin(); (ix != Cx.cend()) && !mobilQ; ++ix) { 
	if (ix->first != y) { 
	  //std::cerr << "MOBIL: " << x << "-" << ix->first << " canForm: " << genome_match::mcolor_to_name(QQ) << std::endl;
	  mobilQ = canformQ(ix->first, QQ);
	} 
      }
	
      if (!mobilQ) { 
	for(auto iy = Cy.cbegin(); (iy != Cy.cend()) && !mobilQ; ++iy) { 
	  //std::cerr << "MOBIL: " << y << "-" << iy->first << " canForm: " << genome_match::mcolor_to_name(QQ) << std::endl;
	  mobilQ = canformQ(iy->first, QQ);
	}
      }
    } 
  }

  return mobilQ; 
} 

template<class graph_t> 
bool Algorithm<graph_t>::stage2() { 
  bool isChanged = false; 
  size_t nr = 0; // number of rearrangements 

  do {
    nr = 0; 
	
    for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {  
      const vertex_t& x = *is;

      if (graph.is_duplication_vertex(x) || graph.is_indel_vertex(x)) { 
	continue; 
      } 

      Mularcs<Mcolor> M;
      Mularcs<Mcolor> Cx;
      {
	Mularcs<Mcolor> oldM = graph.get_adjacent_multiedges(x);
	Mularcs<Mcolor> oldCx = graph.get_adjacent_multiedges(x, split_bad_colors);	
	for(auto im = oldM.cbegin(); im != oldM.cend(); ++im) { 
	  if (im->second != graph.get_complete_color() && (im->first == Infty 
		|| (im->first != Infty && !graph.is_duplication_vertex(im->first) && !graph.is_indel_vertex(im->first)))) { 	
	    auto ins = oldCx.equal_range(im->first);
	    Cx.insert(ins.first, ins.second);
	    M.insert(im->first, im->second);
	  }
	}
	M.erase(Infty);
      }

      bool found = false;
      for(auto im = M.cbegin(); (im != M.cend()) && !found; ++im) {
	const vertex_t& y = im->first; // Q == im->second - color of central edge

	Mularcs<Mcolor> Cy;
	{
	  Mularcs<Mcolor> oldCy = graph.get_adjacent_multiedges(y, split_bad_colors);
	  oldCy.erase(x);
	  std::unordered_set<vertex_t> processed;
	  for(auto iy = oldCy.cbegin(); iy != oldCy.cend(); ++iy) { 
	    if (processed.find(iy->first) == processed.cend()) {
	      if ((iy->first == Infty) || (iy->first != Infty && !graph.is_duplication_vertex(iy->first) && !graph.is_indel_vertex(iy->first))) {
		auto ins = oldCy.equal_range(iy->first);
		Cy.insert(ins.first, ins.second); 
		processed.insert(iy->first);
	      } 
	    }
	  }
	} 

	if (!is_mobil_edge(y, Cx, Cy)) {
	  //std::cerr << "NOT MOBIL" << std::endl;
	  for (auto ix = Cx.cbegin(); ix != Cx.cend(); ++ix) { 
            vertex_t v = Cy.get_vertex(ix->second);
	    if (ix->first != y && graph.is_vec_T_color(ix->second) && !v.empty()) { 
	      //std::cerr << " Sub-multiedge " << genome_match::mcolor_to_name(ix->second) << std::endl;
	      graph.apply_two_break(TwoBreak<Mcolor>(x, ix->first, y, v, ix->second));
	      found = true;
	      ++nr;
	    }
	  } 
	} 
      }
    } 

    if (nr != 0) { 
      isChanged = true;
    } 
  } while (nr > 0); 
   
  return isChanged;
} 
 
#endif
