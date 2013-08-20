#ifndef STAGE_3_H_
#define STAGE_3_H_

// cut the graph into connected components
template<class graph_t>
bool Algorithm<graph_t>::stage5_1() { 
  bool isChanged = false;
  size_t number_rear = 0; // number of rearrangements 

  do {
    number_rear = 0;
    //std::cerr << "Stage 5_1: splitting into connected components" << std::endl;
    // go over all T-consistent multicolors
    for(auto ic = graph->cbegin_T_consistent_color(); ic != graph->cend_T_consistent_color(); ++ic) {
      const Mcolor& Q = *ic;

      bool repeat = true;
      while(repeat) {
	repeat = false;

	equivalence<vertex_t> CC; // connected components
	std::unordered_map<vertex_t, vertex_t> QQ; // multiedges of colors !Q (!*ic)
	for(const auto &x : *graph) {
	  Mularcs<Mcolor> mularcs = graph->get_adjacent_multiedges(x); 
	  if (mularcs.size() != 1) { 
	    for(const auto &arc : mularcs) {    
	      if (arc.second == *ic) { // edges of color Q (*ic)
		QQ.insert(std::make_pair(x, arc.first));
	      } else if (arc.first != Infty) { // reg. edges of color !Q (!*ic)
		CC.addrel(x, arc.first);
	      }
	    }
	  }
	} 
	CC.update();
		    
	// N.B. for regular edges (x, y) of multicolor Q, QQ[x] = y and QQ[y] = x
	// for irregular edge (x,oo), QQ[x] = oo		    
	//std::cerr << genome_match::mcolor_to_name(*ic) << " ~ " << CC.classes() << std::endl;
    
	std::map<vertex_t, std::set<arc_t> > EC; // reg. edges between diff. connected components of color Q
	std::map<vertex_t, std::set<arc_t> > EI; // irreg. edges of color Q
	for(const auto &edge : QQ) {
	  if (edge.second == Infty) {
	    EI[CC[edge.first]].insert(edge); 
	  } else if (CC.isequiv(edge.first, edge.second)) { 
	    continue;
	  } else {			
	    EC[CC[edge.first]].insert(edge); 
	  }
	}
		
	std::unordered_set<vertex_t> processed;
	// reg. edges between diff. connected components of color Q
	for (const auto &reg_edge : EC) {
	  if (processed.count(reg_edge.first) != 0) { 
	    continue;
	  }

	  // connected component with a single external edge
	  if (reg_edge.second.size() == 1) {
	    const arc_t& p = *(reg_edge.second.begin());
	    arc_t q;

	    if (graph->is_indel_vertex(p.first) || graph->is_indel_vertex(p.second)) {
	      continue;
	    }

	    // look for irregular edges
	    bool found = false;
    
	    for(auto ii = EI[reg_edge.first].cbegin(); (ii != EI[reg_edge.first].cend()) && !found; ++ii) {
	      const arc_t& ireg_edge = *ii; // edge
	      // let check what would happen with (be-)edge e=(p.first, irreg.first)
	      Mcolor color = Q;
	      for(size_t i = 0; i < graph->count_local_graphs(); ++i) {
		if (graph->is_exist_edge(i, p.first) && graph->get_adjecent_vertex(i, p.first) == ireg_edge.first) {
		  color.insert(i);
		}
	      } 
	      // if e is enriched to T-consistent color, great!
	      if (color.size() > Q.size() && graph->is_T_consistent_color(color)) {
		//std::cerr << "perfect edge is found" << std::endl;
		q = ireg_edge;
		found = true;
		EI[reg_edge.first].erase(ii);
	      }
	    }

	    // we did not find good edge, but there are some candidates - take any
	    if (!found && EI[reg_edge.first].size() == 1) {
	      //std::cerr << "somewhat good but unique edge is found" << std::endl;
	      q = *(EI[reg_edge.first].cbegin());
	      EI.erase(reg_edge.first);
	      found = true;
	    }

	    if (!found && EI[reg_edge.first].size() == 0) {
	      //std::cerr << "no irregular edges, do fission" << std::endl;
	      q = std::make_pair(Infty, Infty);
	      found = true;
	    }
    
	    if (found) {	      
	      processed.insert(CC[p.first]);
	      if (q.first != Infty) { 
		processed.insert(CC[q.first]);
	      } 
	      processed.insert(CC[p.second]);
	      if (q.second != Infty) { 
		processed.insert(CC[q.second]);
	      } 
	      //std::cerr << "Stage 5_1: " << p.first << " - " << p.second << "\tX\t" << q.first << " - " << q.second << std::endl;
	      graph->apply_two_break(TwoBreak<Mcolor>(p, q, Q));
	      ++number_rear;			   
	      repeat = true;
            }
	  } else if (reg_edge.second.size() == 2 && EI[reg_edge.first].size() == 0) {
	    const arc_t& p = *(reg_edge.second.begin());
	    const arc_t& q = *(reg_edge.second.rbegin());

	    // N.B. we have CC[p.first] == CC[q.first] == ie->first
	    if ((processed.count(CC[p.second]) == 0) && (processed.count(CC[q.second]) == 0) && CC[p.second] == CC[q.second]) { 
	      processed.insert({CC[p.first], CC[q.first], CC[p.second], CC[q.second]});
	      graph->apply_two_break(TwoBreak<Mcolor>(p, q, Q));
	      ++number_rear;
	      repeat = true;
	      //std::cerr << "Stage 222.2: " << p.first << " - " << p.second << "\tX\t" << q.first << " - " << q.second << endl;
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

//search for 4-cycles
template<class graph_t>
bool Algorithm<graph_t>::stage5_2() { 
  bool isChanged = false;
  size_t number_rear = 0; // number of rearrangements

  do {
    number_rear = 0;

    for(const auto &x : *graph) {
      if (graph->is_duplication_vertex(x) || graph->is_indel_vertex(x)) {
	continue;
      }

      Mularcs<Mcolor> mularcs_x = graph->get_adjacent_multiedges(x);
      bool next = false;

      for(auto im = mularcs_x.cbegin(); (im != mularcs_x.cend()) && (!next); ++im) {
	const vertex_t& y = im->first;
	const Mcolor& Q = im->second;

	if (graph->is_vec_T_consistent_color(Q) && y != Infty && !graph->is_duplication_vertex(y) && !graph->is_indel_vertex(y)) { 
	  Mularcs<Mcolor> mularcs_y = graph->get_adjacent_multiedges(y);
	  mularcs_y.erase(x);

	  for(auto jm = mularcs_y.cbegin(); (jm != mularcs_y.cend()) && (!next); ++jm) {
	    const vertex_t& z = jm->first;
	    if (z != Infty) { 
	      const vertex_t& v = graph->get_adjacent_multiedges(z).get_vertex(Q); 
	      if (!v.empty() && mularcs_x.defined(v)) { 
		//std::cerr << "Stage 5_2: " << x << " - " << y << "\tX\t" << v << " - " << z << std::endl;
		graph->apply_two_break(TwoBreak<Mcolor>(x, y, v, z, Q));
		++number_rear;
		next = true;
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
