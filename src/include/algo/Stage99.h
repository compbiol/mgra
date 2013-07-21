#ifndef STAGE_99_H_
#define STAGE_99_H_

// cut hanging free ends
template<class graph_t>
bool Algorithm<graph_t>::cut_free_ends() {
  bool simplified = false;
  size_t nf = 0; // number of fussions/fissions

  do {
	nf = 0;
	for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {  
		const std::string& x = *is;
		Mularcs<Mcolor> M = graph.get_adjacent_multiedges(x);
	 
		const Mcolor& Q1 = M.find(Infty)->second;
                vertex_t y;
		Mcolor Q2;
		if( M.cbegin()->first == Infty) {
                    y = M.crbegin()->first;
		    Q2 = M.crbegin()->second;
		} else {
                    y = M.cbegin()->first;
		    Q2 = M.cbegin()->second;
		}
		if( !member(graph.DiColor,Q1) && member(graph.DiColor,Q2) /* && !member(graph.get_adjacent_multiedges(y),Infty) */ ) {
		//if( member(DiColor,Q2) && !member(graph.get_adjacent_multiedges(y),Infty) ) {
		    outlog << "Unhanging fission:" << endl;
		    //if( TwoBreak<graph_t, Mcolor>(x,y,Infty,Infty,Q2).apply(graph,true) ) nf++;
		    graph.apply_two_break(TwoBreak<Mcolor>(x, y, Infty, Infty, Q2), true); 
		    ++nf;
		}
	    }
	if (nf > 0) simplified = true;
  } while (nf > 0);

  return simplified; 
} 

template<class graph_t>
bool Algorithm<graph_t>::find_reliable_path() { 
  bool simplified = false;
  size_t nr = 0;
  size_t nf = 0; // number of fussions/fissions

  do {
	nr = 0;
	nf = 0; 

	for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {  
	     const std::string& x = *is;
	     Mularcs<Mcolor> M = graph.get_adjacent_multiedges(x);
	    // generalized reliable simple path	
	    if (M.size()>=2 ) {
			for(auto im = M.cbegin();im!=M.cend();++im) {
			  const std::string& y = im->first;
   		          const Mcolor& Q = im->second;
	
	                  if( y==Infty ) continue;
	                    
			  Mularcs<Mcolor> My = graph.get_adjacent_multiedges(y);
			  My.erase(x);
     			  Mularcs<Mcolor> Mx = M;
			  Mx.erase(y);

		    if((Mx.find(Infty) != Mx.cend()) && (My.find(Infty) != My.cend()) && (Mx.find(Infty)->second == My.find(Infty)->second) 
			&& graph.is_T_consistent_color(Mx.find(Infty)->second) ) {
			Mcolor C(Q, Mx.find(Infty)->second, Mcolor::Union);
			if (!graph.is_T_consistent_color(C)) continue;
			for(auto iq = Mx.find(Infty)->second.cbegin(); iq!=Mx.find(Infty)->second.cend(); ++iq) {
			    graph.add_edge(iq->first, x, y);
			}
			outlog << "Stage 22: fusion " << x << " + " << y << std::endl;
			nf++;
			break;
		    }
		    
		    auto Cx = Mx.cend();
		    for(auto jm = Mx.cbegin();jm!=Mx.cend();++jm) {
			if( !graph.is_T_consistent_color(jm->second) ) continue;
			Mcolor C(Q, jm->second, Mcolor::Union);
			if (graph.is_T_consistent_color(C)) {
			    if( Cx!=Mx.cend() ) { Cx=Mx.cend(); break; }
			    Cx = jm;
			}
	    	    }
	    	    if( Cx == Mx.cend() ) continue;
		    
		    auto Cy = My.cend();
		    for(auto jm = My.cbegin(); jm != My.cend(); ++jm) {
			if( !graph.is_T_consistent_color(jm->second) ) continue;
			Mcolor C(Q, jm->second, Mcolor::Union);
			if( graph.is_T_consistent_color(C) ) {
			    if( Cy!=My.cend() ) { Cy=My.cend(); break; }
			    Cy = jm;
			}
	    	    }
	    	    if( Cy == My.cend() ) continue;
	    	    
	    	    if( Cx->second == Cy->second ) {
                        //if( TwoBreak<graph_t, Mcolor>(x,Cx->first,y,Cy->first,Cx->second).apply(graph,true) ) nr++;
			graph.apply_two_break(TwoBreak<Mcolor>(x, Cx->first, y, Cy->first, Cx->second), true);
			nr++;
	    		outlog << "Stage 22: fusion " << x << " + " << y << endl;
	    		break;
	    	    }
		}
	    }
	} 
	if (nf > 0 || nr > 0) simplified = true;
  } while (nf > 0 || nr > 0);

  return simplified; 

} 
#endif
