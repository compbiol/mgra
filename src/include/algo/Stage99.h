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
		mularcs_t M = graph.get_adjacent_multiedges(x);
	 
		const Mcolor& Q1 = M[Infty];
                vertex_t y;
		Mcolor Q2;
		if( M.begin()->first == Infty) {
                    y = M.rbegin()->first;
		    Q2 = M.rbegin()->second;
		} else {
                    y = M.begin()->first;
		    Q2 = M.begin()->second;
		}
		if( !member(graph.DiColor,Q1) && member(graph.DiColor,Q2) /* && !member(graph.get_adjacent_multiedges(y),Infty) */ ) {
		//if( member(DiColor,Q2) && !member(graph.get_adjacent_multiedges(y),Infty) ) {
		    outlog << "Unhanging fission:" << endl;
		    if( TwoBreak(x,y,Infty,Infty,Q2).apply(graph,true) ) nf++;
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
	     mularcs_t M = graph.get_adjacent_multiedges(x);
	    // generalized reliable simple path	
	    if (M.size()>=2 ) {
			for(auto im = M.begin();im!=M.end();++im) {
			  const std::string& y = im->first;
   		          const Mcolor& Q = im->second;
	
	                  if( y==Infty ) continue;
	                    
			  mularcs_t My = graph.get_adjacent_multiedges(y);
			  My.erase(x);
     			  mularcs_t Mx = M;
			  Mx.erase(y);

		    if( member(Mx,Infty) && member(My,Infty) && Mx[Infty]==My[Infty] && graph.is_T_consistent_color(Mx[Infty]) ) {
			Mcolor C(Q, Mx[Infty], Mcolor::Union);
			if (!graph.is_T_consistent_color(C)) continue;
			for(auto iq = Mx[Infty].begin(); iq!=Mx[Infty].end(); ++iq) {
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
			    if( Cx!=Mx.end() ) { Cx=Mx.end(); break; }
			    Cx = jm;
			}
	    	    }
	    	    if( Cx == Mx.end() ) continue;
		    
		    auto Cy = My.cend();
		    for(auto jm = My.cbegin(); jm != My.cend(); ++jm) {
			if( !graph.is_T_consistent_color(jm->second) ) continue;
			Mcolor C(Q, jm->second, Mcolor::Union);
			if( graph.is_T_consistent_color(C) ) {
			    if( Cy!=My.end() ) { Cy=My.end(); break; }
			    Cy = jm;
			}
	    	    }
	    	    if( Cy == My.end() ) continue;
	    	    
	    	    if( Cx->second == Cy->second ) {
                        if( TwoBreak(x,Cx->first,y,Cy->first,Cx->second).apply(graph,true) ) nr++;
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
