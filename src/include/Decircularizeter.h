#ifndef DECIRCULARIZETER_H_
#define DECIRCULARIZETER_H_  

template<class graph_t>
struct Decircularizeter {

  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::edge_t edge_t;
  typedef typename graph_t::twobreak_t twobreak_t; 
  typedef typename graph_t::transform_t transform_t;
  typedef typename graph_t::partgraph_t partgraph_t; 
  
  Decircularizeter(graph_t const & gr, partgraph_t const & b_edges) 
  : graph(gr)
  , bad_edges(b_edges)
  { 
  }

  transform_t decircularize(partgraph_t& local_graph, transform_t& transform);
  size_t count_circular_chromosome(partgraph_t const & local_graph) const;

private: 
  bool is_circular_chromosome(partgraph_t const & local_graph, vertex_t const & x, std::unordered_set<vertex_t>& processed) const;
private: 
  graph_t const & graph;
  partgraph_t const & bad_edges;
}; 

template<class graph_t>
bool Decircularizeter<graph_t>::is_circular_chromosome(partgraph_t const & local_graph, vertex_t const & x, std::unordered_set<vertex_t>& processed) const {
  bool circular = false;
  bool have_deletion = false;
  vertex_t current = x; 
  vertex_t previous = x;
  processed.insert(x);

  do { 
    previous = current; 
    current = graph.get_obverse_vertex(previous);
    if (processed.count(current) == 0) {
      processed.insert(current);
      if (local_graph.defined(current)) {
        previous = current; 
        current = local_graph[previous];
        if (bad_edges.defined(previous, current)) {
          have_deletion = true;
        }
        if (processed.count(current) != 0) {
          circular = true;
        } 
        processed.insert(current); 
      }
    } else { 
      circular = true;
    } 
  } while ((current != Infty) && local_graph.defined(current) && !circular);

  if (!circular && local_graph.defined(x)) {	
    for (vertex_t y = local_graph[x]; local_graph.defined(y) && (y != Infty); y = local_graph[y]) {
      processed.insert(y);
      if (y != Infty) {
      	y = graph.get_obverse_vertex(y);
    	  processed.insert(y);
      }
    }
  } 

  if (have_deletion) { 
    return false;
  } 

  return circular;
}

template<class graph_t>
size_t Decircularizeter<graph_t>::count_circular_chromosome(partgraph_t const & local_graph) const {
  std::unordered_set<vertex_t> processed;
  size_t count_chr = 0;
  
  for (auto const & x : graph) { 
    if (processed.count(x) == 0) { 
      std::unordered_set<vertex_t> chr_set;
      bool circular = is_circular_chromosome(local_graph, x, chr_set); 
      std::copy(chr_set.begin(), chr_set.end(), std::inserter(processed, processed.end()));
      if (circular) {
        ++count_chr; 
      }
    } 
  }
  
  return count_chr;
}

/* Given a non-linear genome PG of a multicolor Q and a transformation into a linear genome,
 * find linearizing fissions in the transformation, move them to beginning, and apply to PG
 * i.e., try to reorder transformation into: PG -> PG' -> linear genome, where PG' is linear
 * and obtained from PG with fission.
 * We replace PG with PG' and return the transformation PG -> PG'
 * Transformation may contain only multicolors Q' with Q'\cap Q = 0 or Q.
*/
template<class graph_t>
typename graph_t::transform_t Decircularizeter<graph_t>::decircularize(partgraph_t& PG, transform_t& TG) {
  // decircularizing sub-transform that is removed
  transform_t D;

  size_t CircSize = count_circular_chromosome(PG);
  if (CircSize == 0) {
	  return D;
  } 

 // std::cerr << "Eliminating " << CircSize << " circular chromosomes in " << genome_match::mcolor_to_name(Q) << std::endl;

  partgraph_t T = PG; // reconstructed genome ("bad")
  auto start = TG.begin();

  // looking for de-circularizig 2-breaks
  for(auto it = start; it != TG.end();) {
    it->apply_single(T);

    size_t ccsize = count_circular_chromosome(T);

		if (ccsize >= CircSize) {
    	++it;
    	continue;
		}

		//std::cerr << "Found problematic 2-break: ";// << *it << "\t";

		// move t over to beginning
		for(auto jt = it; jt != TG.begin();) {

	    auto kt = jt--; // jt, kt are successive, *kt == t

      twobreak_t const & s = *jt;
	    twobreak_t const & t = *kt;

			//outlog << "... trying to swap with " << s << endl;

	    edge_t p1, q1, p2, q2;

	    bool usearc = false;

	    mcolor_t C(t.get_mcolor(), s.get_mcolor(), mcolor_t::Intersection);
	    if (!C.empty()) {

			/*
				 p1=(x1,x2) x (y1,y2)=q1
				 p2=(x1,y1) x (x3,y3)=q2
	    
				 into:
	    
				 (x1,x2) x (x3,y3)
				 (y3,x2) x (y1,y2)
			*/
    
  			for(int j = 0; j < 2; ++j) {    
  			  if (t.get_arc(j) == std::make_pair(jt->get_vertex(0), jt->get_vertex(2))) { 
  					usearc = true;
  	    
  					q2 = t.get_arc(1 - j);
  		    
  					p1 = jt->get_arc(0);
  					q1 = jt->get_arc(1);
  			  } else if (t.get_arc(j) == std::make_pair(jt->get_vertex(2), jt->get_vertex(0))) {
  					usearc = true;
  		    
  					q2 = t.get_arc(1 - j);
  		    
  					p1 = jt->get_arc(1);
  					q1 = jt->get_arc(0);
  			  } else if (t.get_arc(j) == std::make_pair(jt->get_vertex(1), jt->get_vertex(3))) {
  					usearc = true;
  		    
  					q2 = t.get_arc(1 - j);
  		    
  					p1 = std::make_pair(jt->get_vertex(1), jt->get_vertex(0));
  					q1 = std::make_pair(jt->get_vertex(3), jt->get_vertex(2));
  			  } else if (t.get_arc(j) == std::make_pair(jt->get_vertex(3), jt->get_vertex(1))) {
  					usearc = true;
  		    
  					q2 = t.get_arc(1 - j);
  		    
  					p1 = std::make_pair(jt->get_vertex(3), jt->get_vertex(2));
  					q1 = std::make_pair(jt->get_vertex(1), jt->get_vertex(0));
  			  }

  			  if (usearc) { 
  			    break; 
  			  }
  			}
      }

		    // TwoBreak t0 = t;
		    if (usearc) {
  				if (t.get_mcolor() != s.get_mcolor())  { 
  					break;
  				}
          //std::cerr << "Insert " << std::endl;
          //std::cerr << q2.second << " " << p1.second << " " << q1.first << " " << q1.second << genome_match::mcolor_to_name(t.get_mcolor()) << std::endl;
          //std::cerr << p1.first << " " << p1.second << " " << q2.first << " " << q2.second << genome_match::mcolor_to_name(t.get_mcolor()) << std::endl;
          *jt = twobreak_t(p1.first, p1.second, q2.first, q2.second, t.get_mcolor());
  				*kt = twobreak_t(q2.second, p1.second, q1.first, q1.second, t.get_mcolor());
		    } else {
          //std::cerr << "Swap " << std::endl;
  				twobreak_t temp = *kt;
  				*kt = *jt;
          *jt = temp;
        }

		    {
  				mcolor_t C(kt->get_mcolor(), it->get_mcolor(), mcolor_t::Intersection);
  	    
          // N.B. at this point if C is not empty, then C == Q
  				if (!C.empty()) {
  			  	kt->inverse().apply_single(T);
  					ccsize = count_circular_chromosome(T);
  				}
		    }

	    /*
	    if( CC.size() > ccsize ) {
		    outlog << "Cannot pop:" << endl;
		    outlog << *jt << " , " << t0 << "  -->  " << t << " , " << *kt << endl;
	    }
	    */
    }

		if (ccsize < CircSize) {
    	//std::cerr << " SUCCEDED" << std::endl;
    	// move t away from the transformation TG and save it to D
      TG.begin()->apply_single(PG);
    	D.push_back(*TG.begin());

    	TG.erase(TG.begin());

    	CircSize = count_circular_chromosome(PG);

    	if (CircSize == 0) { 
      		break;
    	} 

    	start = TG.begin();
		} else {  // did not succeed
    	start++;
	    //std::cerr << " FAILED" << std::endl;
		}

		T = PG;
		for(it = TG.begin(); it != start; ++it) {
      it->apply_single(T);
		}
  }
   
    //if (CircSize > 0) {
	//std::cerr << "Unable to de-circularize ;-(" << std::endl;
    //}

  return D;
}

#endif
