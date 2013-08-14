#ifndef DECIRCULARIZETER_H_
#define DECIRCULARIZETER_H_  

template<class graph_t>
struct Decircularizeter {
  Decircularizeter(const graph_t& gr) 
  : graph(gr)
  { 
  }

  transform_t decircularize(partgraph_t& PG, transform_t& TG);
private: 
  Chromosome getchr(const partgraph_t& PG, const vertex_t& x, std::set<vertex_t>& getchrset);
  void splitchr(const partgraph_t& PG, Genome& AllChr, std::list<std::set<vertex_t> >& CircChr);
  std::pair<size_t, size_t> numchr(const partgraph_t& PG); 
private: 
  const graph_t& graph;
//partgraph_t local_genome;
//transform_t local_transform;
}; 

template<class graph_t>
Chromosome Decircularizeter<graph_t>::getchr(const partgraph_t& PG, const vertex_t& x, std::set<vertex_t>& getchrset) {
    std::list<std::pair<vertex_t, int> > path;
    bool circular = false;
    getchrset.insert(x);

    auto changer_lambda = [&] (const vertex_t& t, bool flag) -> void {
  	if (flag) { 
	(*t.rbegin() == 't')?path.push_back(std::make_pair(t.substr(0, t.size() - 1), -1)):path.push_back(std::make_pair(t.substr(0, t.size() - 1), 1));
	} else { 
	(*t.rbegin() == 't')?path.push_front(std::make_pair(t.substr(0, t.size() - 1), 1)):path.push_front(std::make_pair(t.substr(0, t.size() - 1), -1));	
	}  
    };

    for(vertex_t y = graph.get_obverse_vertex(x); ; ) {
	if (getchrset.find(y) != getchrset.end()) {
	    circular = true;
	    break; // circ
	}
	getchrset.insert(y);

    	changer_lambda(y, true);

	if (!PG.defined(y)) { 
	  break; // linear
	} 

	y = PG[y];

	if (getchrset.find(y) != getchrset.end()) {
	    circular = true;
	    break; // circ
	}
	getchrset.insert(y);

	if (y == Infty) { 
		break;
	} 
	y = graph.get_obverse_vertex(y);
    }

    if (!circular && PG.defined(x)) {	
	vertex_t y = x;
	while (PG.defined(y) && (y != Infty)) {
	    y = PG[y];
	    getchrset.insert(y);

	    if (y != Infty) {
	    	y = graph.get_obverse_vertex(y);
	    	getchrset.insert(y);

    	    	changer_lambda(y, false);
	    }
	}
    }

    return Chromosome(path, circular);
}

template<class graph_t>
void Decircularizeter<graph_t>::splitchr(const partgraph_t& PG, Genome& AllChr, std::list<std::set<vertex_t> >& CircChr) {
    std::unordered_set<vertex_t> processed;
    std::string name_chr("chr");
    size_t count = 1; 

    for(const auto &x : graph) { 
	if (processed.find(x) == processed.end()) { 
		std::set<vertex_t> getchrset;
	        Chromosome chromosome = getchr(PG, x, getchrset);
	
		AllChr.insert(name_chr + toString(count++), chromosome);

	        std::copy(getchrset.begin(), getchrset.end(), std::inserter(processed, processed.end()));
	
		if (chromosome.is_circular()) {
		    CircChr.push_back(getchrset);
		}
	} 
    }
}

template<class graph_t>
std::pair<size_t, size_t> Decircularizeter<graph_t>::numchr(const partgraph_t& PG) {
    Genome AllChr;
    std::list<std::set<vertex_t> > CircChr;
    splitchr(PG, AllChr, CircChr);
    return std::make_pair(AllChr.count_chromosome(), CircChr.size());
}

/* Given a non-linear genome PG of a multicolor Q and a transformation into a linear genome,
 * find linearizing fissions in the transformation, move them to beginning, and apply to PG
 * i.e., try to reorder transformation into: PG -> PG' -> linear genome, where PG' is linear
 * and obtained from PG with fission.
 * We replace PG with PG' and return the transformation PG -> PG'
 * Transformation may contain only multicolors Q' with Q'\cap Q = 0 or Q.
*/
template<class graph_t>
transform_t Decircularizeter<graph_t>::decircularize(partgraph_t& PG, transform_t& TG) {
    // decircularizing sub-transform that is removed
    transform_t D;

    size_t CircSize = numchr(PG).second;
    if (CircSize == 0) {
	return D;
    } 

    //std::cerr << "Eliminating " << CircSize << " circular chromosomes in " << genome_match::mcolor_to_name(Q) << std::endl;

    partgraph_t T = PG; // reconstructed genome ("bad")
    auto start = TG.begin();

    // looking for de-circularizig 2-breaks
    for(auto it = start; it != TG.end();) {
        it->apply_single(T);

        size_t ccsize = numchr(T).second;

	if (ccsize >= CircSize) {
	    ++it;
	    continue;
	}

	//std::cerr << "Found problematic 2-break: ";// << *it << "\t";

	// move t over to beginning
	for(auto jt = it; jt != TG.begin();) {

	    auto kt = jt--; // jt, kt are successive, *kt == t

	    const TwoBreak<Mcolor>& t = *kt;
	    const TwoBreak<Mcolor>& s = *jt;

//            outlog << "... trying to swap with " << s << endl;

	    std::pair<vertex_t, vertex_t> p1, q1, p2, q2;

	    bool usearc = false;

	    Mcolor C(t.get_mcolor(), s.get_mcolor(), Mcolor::Intersection);
	    if (!C.empty()) {


		/*
			 p1=(x1,x2) x (y1,y2)=q1
			 p2=(x1,y1) x (x3,y3)=q2
    
			 into:
    
			 (x1,x2) x (x3,y3)
			 (y3,x2) x (y1,y2)
		*/
    
		for(int j = 0; j < 2; ++j) {    
		    if (t.get_arc(j) == std::make_pair(jt->get_arc(0).first, jt->get_arc(1).first)) { 
			usearc = true;
    
			p2 = t.get_arc(j);
			q2 = t.get_arc(1 - j);
    
			p1 = jt->get_arc(0);
			q1 = jt->get_arc(1);
		    } else if (t.get_arc(j) == std::make_pair(jt->get_arc(1).first, jt->get_arc(0).first)) {
			usearc = true;
    
			p2 = t.get_arc(j);
			q2 = t.get_arc(1 - j);
    
			p1 = jt->get_arc(1);
			q1 = jt->get_arc(0);
		    } else if (t.get_arc(j) == std::make_pair(jt->get_arc(0).second, jt->get_arc(1).second)) {
			usearc = true;
    
			p2 = t.get_arc(j);
			q2 = t.get_arc(1 - j);
    
			p1 = std::make_pair(jt->get_arc(0).second, jt->get_arc(0).first);
			q1 = std::make_pair(jt->get_arc(1).second, jt->get_arc(1).first);
		    } else if (t.get_arc(j) == std::make_pair(jt->get_arc(1).second, jt->get_arc(0).second)) {
			usearc = true;
    
			p2 = t.get_arc(j);
			q2 = t.get_arc(1 - j);
    
			p1 = std::make_pair(jt->get_arc(1).second, jt->get_arc(1).first);
			q1 = std::make_pair(jt->get_arc(0).second, jt->get_arc(0).first);
		    }
		    if (usearc) break;
		}
	    }

	    // TwoBreak t0 = t;

	    if (usearc) {
		if (t.get_mcolor() != s.get_mcolor()) break;
		*kt = TwoBreak<Mcolor>(q2.second, p1.second, q1.first, q1.second, t.get_mcolor());
		*jt = TwoBreak<Mcolor>(p1.first, p1.second, q2.first, q2.second, t.get_mcolor());
	    } else {
		TwoBreak<Mcolor> temp = *kt;
		*kt = *jt;
                *jt = temp;
	    }

	    {
		Mcolor C(kt->get_mcolor(), it->get_mcolor(), Mcolor::Intersection);
    
                // N.B. at this point if C is not empty, then C == Q
		if (!C.empty()) {
		    kt->inverse().apply_single(T);

		    ccsize = numchr(T).second;
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

	    CircSize = numchr(PG).second;

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
