#ifndef RECOVEREDGENOMES_H_
#define RECOVEREDGENOMES_H_

template<class graph_t>
struct RecoveredGenomes { 
  RecoveredGenomes(const graph_t& gr) 
    : graph(gr)
  { 
  }

  void main_algorithm(std::vector<partgraph_t>& RG);
  std::pair<path_t, bool> getchr(const partgraph_t& PG, const vertex_t& x, std::set<vertex_t>& getchrset);
  std::pair<size_t, size_t> numchr(const partgraph_t& PG);
  void splitchr(const partgraph_t& PG, std::set<std::pair<path_t, bool> >& AllChr, std::list<std::set<vertex_t> >& CircChr);
  void calc_count_2break_type(); 
public: 
  std::list<std::set<vertex_t> > pg_empty;
 
private: 
  const graph_t& graph; 
};

/*void calc_count_2break_type() {
    // number of reversals, interchromosomal translocations, and fissions/fusions
    std::array<size_t, 3> RTF;
    RTF.fill(0);

    for(auto it = graph.crbegin_2break_history(); it != graph.crend_2break_history(); ++it) {
	const Mcolor& Q = it->get_mcolor();
	size_t i = 0;
	for(auto im = graph.cbegin_T_consistent_color(); im != graph.cend_T_consistent_color(); ++im) {
	   if (!Q.includes(*im)) { 
	        ++i;
		continue;
	    }
            //size_t nchr_old = 0;
	    //if (Q == *im) {
		//nchr_old = numchr(graph, RG[i]).first;
	    //}
 it->inverse().apply_single(RG[i]);

	    if (Q == *im) {
		//std::cerr << " " << genome_match::mcolor_to_name(*im);
		
		std::unordered_set<vertex_t> vert;
		if (it->get_arc(0).first != Infty) vert.insert(it->get_arc(0).first);
		if (it->get_arc(0).second != Infty) vert.insert(it->get_arc(0).second);
		if (it->get_arc(1).first != Infty) vert.insert(it->get_arc(1).first);
		if (it->get_arc(1).second != Infty) vert.insert(it->get_arc(1).second);
    
		std::set<vertex_t> getchrset;
		//getchr(RG[i], *vert.begin(), getchrset);
    
		//bool samechr = true;
		vert.erase(vert.begin());
		for(const auto &iv : vert) {
		    if (getchrset.find(iv) == getchrset.end()) {
			samechr = false;
			break;
		    }
		}
		//size_t nchr_new = numchr(graph, RG[i]).first;
		if (nchr_new != nchr_old) {
		    ++RTF[2];
		} else {
		    if (samechr) {
			++RTF[0];
		    } else { 
			++RTF[1];
		    } 
		}
	    }  
std::vector<size_t> tot(3);
    std::cerr << "% Number of reversals / translocations / fissions+fusions: " << std::endl;
    std::cerr << "Total\t&\t" << tot[0] << " & " << tot[1] << " & " << tot[2] << " &\t" << tot[0]+tot[1]+tot[2] << " \\\\" << std::endl;
}*/
 
template<class graph_t>
void RecoveredGenomes<graph_t>::main_algorithm(std::vector<partgraph_t>& RG) {
    for(auto it = graph.crbegin_2break_history(); it != graph.crend_2break_history(); ++it) {
	size_t i = 0;
	for(auto im = graph.cbegin_T_consistent_color(); im != graph.cend_T_consistent_color(); ++im, ++i) {
	    if (it->get_mcolor().includes(*im)) { 
	        it->inverse().apply_single(RG[i]);
	    } 
	}
    }
}

template<class graph_t>
std::pair<size_t, size_t> RecoveredGenomes<graph_t>::numchr(const partgraph_t& PG) {
    std::set<std::pair<path_t, bool> > AllChr;
    std::list<std::set<vertex_t> > CircChr;
    splitchr(PG, AllChr, CircChr);
    return std::make_pair(AllChr.size(), CircChr.size());
}

template<class graph_t>
std::pair<path_t, bool> RecoveredGenomes<graph_t>::getchr(const partgraph_t& PG, const vertex_t& x, std::set<vertex_t>& getchrset) {
    path_t path;
    bool circular = false;
    getchrset.insert(x);

    auto changer_lambda = [&] (const vertex_t& t, bool flag) -> void {
  	if (flag) { 
	((*t.rbegin() == 't')?path.push_back("-" + t.substr(0, t.size() - 1)):path.push_back("+" + t.substr(0, t.size() - 1)));
	} else { 
	((*t.rbegin() == 't')?path.push_front("+" + t.substr(0, t.size() - 1)):path.push_front("-" + t.substr(0, t.size() - 1)));	
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

    return std::make_pair(path, circular);
}

template<class graph_t>
void RecoveredGenomes<graph_t>::splitchr(const partgraph_t& PG, std::set<std::pair<path_t, bool> >& AllChr, std::list<std::set<vertex_t> >& CircChr) {

    if (&CircChr != &pg_empty) { 
	CircChr.clear();
    } 
    AllChr.clear();
    std::unordered_set<vertex_t> processed;

    for(const auto &x : graph) { 
	if (processed.find(x) == processed.end()) { 
		std::set<vertex_t> getchrset;
	        std::pair<path_t, bool> pathb = getchr(PG, x, getchrset);
	
		AllChr.insert(pathb);

	        std::copy(getchrset.begin(), getchrset.end(), std::inserter(processed, processed.end()));
	
		if (pathb.second && (&CircChr != &pg_empty)) {
		    CircChr.push_back(getchrset);
		}
	} 
    }
}



#endif 

