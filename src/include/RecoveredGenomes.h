#ifndef RECOVEREDGENOMES_H_
#define RECOVEREDGENOMES_H_

template<class graph_t>

struct RecoveredGenomes { 

  RecoveredGenomes(const graph_t& gr, const Mcolor& target) 
  : graph(gr)
  {
    if (!target.empty()) { 
      recovered_graphs.resize(1);
      for(const auto &x : graph) {
	    if (recovered_graphs[0].defined(x)) { 
		continue;
	    }

	    vertex_t y = Infty;
	    bool good = true;
	    size_t def = 0;
	    for (const auto &col : target) { 
		size_t numb_color = col.first; 
		if (graph.is_exist_edge(numb_color, x)) {
		    ++def;
		    if (y == Infty) { 
			y = graph.get_adjecent_vertex(numb_color, x);
		    } 

		    if (y != graph.get_adjecent_vertex(numb_color, x)) { 
			good = false;
		    }
		}
	    }

	    if (good && def == target.size() && y != Infty) {
                recovered_graphs[0].insert(x, y);
	    }
	}
    } else {
      recovered_graphs.resize(graph.count_vec_T_consitent_color(), *(graph.cbegin_local_graphs())); 
      for(auto it = graph.crbegin_2break_history(); it != graph.crend_2break_history(); ++it) {
	size_t i = 0;
	for(auto im = graph.cbegin_T_consistent_color(); im != graph.cend_T_consistent_color(); ++im, ++i) {
	  if (it->get_mcolor().includes(*im)) { 
	    it->inverse().apply_single(recovered_graphs[i]);
	  } 
	}
      } 
    }
  }

  void main_algorithm(std::vector<partgraph_t>& RG);

  Chromosome getchr(const partgraph_t& PG, const vertex_t& x, std::set<vertex_t>& getchrset);
  void splitchr(const partgraph_t& PG, Genome& AllChr, std::list<std::set<vertex_t> >& CircChr);

  std::pair<size_t, size_t> numchr(const partgraph_t& PG);
  void calc_count_2break_type(); 

  std::vector<Genome> get_genomes();

private: 
  Genome get_genome(size_t index);
  Chromosome get_chromosome(size_t index, const vertex_t& x, std::unordered_set<vertex_t>& processed);

public: 
  std::list<std::set<vertex_t> > pg_empty;
 
private: 
  const graph_t& graph;
  std::vector<partgraph_t> recovered_graphs;
  std::vector<partgraph_t> recovered_transformation;
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
std::vector<Genome> RecoveredGenomes<graph_t>::get_genomes() { 
  std::vector<Genome> genomes;
  for (size_t i = 0; i < recovered_graphs.size(); ++i) {
    genomes.push_back(get_genome(i)); 
  } 
  return genomes;
} 

template<class graph_t>
Genome RecoveredGenomes<graph_t>::get_genome(size_t index) { 
  Genome genome; 
  std::unordered_set<vertex_t> processed;
  std::string name_chr("chr");
  size_t count = 1; 

  for(const auto &x : graph) { 
    if (processed.find(x) == processed.end()) { 
      std::set<vertex_t> getchrset;
      Chromosome chromosome = get_chromosome(index, x, processed);
      genome.insert(name_chr + toString(count++), chromosome);
    }
  } 

  return genome;
} 

template<class graph_t>
Chromosome RecoveredGenomes<graph_t>::get_chromosome(size_t index, const vertex_t& x, std::unordered_set<vertex_t>& processed) { 
  std::list<std::pair<vertex_t, int> > path;
  bool circular = false;
  processed.insert(x);

  auto changer_lambda = [&] (const vertex_t& t, bool flag) -> void {
    if (flag) { 
      (*t.rbegin() == 't')?path.push_back(std::make_pair(t.substr(0, t.size() - 1), -1)):path.push_back(std::make_pair(t.substr(0, t.size() - 1), 1));
    } else { 
      (*t.rbegin() == 't')?path.push_front(std::make_pair(t.substr(0, t.size() - 1), 1)):path.push_front(std::make_pair(t.substr(0, t.size() - 1), -1));	
    }  
  };

  for(vertex_t y = graph.get_obverse_vertex(x); ;) {
    if (processed.find(y) != processed.end()) {
      circular = true;
      break; // circ
    }
  
    processed.insert(y);

    changer_lambda(y, true);

    if (!recovered_graphs[index].defined(y)) { 
      break; // linear
    } 

    y = recovered_graphs[index][y]; //FIXME

    if (processed.find(y) != processed.end()) {
       circular = true;
       break; // circ
    }
    processed.insert(y);

    if (y == Infty) { 
      break;
    } 
    y = graph.get_obverse_vertex(y);
  }

  if (!circular && recovered_graphs[index].defined(x)) {	
    vertex_t y = x;
    while (recovered_graphs[index].defined(y) && (y != Infty)) {
      y = recovered_graphs[index][y]; //FIXME
      processed.insert(y);

      if (y != Infty) {
        y = graph.get_obverse_vertex(y);
        processed.insert(y);
        changer_lambda(y, false);
      }
    }
  }

  return Chromosome(path, circular);
} 

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
Chromosome RecoveredGenomes<graph_t>::getchr(const partgraph_t& PG, const vertex_t& x, std::set<vertex_t>& getchrset) {
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
void RecoveredGenomes<graph_t>::splitchr(const partgraph_t& PG, Genome& AllChr, std::list<std::set<vertex_t> >& CircChr) {

    if (&CircChr != &pg_empty) { 
	CircChr.clear();
    } 

    std::unordered_set<vertex_t> processed;
    std::string name_chr("chr");
    size_t count = 1; 

    for(const auto &x : graph) { 
	if (processed.find(x) == processed.end()) { 
		std::set<vertex_t> getchrset;
	        Chromosome chromosome = getchr(PG, x, getchrset);
	
		AllChr.insert(name_chr + toString(count++), chromosome);

	        std::copy(getchrset.begin(), getchrset.end(), std::inserter(processed, processed.end()));
	
		if (chromosome.is_circular() && (&CircChr != &pg_empty)) {
		    CircChr.push_back(getchrset);
		}
	} 
    }
}

template<class graph_t>
std::pair<size_t, size_t> RecoveredGenomes<graph_t>::numchr(const partgraph_t& PG) {
    Genome AllChr;
    std::list<std::set<vertex_t> > CircChr;
    splitchr(PG, AllChr, CircChr);
    return std::make_pair(AllChr.count_chromosome(), CircChr.size());
}

#endif 

