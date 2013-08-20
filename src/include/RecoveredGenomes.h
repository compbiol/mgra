#ifndef RECOVEREDGENOMES_H_
#define RECOVEREDGENOMES_H_

#include "Decircularizeter.h"

template<class graph_t>
struct RecoveredGenomes {
  typedef std::list<TwoBreak<Mcolor> > transform_t;

  RecoveredGenomes(const graph_t& gr, const Mcolor& target, const std::map<Mcolor, std::set<arc_t> >& b_edges);

  std::vector<Genome> get_genomes();

  std::vector<transform_t> get_history() const { 
    return recovered_transformation;
  }

private: 
  Genome get_genome(size_t index);
  Chromosome get_chromosome(size_t index, const vertex_t& x, std::unordered_set<vertex_t>& chromosome_set);

private: 
  const graph_t& graph;
  const std::map<Mcolor, std::set<arc_t> >& bad_edges;
  std::vector<std::string> name_genomes;
  std::vector<partgraph_t> recovered_graphs;
  std::vector<transform_t> recovered_transformation;
};

template<class graph_t>
RecoveredGenomes<graph_t>::RecoveredGenomes(const graph_t& gr, const Mcolor& target, const std::map<Mcolor, std::set<arc_t> >& b_edges)  
: graph(gr)
, bad_edges(b_edges)
{
 if (!target.empty()) { 
   name_genomes.push_back(genome_match::mcolor_to_name(target));
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
      recovered_transformation.resize(graph.count_vec_T_consitent_color());
      recovered_graphs.resize(graph.count_vec_T_consitent_color(), *(graph.cbegin_local_graphs())); 

      for(auto it = graph.crbegin_2break_history(); it != graph.crend_2break_history(); ++it) {
	size_t i = 0;
	for(auto im = graph.cbegin_T_consistent_color(); im != graph.cend_T_consistent_color(); ++im, ++i) {
	  if (it->get_mcolor().includes(*im)) { 
	    it->inverse().apply_single(recovered_graphs[i]);
	  }
	  if (it->get_mcolor() == *im) {
	    recovered_transformation[i].push_front(*it);
          }
	}
      }

      size_t i = 0; 
      for (auto im = graph.cbegin_T_consistent_color(); im != graph.cend_T_consistent_color(); ++im, ++i) {    
        name_genomes.push_back(genome_match::mcolor_to_name(*im));
        std::set<arc_t> edges; 
        if (bad_edges.count(*im) != 0) {
          edges = bad_edges.find(*im)->second;
        } 

        Decircularizeter<graph_t> dec(graph, edges);        
        //std::cerr << "Initial we have " << dec.count_circular_chromosome(recovered_graphs[i]) << std::endl;
	transform_t T = dec.decircularize(recovered_graphs[i], recovered_transformation[i]);
        //std::cerr << "After we have " << dec.count_circular_chromosome(recovered_graphs[i]) << std::endl;

	// move to adjacent branches
	for (const auto &it : T) {
	  size_t j = 0; 
	  for (auto jt = graph.cbegin_T_consistent_color(); jt != graph.cend_T_consistent_color(); ++jt, ++j) {
	    if ((j != i) && includes(im->cbegin(), im->cend(), jt->cbegin(), jt->cend()) && graph.are_adjacent_branches(*im, *jt)) {
	      recovered_transformation[j].push_back(it);
	    }
	  }
	}
      } 
   }
}

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
  Genome genome(name_genomes[index]); 
  std::unordered_set<vertex_t> processed;
  std::string name_chr("chr");
  size_t count = 1; 

  for(const auto &x : graph) { 
    if (processed.find(x) == processed.end()) { 
      std::unordered_set<vertex_t> chromosome_set;
      Chromosome chromosome = get_chromosome(index, x, chromosome_set);
      if (chromosome.size() != 0) {
        genome.insert(name_chr + toString(count++), chromosome);
      } 
      std::copy(chromosome_set.cbegin(), chromosome_set.cend(), std::inserter(processed, processed.end()));
    }
  } 

  return genome;
} 

template<class graph_t>
Chromosome RecoveredGenomes<graph_t>::get_chromosome(size_t index, const vertex_t& x, std::unordered_set<vertex_t>& chromosome_set) {
    std::list<std::pair<vertex_t, int> > path;
    bool circular = false;
    auto changer_lambda = [&] (const vertex_t& t, bool flag) -> void {
  	if (flag) { 
	(*t.rbegin() == 't')?path.push_back(std::make_pair(t.substr(0, t.size() - 1), -1)):path.push_back(std::make_pair(t.substr(0, t.size() - 1), 1));
	} else { 
	(*t.rbegin() == 't')?path.push_front(std::make_pair(t.substr(0, t.size() - 1), 1)):path.push_front(std::make_pair(t.substr(0, t.size() - 1), -1));	
	}  
    };

  /*do { 
    previous = current; 
    current = graph.get_obverse_vertex(previous);
    if (bad_edges.count(std::make_pair(previous, current)) != 0 || bad_edges.count(std::make_pair(current, previous)) != 0) {
      have_deletion = true;
    }
    if (chromosome_set.count(current) == 0) {
      chromosome_set.insert(current);
      if (recovered_graphs[index].defined(current)) {
        previous = current; 
        current = recovered_graphs[index][previous];
        changer_lambda(current, true);
        if (bad_edges.count(std::make_pair(previous, current)) != 0 || bad_edges.count(std::make_pair(current, previous)) != 0) {
          have_deletion = true;
        }
        if (chromosome_set.count(current) != 0) {
          circular = true;
        } 
        chromosome_set.insert(current); 
      }
    } else { 
      circular = true;
    } 
  } while ((current != Infty) && recovered_graphs[index].defined(current) && !circular);*/
  vertex_t current = graph.get_obverse_vertex(x); 
  vertex_t previous = x;
  chromosome_set.insert(x);
  
    for(;;) {
	if (chromosome_set.count(current) != 0) {
	    circular = true;
	    break; // circ
	}
	chromosome_set.insert(current);

    	changer_lambda(current, true);

	if (!recovered_graphs[index].defined(current)) { 
	  break; // linear
	} 

        previous = current;
	current = recovered_graphs[index][previous]; 

        /*if (bad_edges.count(std::make_pair(previous, current)) != 0 || bad_edges.count(std::make_pair(current, previous)) != 0) {
          have_deletion = true;
        }*/

	if (chromosome_set.count(current) != 0) {
	    circular = true;
	    break; // circ
	}
	chromosome_set.insert(current);

	if (current == Infty) { 
		break;
	}

        previous = current;
	current = graph.get_obverse_vertex(previous); 
        
        /*if (bad_edges.count(std::make_pair(previous, current)) != 0 || bad_edges.count(std::make_pair(current, previous)) != 0) {
          have_deletion = true;
        }*/
        
    }

    if (!circular && recovered_graphs[index].defined(x)) {	
	vertex_t y = x;
	while (recovered_graphs[index].defined(y) && (y != Infty)) {
	    y = recovered_graphs[index][y];
	    chromosome_set.insert(y);

	    if (y != Infty) {
	    	y = graph.get_obverse_vertex(y);
	    	chromosome_set.insert(y);
    	    	changer_lambda(y, false);
	    }
	}
    } 
#ifdef VERSION2
    else if (current == graph.get_obverse_vertex(previous) && circular) {
      return Chromosome();
    } 
#endif

    return Chromosome(path, circular);
}
#endif 

