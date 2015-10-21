#ifndef LINEARIZE_UTILS_HPP
#define LINEARIZE_UTILS_HPP

namespace algo { 

namespace linearize { 

/**
 *
 */
template<class partgraph_t, class transformation_t>
partgraph_t apply_transformation(partgraph_t traverse_graph, transformation_t const & transform) {
  for (auto const & twobreak : transform) { 
    twobreak.apply_single(traverse_graph);
  } 
  return traverse_graph;
} 

/**
 * Function get one chromosome from LOCAL_GRAPH, where X - is vertex which belong our chromosome.
 */
template<class graph_pack_t, class partgraph_t>
structure::Chromosome get_chromosome(graph_pack_t const & graph_pack, partgraph_t const & local_graph, vertex_t const & x, std::unordered_set<vertex_t>& processed) {
  std::list<std::pair<vertex_t, int> > path;
  bool circular = false;
  bool have_deletion = false;
  vertex_t previous = x;
  vertex_t current = x; 
  processed.insert(x);

  auto const changer_lambda = [&] (vertex_t const & t, bool flag) -> void {
    if (flag) { 
      (*t.rbegin() == 't')?path.push_back(std::make_pair(t.substr(0, t.size() - 1), -1)):path.push_back(std::make_pair(t.substr(0, t.size() - 1), 1));
    } else { 
      (*t.rbegin() == 't')?path.push_front(std::make_pair(t.substr(0, t.size() - 1), 1)):path.push_front(std::make_pair(t.substr(0, t.size() - 1), -1));  
    }  
  };

  do { 
    previous = current; 
    current = graph_pack.graph.get_obverse_vertex(previous); 
  
    if (processed.count(current) == 0) { 
      processed.insert(current);
      changer_lambda(current, true);

      if (local_graph.defined(current)) { 
        previous = current; 
        current = local_graph.find(previous)->second; 

        //if (graph_pack.is_prosthetic_chromosome(previous, current)) {
        if (current == graph_pack.graph.get_obverse_vertex(previous)) {
          have_deletion = true;
        }

        if (processed.count(current) != 0) {
          circular = true;
        } else if (current != Infty) { 
          processed.insert(current);  
        }
      } 
    } else { 
      circular = true;
    } 
  } while ((current != Infty) && local_graph.defined(current) && !circular);

  if (!circular && local_graph.defined(x)) {  
    for (vertex_t y = local_graph.find(x)->second; local_graph.defined(y) && (y != Infty); y = local_graph.find(y)->second) {
      processed.insert(y);
      if (y != Infty) {
        y = graph_pack.graph.get_obverse_vertex(y);
        processed.insert(y);
        changer_lambda(y, false);
      }
    } 
  } 
    
  if (have_deletion) {
    return structure::Chromosome();
  }

  return structure::Chromosome(path, circular);
}

/**
 * create genome from local graph with specified color and name. 
 * return genome in structure::Genome class. 
 */
template<class graph_pack_t, class partgraph_t>
structure::Genome get_genome(graph_pack_t const & graph_pack, std::string const & name, partgraph_t const & local_graph) {
	std::string const name_chr("chr");
  std::unordered_set<vertex_t> processed;
  structure::Genome genome(name); 
  size_t count = 1; 

  for (vertex_t const & x : graph_pack.graph) { 
    if (processed.find(x) == processed.end()) { 
      structure::Chromosome chromosome = algo::linearize::get_chromosome(graph_pack, local_graph, x, processed);
      if (chromosome.size() != 0) {
        genome.insert(name_chr + std::to_string(count++), chromosome);
      } 
    }
  } 

  return genome;
}

template<class graph_pack_t, class partgraph_t>
structure::Chromosome get_chromosome_by_edge(graph_pack_t const & graph_pack, partgraph_t const & local_graph, 
                typename graph_pack_t::edge_t const & e, std::unordered_set<vertex_t>& processed) {
  
  if (e.first != Infty) { 
    return get_chromosome(graph_pack, local_graph, e.first, processed);
  } else if (e.second != Infty) { 
    return get_chromosome(graph_pack, local_graph, e.second, processed);
  } 
  return structure::Chromosome();
} 


/**
 * count number of circular chromosomes from genome by graph with specified color. 
 * in paper about linearization algorithm this functon corresponding c(.)
 * return number of cicular chromosomes in local graph (genome). 
 */
template<class graph_pack_t, class partgraph_t>
size_t count_circular_chromosome(graph_pack_t const & graph_pack, partgraph_t const & local_graph) { 
	std::unordered_set<vertex_t> processed;
  size_t count_circular_chr = 0;
  
  for (vertex_t const & x : graph_pack.graph) { 
    if (processed.count(x) == 0) { 
      structure::Chromosome chromosome = get_chromosome(graph_pack, local_graph, x, processed); 
      if (chromosome.is_circular()) {
        ++count_circular_chr; 
      }
    } 
  }
  
  return count_circular_chr;
}

}

}

#endif 