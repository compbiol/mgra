#ifndef RECOVERED_INFO_HPP
#define RECOVERED_INFO_HPP

#include "Linearizator.hpp"

template<class graph_t>
struct RecoveredInfo {
  
  typedef typename graph_t::mcolor_type mcolor_t;
  
  typedef typename graph_t::partgraph_t partgraph_t;
  typedef typename graph_t::twobreak_t twobreak_t; 
  typedef typename graph_t::transform_t transform_t;

  typedef structure::Genome genome_t;
  typedef structure::Chromosome chromosome_t; 
  typedef std::tuple<transform_t, transform_t, transform_t> history_t;

  RecoveredInfo(graph_t const & graph);

  DECLARE_GETTER(std::vector<genome_t>, genomes, genomes);
  DECLARE_GETTER(std::vector<transform_t>, transformations, history);
   
private: 
  void get_ugly_history(std::map<mcolor_t, partgraph_t> & local_graphs, std::map<mcolor_t, transform_t> & transformations) const;

  genome_t get_genome(std::string const & name, partgraph_t const & recovered_graph);
  chromosome_t get_chromosome(partgraph_t const & recovered_graph, vertex_t const & x, std::unordered_set<vertex_t>& chromosome_set);

private: 
  graph_t const & m_graph;
  
  partgraph_t bad_edges;

  std::vector<genome_t> genomes;
  std::vector<transform_t> transformations;
}; 

template<class graph_t>
void RecoveredInfo<graph_t>::get_ugly_history(
    std::map<mcolor_t, partgraph_t> & local_graphs,
    std::map<mcolor_t, transform_t> & transformations) const { 

  for(auto im = m_graph.cbegin_vec_T_consistent_color(); im != m_graph.cend_vec_T_consistent_color(); ++im) {
    local_graphs.insert(std::make_pair(*im, m_graph.get_partgraph(0)));
    transformations.insert(std::make_pair(*im, transform_t())); 
  } 
  local_graphs.insert(std::make_pair(m_graph.get_root_color(), m_graph.get_partgraph(0)));
  transformations.insert(std::make_pair(m_graph.get_root_color(), transform_t())); 

  for(auto it = m_graph.crbegin_2break_history(); it != m_graph.crend_2break_history(); ++it) {
    for(auto im = m_graph.cbegin_vec_T_consistent_color(); im != m_graph.cend_vec_T_consistent_color(); ++im) {
      if (it->get_mcolor().includes(*im)) { 
        it->inverse().apply_single(local_graphs[*im]);
      }

      if (it->get_mcolor() == *im) {
        transformations[*im].push_back(it->inverse());
      }
    }
  }
}

template<class graph_t>
RecoveredInfo<graph_t>::RecoveredInfo(graph_t const & graph) 
: m_graph(graph)
, bad_edges(graph.get_bad_edges())
{ 
  assert(!cfg::get().is_target_build);

  /*Get transformation and graphs for linearization*/
  std::map<mcolor_t, partgraph_t> recovered_graphs;
  std::map<mcolor_t, transform_t> recovered_transformations;
  get_ugly_history(recovered_graphs, recovered_transformations);
  
  /*Algorithm for linearization*/
  Linearizator<graph_t> linearizator(graph);
    
  for (auto const & tree : cfg::get().phylotrees) {
    tree.walk_and_linearizeate(linearizator, recovered_graphs, recovered_transformations); 
  }

  /*Recover genomes and transformation*/
  for (auto const & local_graph : recovered_graphs) {
    genomes.push_back(get_genome(cfg::get().mcolor_to_name(local_graph.first), local_graph.second)); 
    transformations.push_back(recovered_transformations[local_graph.first]);
  }
}

template<class graph_t>
structure::Genome RecoveredInfo<graph_t>::get_genome(std::string const & name, partgraph_t const & recovered_graph) { 
  genome_t genome(name); 
  
  std::unordered_set<vertex_t> processed;
  std::string name_chr("chr");
  size_t count = 1; 

  for(vertex_t const & x : m_graph) { 
    if (processed.find(x) == processed.end()) { 
      std::unordered_set<vertex_t> chromosome_set;
      chromosome_t chromosome = get_chromosome(recovered_graph, x, chromosome_set);
      if (chromosome.size() != 0) {
        genome.insert(name_chr + std::to_string(count++), chromosome);
      } 
      std::copy(chromosome_set.cbegin(), chromosome_set.cend(), std::inserter(processed, processed.end()));
    }
  } 

  return genome;
}

template<class graph_t>
structure::Chromosome RecoveredInfo<graph_t>::get_chromosome(partgraph_t const & recovered_graph, vertex_t const & x, std::unordered_set<vertex_t>& chromosome_set) {
  std::list<std::pair<vertex_t, int> > path;
  bool circular = false;
  bool have_deletion = false;

  auto const changer_lambda = [&] (vertex_t const & t, bool flag) -> void {
    if (flag) { 
      (*t.rbegin() == 't')?path.push_back(std::make_pair(t.substr(0, t.size() - 1), -1)):path.push_back(std::make_pair(t.substr(0, t.size() - 1), 1));
    } else { 
      (*t.rbegin() == 't')?path.push_front(std::make_pair(t.substr(0, t.size() - 1), 1)):path.push_front(std::make_pair(t.substr(0, t.size() - 1), -1));  
    }  
  };

  vertex_t current = m_graph.get_obverse_vertex(x); 
  vertex_t previous = x;
  chromosome_set.insert(x);
  
  while ((current != Infty) && !circular) { 
    if (chromosome_set.count(current) != 0) {
      circular = true;
    } else { 
      chromosome_set.insert(current);
      changer_lambda(current, true);

      if (!recovered_graph.defined(current)) { 
        break; // linear
      } 

      previous = current;
      current = recovered_graph.find(previous)->second; 

      if (bad_edges.defined(previous, current)) {
        have_deletion = true;
      }

      if (chromosome_set.count(current) != 0) {
        circular = true;
      } else { 
        chromosome_set.insert(current);
        if (current != Infty) { 
          previous = current;
          current = m_graph.get_obverse_vertex(previous); 
        } 
      } 
    }
  }

  if (!circular && recovered_graph.defined(x)) {  
    vertex_t y = x;
    while (recovered_graph.defined(y) && (y != Infty)) {
      y = recovered_graph.find(y)->second;
      chromosome_set.insert(y);

      if (y != Infty) {
        y = m_graph.get_obverse_vertex(y);
        chromosome_set.insert(y);
        changer_lambda(y, false);
      }
    }
  } 
    
  if (have_deletion) {
    return chromosome_t();
  }

  return chromosome_t(path, circular);
}

#endif