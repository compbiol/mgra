#ifndef MULTIGRAPH_HPP
#define MULTIGRAPH_HPP

#include "defined.h" 

struct MultiGraph {
  typedef std::string orf_t;
  typedef structure::Genome genome_t;
  typedef utility::sym_multihashmap<vertex_t> partgraph_t;

  explicit MultiGraph(std::vector<genome_t> const & genomes) 
  : m_local_graphs(genomes.size()) 
  { 
    std::unordered_set<orf_t> blocks;

    for (auto const & genome : genomes) { 
      for(auto const & chromosome : genome) {
        for(auto const & orf : chromosome.second) {
      	  if (blocks.count(orf.second.first) == 0) { 
      	    obverse_edges.insert(std::make_pair(orf.second.first + "t", orf.second.first + "h"));
      	    obverse_edges.insert(std::make_pair(orf.second.first + "h", orf.second.first + "t"));
      	    blocks.insert(orf.second.first);
      	    vertex_set.insert(orf.second.first + "t"); 
      	    vertex_set.insert(orf.second.first + "h"); 
      	  }
        } 
      }
    }

    for(size_t i = 0; i < genomes.size(); ++i) {
      add_edges(i, genomes[i]); 
    }	
  } 

  inline void add_vertex(vertex_t const & v) { 
    vertex_set.insert(v);
  } 

  inline void erase_vertex(vertex_t const & v) { 
    vertex_set.erase(v);
  } 
    
  inline void add_edge(size_t index, vertex_t const & u, vertex_t const & v) { 
    assert(index < m_local_graphs.size());
    assert(u != Infty || v != Infty);
    m_local_graphs[index].insert(u, v);
  }

  inline void erase_edge(size_t index, vertex_t const & u, vertex_t const & v) {
    assert(index < m_local_graphs.size());
    assert(u != Infty || v != Infty);
    return m_local_graphs[index].erase(u, v);
  } 

  template<class mularcs_t>
  mularcs_t get_all_adjacent_multiedges(vertex_t const & u) const; 

  template<class mcolor_t>
  mcolor_t get_all_multicolor_edge(vertex_t const & u, vertex_t const & v) const;

  size_t degree_vertex(vertex_t const & u) const { 
    std::unordered_set<vertex_t> processed;
    for (size_t i = 0; i < m_local_graphs.size(); ++i) {
      auto iters = m_local_graphs[i].equal_range(u);
      for (auto it = iters.first; it != iters.second; ++it) { 
        processed.insert(it->second); 
      }
    }
    return processed.size();
  }

  inline vertex_t get_obverse_vertex(vertex_t const & v) const {
    assert(v != Infty);
    if (obverse_edges.count(v) != 0) {
      return obverse_edges.find(v)->second;
    } else { 
      return vertex_t();
    }
  }   

  inline bool is_identity() const { 
    for(auto it = m_local_graphs.cbegin(); it != m_local_graphs.cend() - 1; ++it) { 
      if (*it != *(it + 1)) {
        return false;
      }
    }
    return true;
  }

  inline partgraph_t const & get_partgraph(size_t index) const { 
    assert(index < m_local_graphs.size());
    return m_local_graphs[index];
  } 

  DECLARE_DELEGATE_CONST_METHOD( size_t, vertex_set, size, size )
  DECLARE_DELEGATE_CONST_METHOD( size_t, m_local_graphs, count_local_graphs, size )

  typedef std::set<vertex_t>::const_iterator citer;
  DECLARE_CONST_ITERATOR( citer, vertex_set, begin, cbegin )  
  DECLARE_CONST_ITERATOR( citer, vertex_set, end, cend )
  DECLARE_CONST_ITERATOR( citer, vertex_set, cbegin, cbegin )  
  DECLARE_CONST_ITERATOR( citer, vertex_set, cend, cend )
  
private:
  void add_edges(size_t index, genome_t const & genome) {
    auto const rear_lambda = [] (std::pair<orf_t, int> const & orf) -> vertex_t { 
      return ((orf.second > 0)?(orf.first + "h"):(orf.first + "t"));
    };

    auto const front_lambda = [] (std::pair<orf_t, int> const & orf) -> vertex_t { 
      return ((orf.second > 0)?(orf.first + "t"):(orf.first + "h"));
    };

    for (auto const & chromosome : genome) {
      vertex_t current_vertex = rear_lambda(chromosome.second.begin()->second);
      for (auto gene = (++chromosome.second.begin()); gene != chromosome.second.end(); ++gene) { 
        m_local_graphs[index].insert(current_vertex, front_lambda(gene->second));	
        current_vertex = rear_lambda(gene->second);
      }

      if (chromosome.second.is_circular()) {
      	m_local_graphs[index].insert(front_lambda(chromosome.second.begin()->second), rear_lambda((--chromosome.second.end())->second)); 
      } else { 
      	m_local_graphs[index].insert(front_lambda(chromosome.second.begin()->second), Infty);
      	m_local_graphs[index].insert(rear_lambda((--chromosome.second.end())->second), Infty); 
      }
    }
  }

protected:
  std::set<vertex_t> vertex_set; //set of vertice
  std::unordered_map<vertex_t, vertex_t> obverse_edges; //obverse relation 
  std::vector<partgraph_t> m_local_graphs; //local graphs of each color 
};	

template<class mularcs_t>
mularcs_t MultiGraph::get_all_adjacent_multiedges(vertex_t const & u) const {
  assert(u != Infty);
  
  mularcs_t output;

  for (size_t i = 0; i < m_local_graphs.size(); ++i) {
    auto iters = m_local_graphs[i].equal_range(u);
    for (auto it = iters.first; it != iters.second; ++it) { 
      output.insert(it->second, i); 
    }
  }

  return output;
} 

template<class mcolor_t>
mcolor_t MultiGraph::get_all_multicolor_edge(vertex_t const & u, vertex_t const & v) const {
  assert(u != Infty || v != Infty);

  mcolor_t result;

  for (size_t i = 0; i < m_local_graphs.size(); ++i) {
    auto iters = m_local_graphs[i].equal_range(u);
    for (auto it = iters.first; it != iters.second; ++it) { 
      if (it->second == v) { 
        result.insert(i);
      }
    } 
  } 

  return result;
}

  //FIXME IF WE RECONSTRUCT ANCESTORS WITH DUPLICATION EQUAL RANGE
  /*inline vertex_t get_adjecent_vertex(size_t index, vertex_t const & first) const {  
    assert(index < m_local_graphs.size() && (m_local_graphs[index].count(first) != 0));
    return m_local_graphs[index].find(first)->second;
  }    

  inline bool is_exist_edge(size_t index, vertex_t const & first) const { 
    assert(index < m_local_graphs.size());
    auto edge = m_local_graphs[index].find(first);
    if (edge != m_local_graphs[index].end() && edge->second != Infty)  {
      return true;
    }
    return false; 
  }*/

#endif
