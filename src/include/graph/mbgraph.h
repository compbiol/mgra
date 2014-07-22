#ifndef MBGRAPH_H_
#define MBGRAPH_H_

#include "defined.h" 

struct MBGraph {
  typedef structure::Genome genome_t;
  typedef std::string orf_t;

  explicit MBGraph(std::vector<genome_t> const & genomes) 
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

  inline vertex_t get_obverse_vertex(vertex_t const & v) const {
    //assert(obverse_edges.count(v) != 0);
    if (obverse_edges.count(v) != 0) {
      return obverse_edges.find(v)->second;
    } else { 
      return vertex_t();
    }
  } 

  //FIXME IF WE RECONSTRUCT ANCESTORS WITH DUPLICATION EQUAL RANGE
  inline vertex_t get_adjecent_vertex(size_t index, vertex_t const & first) const {  
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
  }
 
  inline bool is_identity() { 
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

  inline size_t size() const { 
    return vertex_set.size(); 
  }
 
  inline size_t count_local_graphs() const { 
    return m_local_graphs.size();
  } 

  typedef std::set<vertex_t>::const_iterator citer;
  DECLARE_CONST_ITERATOR( citer, vertex_set, begin, cbegin )  
  DECLARE_CONST_ITERATOR( citer, vertex_set, end, cend )
  DECLARE_CONST_ITERATOR( citer, vertex_set, cbegin, cbegin )  
  DECLARE_CONST_ITERATOR( citer, vertex_set, cend, cend )
  
protected: 
  inline void add_edge(size_t index, vertex_t const & first, vertex_t const & second) { 
    assert(index < m_local_graphs.size());
    assert(first != Infty || second != Infty);
    m_local_graphs[index].insert(first, second);
  }

  inline void erase_edge(size_t index, vertex_t const & first, vertex_t const & second) {
    assert(index < m_local_graphs.size());
    assert(first != Infty || second != Infty);
    return m_local_graphs[index].erase(first, second);
  } 

private:
  void add_edges(size_t index, genome_t const & genome) {
    auto const rearLambda = [] (std::pair<orf_t, int> const & orf) -> vertex_t { 
      return ((orf.second > 0)?(orf.first + "h"):(orf.first + "t"));
    };

    auto const frontLambda = [] (std::pair<orf_t, int> const & orf) -> vertex_t { 
      return ((orf.second > 0)?(orf.first + "t"):(orf.first + "h"));
    };

    for (auto const & chromosome : genome) {
      vertex_t current_vertex = rearLambda(chromosome.second.begin()->second);
      for (auto gene = (++chromosome.second.begin()); gene != chromosome.second.end(); ++gene) { 
        m_local_graphs[index].insert(current_vertex, frontLambda(gene->second));	
        current_vertex = rearLambda(gene->second);
      }

      if (chromosome.second.is_circular()) {
      	m_local_graphs[index].insert(frontLambda(chromosome.second.begin()->second), rearLambda((--chromosome.second.end())->second)); 
      } else { 
      	m_local_graphs[index].insert(frontLambda(chromosome.second.begin()->second), Infty);
      	m_local_graphs[index].insert(rearLambda((--chromosome.second.end())->second), Infty); 
      }
    }
  }

protected:
  std::set<vertex_t> vertex_set; //set of vertice
  std::unordered_map<vertex_t, vertex_t> obverse_edges; //obverse relation 
  std::vector<partgraph_t> m_local_graphs; //local graphs of each color 
};	

#endif
