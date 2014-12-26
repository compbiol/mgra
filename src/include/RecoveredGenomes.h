#ifndef RECOVEREDGENOMES_H_
#define RECOVEREDGENOMES_H_

#include "Decircularizeter.h"

template<class graph_t>
struct RecoveredGenomes {
  
  typedef typename graph_t::mcolor_type mcolor_t;
  
  typedef typename graph_t::partgraph_t partgraph_t;
  typedef typename graph_t::twobreak_t twobreak_t; 
  typedef typename graph_t::transform_t transform_t;
 
  typedef structure::Genome genome_t;
  typedef structure::Chromosome chromosome_t;

  RecoveredGenomes(graph_t const & gr, partgraph_t const & b_edges);

  inline std::vector<genome_t> get_genomes() { 
    std::vector<genome_t> genomes;
    for (size_t i = 0; i < recovered_graphs.size(); ++i) {
      genomes.push_back(get_genome(i)); 
    } 
    return genomes;
  } 

  std::vector<transform_t> get_history() const { 
    return recovered_transformation;
  }

private: 
  bool are_adjacent_branches(mcolor_t const & A, mcolor_t const & B) const; 

  genome_t get_genome(size_t index);
  chromosome_t get_chromosome(size_t index, vertex_t const & x, std::unordered_set<vertex_t>& chromosome_set);

private: 
  graph_t const & graph;
  partgraph_t bad_edges;
  std::vector<std::string> name_genomes;
  std::vector<partgraph_t> recovered_graphs;
  std::vector<transform_t> recovered_transformation;
};

template<class graph_t>
bool RecoveredGenomes<graph_t>::are_adjacent_branches(mcolor_t const & mcolor_a, mcolor_t const & mcolor_b) const { 
  if (!graph.is_T_consistent_color(mcolor_a) || !graph.is_T_consistent_color(mcolor_b)) { 
    return false;
  } 

  mcolor_t color1; 
  mcolor_t color2;
  
  if (mcolor_a.size() >= mcolor_b.size()) {
    color1 = mcolor_a;
    color2 = mcolor_b;
  } else {
    color1 = mcolor_b;
    color2 = mcolor_a;
  }
  
  mcolor_t diff_color(color1, color2, mcolor_t::Difference);
  
  if (diff_color.size() == color1.size() - color2.size() && graph.is_T_consistent_color(diff_color)) {    
    return true;
  } 
  
  mcolor_t union_color(color1, color2, mcolor_t::Union); 
  if (union_color.size() == color1.size() + color2.size() && graph.is_T_consistent_color(union_color)) {  
    return true;
  } 
  
  return false;
}

template<class graph_t>
RecoveredGenomes<graph_t>::RecoveredGenomes(graph_t const & gr, partgraph_t const & b_edges)  
: graph(gr)
, bad_edges(b_edges)
{
  if (cfg::get().is_target_build) { 
    ;
    /*mcolor_t const & target = conf.get_target();
    name_genomes.push_back(conf.mcolor_to_name(target));
    recovered_graphs.resize(1);
    for(auto const &x : graph) {
      if (recovered_graphs[0].defined(x)) { 
        continue;
      }

      vertex_t y = Infty;
      bool good = true;
      size_t def = 0;
      std::for_each(target.cbegin(), target.cend(), [&] (std::pair<size_t, size_t> const & col) -> void {
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
      }); 
     
      if (good && def == target.size() && y != Infty) {
        recovered_graphs[0].insert(x, y);
      }
    }*/
  } else {
    recovered_transformation.resize(graph.count_vec_T_consitent_color());
    recovered_graphs.resize(graph.count_vec_T_consitent_color() + 1, graph.get_partgraph(0));

    for(auto it = graph.crbegin_2break_history(); it != graph.crend_2break_history(); ++it) {
      size_t i = 0;
      for(auto im = graph.cbegin_vec_T_consistent_color(); im != graph.cend_vec_T_consistent_color(); ++im, ++i) {
        if (it->get_mcolor().includes(*im)) { 
          //std::cerr << conf.mcolor_to_name(*im) << std::endl; 
          //std::cerr << it->get_vertex(0) << " " << it->get_vertex(1) << " " << it->get_vertex(2) << " " << it->get_vertex(3) << std::endl;
          it->inverse().apply_single(recovered_graphs[i]);
        }

        if (it->get_mcolor() == *im) {
          recovered_transformation[i].push_front(*it);
        }
      }
    }

    size_t i = 0; 
    Decircularizeter<graph_t> dec(graph, bad_edges);          
    for (auto im = graph.cbegin_vec_T_consistent_color(); im != graph.cend_vec_T_consistent_color(); ++im, ++i) {    
      name_genomes.push_back(cfg::get().mcolor_to_name(*im));
      //std::cerr << "Start decircularize procedure for " << conf.mcolor_to_name(*im) << std::endl;
      //std::cerr << "Initial we have " << dec.count_circular_chromosome(recovered_graphs[i]) << std::endl;
      transform_t T = dec.decircularize(recovered_graphs[i], recovered_transformation[i]);
      //std::cerr << "After we have " << dec.count_circular_chromosome(recovered_graphs[i]) << std::endl;

      // move to adjacent branches
      for (auto const &it : T) {
        size_t j = 0; 
        for (auto jt = graph.cbegin_vec_T_consistent_color(); jt != graph.cend_vec_T_consistent_color(); ++jt, ++j) {
          if ((j != i) && includes(im->cbegin(), im->cend(), jt->cbegin(), jt->cend()) && are_adjacent_branches(*im, *jt)) {
            recovered_transformation[j].push_back(it);
          }
        }
      }
    }
    name_genomes.push_back(cfg::get().mcolor_to_name(graph.get_root_color()));
  }
}

template<class graph_t>
structure::Genome RecoveredGenomes<graph_t>::get_genome(size_t index) { 
  genome_t genome(name_genomes[index]); 
  std::unordered_set<vertex_t> processed;
  std::string name_chr("chr");
  size_t count = 1; 

  for(auto const & x : graph) { 
    if (processed.find(x) == processed.end()) { 
      std::unordered_set<vertex_t> chromosome_set;
      chromosome_t chromosome = get_chromosome(index, x, chromosome_set);
      if (chromosome.size() != 0) {
        genome.insert(name_chr + std::to_string(count++), chromosome);
      } 
      std::copy(chromosome_set.cbegin(), chromosome_set.cend(), std::inserter(processed, processed.end()));
    }
  } 

  return genome;
} 

template<class graph_t>
structure::Chromosome RecoveredGenomes<graph_t>::get_chromosome(size_t index, vertex_t const & x, std::unordered_set<vertex_t>& chromosome_set) {
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

  vertex_t current = graph.get_obverse_vertex(x); 
  vertex_t previous = x;
  chromosome_set.insert(x);
  
  while ((current != Infty) && !circular) { 
    if (chromosome_set.count(current) != 0) {
      circular = true;
    } else { 
      chromosome_set.insert(current);
      changer_lambda(current, true);

      if (!recovered_graphs[index].defined(current)) { 
        break; // linear
      } 

      previous = current;
      current = recovered_graphs[index][previous]; 

      if (bad_edges.defined(previous, current)) {
        have_deletion = true;
      }

      if (chromosome_set.count(current) != 0) {
        circular = true;
      } else { 
        chromosome_set.insert(current);
        if (current != Infty) { 
          previous = current;
          current = graph.get_obverse_vertex(previous); 
        } 
      } 
    }
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
    
  if (have_deletion) {
    return chromosome_t();
  }

  return chromosome_t(path, circular);
}
#endif 

